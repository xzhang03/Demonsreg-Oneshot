function data_reg = sbxDemonsRegOneshotTiffCore(data, savepath, ref,...
    meshxx, meshyy, varargin)
%sbxDemonsRegOneshotTiffCore is the core of the namesake function, aiming to improve speed
% This is the tiff version of the code. - SZ

%% Parsing
    
p = inputParser;
addOptional(p, 'ref_downsample_xy', 1, @isnumeric);
addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric); % Sigma for gaussian fit
addOptional(p, 'savewarp', true);
addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter
addOptional(p, 'highpassnorm', true); % Use highpass filter and local normalization before registration.
addOptional(p, 'binbeforehighpassnorm', false); % Bin image before applying spatial filters.
                                                % Recommended for
                                                % imaging systems with
                                                % magnification (e.g.,
                                                % GRIN doublets).
addOptional(p, 'edges', [0 0 0 0]); % Use edges to the processing (needed for bidirectional scanning).

% imdemonreg variables
addOptional(p, 'itr', [32 16 8 4]); % Iterations at each level
addOptional(p, 'PyramidLevels', 4); % Number of levels
addOptional(p, 'AccumulatedFieldSmoothing', 2.5); % Gaussian size for smoothing

parse(p, varargin{:});

% Unpack if needed
if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
    varargin = varargin{:};
end

parse(p, varargin{:});
p = p.Results;


%% Read in data
% Class of data
c = class(data);

%% Make a copy for pre-processing while keeping the original
% For some reason uint16 gives crappy results
% Apply edges if needed.
data2 = single(data(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :));

%% Filter and bin data
% Bin if necessary
if p.binbeforehighpassnorm
    if p.ref_downsample_xy > 1
        data2 = binxy(data2, p.ref_downsample_xy);
    end
end

if p.highpassnorm
    for  i = 1:size(data2,3)
        % Prepare the data by median filtering and local normalizing
        if ~isempty(p.medfilt2size)
            data2(:,:,i) = medfilt2(data2(:,:,i), p.medfilt2size, 'symmetric');
        end
        f_prime = data2(:,:,i) - imgaussfilt(data2(:,:,i),p.hp_norm_sigmas(1));
        g_prime = f_prime ./ (imgaussfilt(f_prime.^2,p.hp_norm_sigmas(2)).^(1/2));

        g_prime(isnan(g_prime)) = 0;

        data2(:,:,i)=g_prime;

    end
end

% Bin, if necessary
if ~p.binbeforehighpassnorm
    if p.ref_downsample_xy > 1
        data2 = binxy(data2, p.ref_downsample_xy);
    end
end

%% Regsitration
% Initialize D_combined
D_combined = zeros(size(data2,1), size(data2,2), size(data2,3) * 2);

% Get the D tensor
for  i = 1:size(data2,3)
    % Demons reg
    [D,~] = imregdemons(data2(:,:,i), ref, p.itr,...
            'AccumulatedFieldSmoothing', p.AccumulatedFieldSmoothing,...
            'PyramidLevels',p.PyramidLevels,'DisplayWaitbar',false);
    
    D_combined(:, :, i*2-1 : i*2) = D; 
end

% resize back if necessary
if p.ref_downsample_xy > 1
    D_combined = imresize(p.ref_downsample_xy * D_combined, p.ref_downsample_xy); % resize to bring back to full size of movie
end

% Re-embed D-combined into the full-size version if using edges
if any(p.edges ~= 0)
    D_combined_full = zeros(size(data,1), size(data,2), size(data,3) * 2);
    D_combined_full(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2), :) = D_combined;
    
    % Use the new D_combined for subsequent processing
    D_combined = D_combined_full;
    clear D_combined_full
end

% Data2 is a large variable
if p.ref_downsample_xy > 1
    clear data2;
end

%% Apply the warp to original movie
% Initialize
data_reg = zeros(size(data));

% Apply
for i = 1 : size(data, 3)
    data_reg(:,:,i) = ...
        interp2(meshxx, meshyy, single(data(:,:,i)),...
        meshxx + D_combined(:,:,2*i-1), meshyy + D_combined(:,:,2*i));
end

% Cast if needed
if ~strcmp(class(data_reg), c)
    data_reg = cast(data_reg, c);
end

%% Save warp if needed
if p.savewarp
    writetiff(D_combined(:,:,1:2:end),[savepath '_X.tif'],'single'); 
    writetiff(D_combined(:,:,2:2:end),[savepath '_Y.tif'],'single'); 
end
end