function data_reg = ...
    sbxDemonsRegOneshotCore(mov_path, savepath, startframe, nframes, ref,...
    meshxx, meshyy, varargin)
%% sbxDemonsRegOneshotCore is the core of the namesake function, aiming to improve speed
    
    p = inputParser;
    addOptional(p, 'pmt', 0, @isnumeric);  % REMEMBER, PMT is 0-indexed   
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
    parse(p, varargin{:});
    
    % Unpack if needed
    if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
        varargin = varargin{:};
    end
    
    parse(p, varargin{:});
    p = p.Results;

%% Filter ref
% Parameters
n = p.hp_norm_sigmas(1);
m = p.hp_norm_sigmas(2);

% Median filter
if ~isempty(p.medfilt2size)
    ref = medfilt2(ref, p.medfilt2size, 'symmetric');
end

if p.highpassnorm
    % Highpass and normalized
    ref_prime = single(ref)-single(imgaussfilt(double(ref),n));
    ref = ref_prime ./ (imgaussfilt(ref_prime.^2,m) .^ (1/2));

    ref(isnan(ref)) = 0;
end

%% Read in data
% Read in
data = sbxReadPMT(mov_path, startframe - 1, nframes, p.pmt);

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
        f_prime = data2(:,:,i) - imgaussfilt(single(data2(:,:,i)),n);
        g_prime = f_prime ./ (imgaussfilt(f_prime.^2,m).^(1/2));

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
    [D,~] = imregdemons(data2(:,:,i), ref, [32 16 8 4],...
            'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitbar',false);
    
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
if ~strcmp(class(data_reg), class(data))
    data_reg = cast(data_reg, c);
end

%% Save warp if needed
if p.savewarp
    writetiff(D_combined,[savepath '.tif'],'single'); 
end
end