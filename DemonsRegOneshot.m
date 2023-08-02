function [data_reg, toseg] = DemonsRegOneshot(im, varargin)
% sbxDemonsRegOneshot is a onestep application of demonsreg, with the
% intention of improving speed. This is the RAM version of the code where
% you pass the movie as a variable.

    p = inputParser;

    % Lab variables variables
    addOptional(p, 'ref', []);  % Ref target
    
    % Image-processing variables
    addOptional(p, 'edges', [0 0 0 0]); % Edges of registration. Pixels outside the edges will not be registered
    addOptional(p, 'refsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
    addOptional(p, 'refoffset', 500, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames
    addOptional(p, 'binxy', 1, @isnumeric);  % Pixels to downsample in xy, will be set to downsample_xy if empty
    
    % Rarely changed image-processing variables
    addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter. Leave empty if no median filter
    addOptional(p, 'highpassnorm', true); % Use highpass filter and local normalization before registration.
    addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric);  % sigmas for highpass filter and local normalization
    addOptional(p, 'binbeforehighpassnorm', true); % Bin image before applying spatial filters. Recommended for imaging
                                                    % systems with magenification (e.g., GRIN doublets)                                                    % magnification
  
    % Efficiency variables
    addOptional(p, 'iterations', 1); % Enable multi-round registration
    addOptional(p, 'parfor', true);
    addOptional(p, 'nworkers', 32); % Max number of workers
    addOptional(p, 'chunksize', 1000); % Chunk size for parallel processing. Decrease if RAM is an issue
        
    % Post processing variables
    addOptional(p, 'posthocmedian', false); % Perform a sliding-window median after registration
    addOptional(p, 'posthocmedianwindow', 10); % Number of frames as the window for posthoc median filter
    
    % imdemonreg variables
    addOptional(p, 'itr', [32 16 8 4]); % Iterations at each level
    addOptional(p, 'PyramidLevels', 4); % Number of levels
    addOptional(p, 'AccumulatedFieldSmoothing', 2.5); % Gaussian size for smoothing
    
    % Export segmentation file
    addOptional(p, 'toseg', false);
    addOptional(p, 'movsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
    addOptional(p, 'movoffset', 500, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames

    % Unpack if needed
    if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
        varargin = varargin{:};
    end

    parse(p, varargin{:});
    p = p.Results;

%% Input cleanup
% Parpool
if p.parfor
    if isempty(gcp('nocreate'))
        parpool();
    end
end

% sizes
sz = size(im);
if length(sz) == 2
    sz(3) = 1;
end

%% Registration
for iteration = 1 : p.iterations    
    fprintf('Iteration: %i\n', iteration);

    % get ref
    if ~isempty(p.ref)
        fprintf('Using external reference');
        
        if p.iterations > 1
            % Set iteration to 1
            fprintf(', recommend changing iterations to 1 in the future\n')
        else
            fprintf('\n')
        end

        % Check ref is compatible
        ref = single(p.ref);
        refsize = size(ref);
        if (refsize(1) ~= sz(1)) || (refsize(2) ~= sz(2))
            disp('Ref size does not match');
            return;
        end

    else 
        tic;
        
        % Ref cauclation
        ref = single(median(im(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2),...
            p.refoffset:p.refoffset+p.refsize-1),3));
        
        t = toc;
        fprintf('Ref calculation done. Elapsed time: %i seconds.\n', round(t));
    end

    % Show ref
    figure
    imshow(ref, []);

    % Bin ref (pre)
    if p.binbeforehighpassnorm && p.binxy > 1
        ref = binxy(ref, p.binxy);
    end
    
    % Local normalize ref
    ref = medfilt2(ref, p.medfilt2size, 'symmetric');
    f_prime = ref - imgaussfilt(ref, p.hp_norm_sigmas(1));
    ref = f_prime ./ (imgaussfilt(f_prime.^2, p.hp_norm_sigmas(2)).^(1/2));
    ref(isnan(ref)) = 0;
    
    % Bin ref (post)
    if ~p.binbeforehighpassnorm && p.binxy > 1
        ref = binxy(ref, p.binxy);
    end
    
    % Parallel processing
    nchunks = ceil(sz(3)/p.chunksize);
    
    % Get the meshgrid to pass to all the workers
    [xx, yy] = meshgrid(1 : sz(2), 1 : sz(1));
    
    % Initialize registered data (using the grid file xx to get size)
    fprintf('Initializing array...')
    tic;
    data_reg = uint16(zeros(sz(1), sz(2)));
    data_reg = repmat(data_reg, [1, 1, sz(3)]);
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    % Initialize a cell for registered data
    data_reg_cell = cell(nchunks, 1);
    
    % Slice data
    fprintf('Slicing data...')
    tic;
    for c = 1 : nchunks
        data_reg_cell{c} = im(:,:,(c-1) * p.chunksize + 1 :...
            min(c * p.chunksize, sz(3)));
    end
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    fprintf('Parallel registration...')
    tic;
    
    % Parallel or not
    if p.parfor
        M = p.nworkers;
    else
        M = 0;
    end
    
    % Parallel processing
    parfor (c = 1 : nchunks, M)
        % Process data
        data_reg_cell{c} = ...
            uint16(sbxDemonsRegOneshotTiffCore(data_reg_cell{c}, '', ref, xx,...
            yy, 'ref_downsample_xy', p.binxy,...
            'hp_norm_sigmas', p.hp_norm_sigmas, 'savewarp', false, ...
            'medfilt2size', p.medfilt2size, 'binbeforehighpassnorm', p.binbeforehighpassnorm,...
            'highpassnorm', p.highpassnorm, 'edges', p.edges, 'itr', p.itr, 'PyramidLevels', p.PyramidLevels,...
            'AccumulatedFieldSmoothing', p.AccumulatedFieldSmoothing));
        
    end
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    % Reconstruct image stack
    fprintf('Reconstruct image stack...')
    tic;
    for c = 1 : nchunks
        % time indices
        i_start = (c-1) * p.chunksize + 1;
        i_end = min(c * p.chunksize, sz(3));
        
        % Reconstruct
        data_reg(:,:,i_start : i_end) = data_reg_cell{c};
        
    end

    % Free memory
    clear data_reg_cell
    t = toc;

    % Prepare for next iteration
    if iteration < p.iterations
        im = data_reg;
    end
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));

end

%% Post hoc median filter
if p.posthocmedian
    fprintf('Posthoc median filter...')
    tic;
    
    % filter
    data_reg = movmedian(data_reg, p.posthocmedianwindow, 3);
    
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
end

%% Segmentation file
if p.toseg
    % mean image
    meanim = median(data_reg(:,:,p.movoffset : p.movoffset + p.movsize - 1),3);
    meanim = double(meanim);
    meanim = imresize(meanim,1/p.binxy);
    
    meanim = medfilt2(meanim, p.medfilt2size, 'symmetric');
    meanim_prime = meanim - imgaussfilt(double(meanim),p.hp_norm_sigmas(1));
    ln_meanim = meanim_prime ./ (imgaussfilt(meanim_prime.^2,p.hp_norm_sigmas(2)) .^ (1/2));
    ln_meanim(isnan(ln_meanim)) = 0;

    toseg = ln_meanim;
else
    toseg = [];
end

end