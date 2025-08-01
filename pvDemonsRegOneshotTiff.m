function pvDemonsRegOneshotTiff(mouse, date, varargin)
% pvDemonsRegOneshot is a onestep application of demonsreg, with the
% intention of improving speed. This is the tiff version of the code.
% -SZ

    p = inputParser;

    % Lab variables variables
    addOptional(p, 'target', []);  % Target run to align to, default runs(1)
    addOptional(p, 'server', '');  % Add in the server name as a string
    addOptional(p, 'runs', []);  % runs of the movie. Can be multiple
    addOptional(p, 'force', false);  % Overwrite files if they exist
    addOptional(p, 'pmt', [], @isnumeric);  % Which PMT to use for analysis, 0-green, 1-red
    addOptional(p, 'chunksize', 300); % Chunk size for parallel processing. Decrease if RAM is an issue
    
    % Plane level
    addOptional(p, 'planes', []); % Plane level
    
    % IO variables
    addOptional(p, 'movtype', 'xyreg');  % input type, can be xyreg, sbxreg, or sbx.
    addOptional(p, 'iterations', 1); % Multi-iteration registration
    
    % Image-processing variables
    addOptional(p, 'edges', [0 0 0 0]); % Edges of registration. Pixels outside the edges will not be registered
    addOptional(p, 'refsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
    addOptional(p, 'refoffset', 100, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames
    addOptional(p, 'binxy', [], @isnumeric);  % Pixels to downsample in xy, will be set to downsample_xy if empty
    addOptional(p, 'reuseref', true); % Reuse reference if can find it (must be the same size).
    addOptional(p, 'showref', true);

    % Rarely changed image-processing variables
    addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter. Leave empty if no median filter
    addOptional(p, 'highpassnorm', true); % Use highpass filter and local normalization before registration.
    addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric);  % sigmas for highpass filter and local normalization
    addOptional(p, 'binbeforehighpassnorm', false); % Bin image before applying spatial filters. Recommended for imaging
                                                    % systems with magenification (e.g., GRIN doublets)                                                    % magnification
    addOptional(p, 'previewlocalnorm', false); % Preview local normalization results

    % Efficiency variables
    addOptional(p, 'savewarp', false); % Save warp parameteres (only useful for applying the same parameter to a different channge)
    addOptional(p, 'parfor', true);
    addOptional(p, 'nworkers', 32); % Max number of workers
    
    % Post processing variables
    addOptional(p, 'posthocmedian', false); % Perform a sliding-window median after registration
    addOptional(p, 'posthocmedianwindow', 10); % Number of frames as the window for posthoc median filter
    
    % PCA cleaning
    addOptional(p, 'PCAclean', false);
    addOptional(p, 'PCs_clean', 2:30);

    % imdemonreg variables
    addOptional(p, 'itr', [32 16 8 4]); % Iterations at each level
    addOptional(p, 'PyramidLevels', 4); % Number of levels
    addOptional(p, 'AccumulatedFieldSmoothing', 2.5); % Gaussian size for smoothing
    
    % Export segmentation file
    addOptional(p, 'toseg', false);
    addOptional(p, 'movsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
    addOptional(p, 'movoffset', 500, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames
    
    % Expoert to preview
    addOptional(p, 'topreview', false);
    addOptional(p, 'previewbin_xyt', [2 4]); %[xy t]

    % Unpack if needed
    if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
        varargin = varargin{:};
    end

    parse(p, varargin{:});
    p = p.Results;

%% Input cleanup
if isempty(p.target), p.target = p.runs(1); end
    
% Fix pmts
if isempty(p.pmt)
    p.pmt = pvPMT(mouse, date, p.runs, 'server', p.server);
    if length(p.pmt) > 1
        error('More than 1 pmt detected, please specify');
    end
end

% Parpool
if isempty(gcp('nocreate')) && p.parfor
    parpool();
end

%% IO paths
% sbx paths
tiffpaths = cell(length(p.runs),1);
for i = 1 : length(p.runs)
    tiffpaths{i} = pvPaths(mouse, date, p.runs(i), p.movtype, 'server', p.server, 'pmt', p.pmt, 'planes', p.planes);
end

% ref paths
istarget = p.runs == p.target;

% Rearrange so that the target gets processed first
p.runs = [p.runs(istarget), p.runs(~istarget)];
[rpath, fname] = fileparts(tiffpaths{p.runs == p.target});
fname = pvParsefn(fname);

regpath = fullfile(fileparts(rpath), 'demonsreg');
refname = fullfile(regpath, sprintf('%s_demonsregref.tif', fname));
matname = fullfile(regpath, sprintf('%s_demonsreg.mat', fname));

% Make demonsreg folder if necessary
if ~exist(regpath, 'dir')
    mkdir(regpath);
end 

% output paths
outpaths = cell(length(p.runs),1);
for i = 1 : length(p.runs)
    switch p.movtype
        otherwise
            outpaths{i} = sprintf('%s_demonsreg.tif', tiffpaths{i}(1:end-4));
    end
end


%% Target
% get ref
fprintf('======= Demonsreg ========\n')
if exist(refname, 'file') && p.reuseref
    % Re-use a ref file
    fprintf('Using existing reference image...\n')
    ref = readtiff(refname);
elseif ~any(istarget)
    % Calculate ref fresh from a differen run (why man)
    fprintf('Calculating reference image from a different run...\n')
    tic;
    
    % Read tiff
    refstack = readtiff(fullfile(rpath, fname));
    ref = single(median(refstack(p.edges(3)+1:end-p.edges(4), p.edges(1)+1:end-p.edges(2),...
        p.refoffset:p.refoffset+p.refsize-1),3));
    
    clear refstack
    
    t = toc;
    fprintf('Ref calculation done. %i seconds.\n', round(t));
    writetiff(ref, refname); % Write unbinned because sometimes binning happens after local normalization
end

%% Registration
for i = 1:length(tiffpaths)
    % Get the path and info file
    movpath = tiffpaths{i};

    if exist(outpaths{i}, 'file') && ~p.force
        fprintf('%s\n', outpaths{i});
        fprintf('Demonsreg files already exist, skipping. \n');
        continue;
    end
    
    for iter = 1 : p.iterations
        if ~p.parfor
            fprintf('[Iteration %i]\n', iter);
        end

        if iter == 1
            % Read
            im = readtiff(movpath);

            % sizes
            sz = size(im);
            
        else
            im = data_reg;
        end

        % Make ref
        if (~exist('ref', 'var') && (i == 1)) || (iter > 1)
            ref = single(median(im(p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),...
                p.refoffset:p.refoffset+p.refsize-1),3));

            % Write ref
            writetiff(ref, refname); % Write unbinned because sometimes binning happens after local normalization
        end

        % Bin ref (pre)
        if p.binbeforehighpassnorm && p.binxy > 1
            ref = binxy(ref, p.binxy);
        end
        
        % Preview local normalization
        if iter == 1 && p.previewlocalnorm
            lnparas = localnormalize_ui('im', ref, 'gausssizes', p.hp_norm_sigmas, 'parametermode', true);
            p.hp_norm_sigmas = lnparas([2,1]);
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
        
        if p.showref
            figure
            imshow(ref, []);
            title(sprintf('Iteration %i', iter));
        end

        % Parallel processing
        nchunks = ceil(sz(3)/p.chunksize);

        % Get the meshgrid to pass to all the workers
        [xx, yy] = meshgrid(1 : sz(2), 1 : sz(1));

        % registration paths
        [rpath, fname] = fileparts(movpath);
        fname = pvParsefn(fname);
        regpath = fullfile(fileparts(rpath), 'demonsreg');

        % Make demonsreg folder if necessary
        if ~exist(regpath, 'dir')
            mkdir(regpath);
        end 
        
        if iter == 1
            % Initialize registered data (using the grid file xx to get size)
            fprintf('Initializing array...')
            tic;
            data_reg = uint16(zeros(sz(1), sz(2)));
            data_reg = repmat(data_reg, [1, 1, sz(3)]);
            t = toc;
            fprintf(' Done. %i seconds.\n', round(t));

            % Initialize a cell for registered data
            data_reg_cell = cell(nchunks, 1);
        end
        
        % Parallel or not
        if p.parfor
            M = p.nworkers;
        else
            M = 0;
        end

        % Slice data
        fprintf('Slicing data...')
        tic;
        for c = 1 : nchunks
            data_reg_cell{c} = im(:,:,(c-1) * p.chunksize + 1 :...
                min(c * p.chunksize, sz(3)));
        end
        t = toc;
        fprintf(' Done. %i seconds.\n', round(t));
        
        fprintf('Parallel registration...')
        tic
        
        % Parallel processing
        parfor (c = 1 : nchunks, M)
            % savepath for warp files
            savepath = [regpath, '\Transforms\', fname, sprintf('_itr%i', iter), '_Dtransform_', sprintf('%02d', c)];

            % Process data
            data_reg_cell{c} = ...
                uint16(pvDemonsRegOneshotTiffCore(data_reg_cell{c}, savepath,...
                ref, xx, yy, 'ref_downsample_xy', p.binxy,...
                'hp_norm_sigmas', p.hp_norm_sigmas, 'savewarp', p.savewarp, ...
                'medfilt2size', p.medfilt2size, 'binbeforehighpassnorm', p.binbeforehighpassnorm,...
                'highpassnorm', p.highpassnorm, 'edges', p.edges, 'itr', p.itr, 'PyramidLevels', p.PyramidLevels,...
                'AccumulatedFieldSmoothing', p.AccumulatedFieldSmoothing));

            if ~p.parfor
                fprintf('Chunk %i/%i registered.\n', c, nchunks);
            end
        end
        t = toc;
        fprintf(' Done. %i seconds.\n', round(t));

        % Reconstruct image stack
        fprintf('Reconstruct image stack...')
        tic;
        for c = 1 : nchunks
            % time indices
            i_start = (c-1) * p.chunksize + 1;
            i_end = min(c * p.chunksize, sz(3));

            % Reconstruct
            data_reg(:,:,i_start : i_end) = data_reg_cell{c};

            if ~p.parfor
                fprintf('Chunk %i/%i reconstructed.\n', c, nchunks);
            end
        end
    end
    
    % Free memory
    clear data_reg_cell
    t = toc;
    fprintf(' Done. %i seconds.\n', round(t));
    
    fprintf('[Output]\n');

    % Post hoc median filter
    if p.posthocmedian
        fprintf('Posthoc median filter...')
        tic;
        
        % filter
        data_reg = movmedian(data_reg, p.posthocmedianwindow, 3);
        
        t = toc;
        fprintf(' Done. %i seconds.\n', round(t));
    end
    
    % PCA clean
    if p.PCAclean
        fprintf('PCA cleaning...')
        tic;
        
        % PCA clean
        dataclass = class(data_reg);
        data_reg = double(data_reg);
        data_reg = PCAClean(data_reg, p.PCs_clean);
        data_reg = cast(data_reg, dataclass);
        
        t = toc;
        fprintf(' Done. %i seconds.\n', round(t));
    end
    fprintf('Saving final stack...')
    
    % Save
    writetiff(data_reg, outpaths{i});
    save(matname, '-v7.3', 'sz', 'p', 'ref');
    
    % Preview
    if p.topreview
        fprintf('Binning preview xy... ');
        tic;
        im_preview = binxy(data_reg, p.previewbin_xyt(1));
        fprintf('Done. Elapsed time (s) = %i\n', round(toc));
        
        fprintf('Binning preview t... ');
        tic;
        im_preview = bint(im_preview, p.previewbin_xyt(2));
        fprintf('Done. Elapsed time (s) = %i\n', round(toc));
        
        % File name
        previewfn_out = sprintf('%s_preview.tif', outpaths{i}(1:end-4));
        writetiff(im_preview, previewfn_out);
    end

    % Segmentation file
    if p.toseg
        % meantifpath
        [outputfolder, ~, ~] = fileparts(outpaths{i});
        meanpath = fullfile(outputfolder, sprintf('%s_%s-%03d_Ch%i_toseg.tif', mouse, date, p.runs(i), p.pmt));
        meanim = median(data_reg(:,:,p.refoffset : p.refoffset + p.refsize - 1),3);
        meanim = double(meanim);
        meanim = imresize(meanim,1/p.binxy);
        
        meanim = medfilt2(meanim, p.medfilt2size, 'symmetric');
        meanim_prime = meanim - imgaussfilt(double(meanim),p.hp_norm_sigmas(1));
        ln_meanim = meanim_prime ./ (imgaussfilt(meanim_prime.^2,p.hp_norm_sigmas(2)) .^ (1/2));
        ln_meanim(isnan(ln_meanim)) = 0;
    
        writetiff(ln_meanim, meanpath);
    end
    fprintf('Saving done. \n');
end


end