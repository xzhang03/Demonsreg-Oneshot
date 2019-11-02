function sbxDemonsRegOneshot(mouse, date, varargin)
% sbxDemonsRegOneshot is a onestep application of demonsreg, with the
% intention of improving speed.

    p = inputParser;

    % Lab variables variables
    addOptional(p, 'target', []);  % Target run to align to, default runs(1)
    addOptional(p, 'server', []);  % Add in the server name as a string
    addOptional(p, 'runs', []);  % Defaults to sbxRuns
    addOptional(p, 'force', false);  % Overwrite files if they exist
    addOptional(p, 'pmt', 0, @isnumeric);  % Which PMT to use for analysis, 0-green, 1-red
    addOptional(p, 'chunksize', 1000); % Chunk size for parallel processing. Decrease if RAM is an issue
    
    % IO variables
    addOptional(p, 'movtype', 'xyreg');  % input type, can be xyreg, sbxreg, or sbx.
    addOptional(p, 'saveastiff', false); % Save as tiff
    
    % Image-processing variables
    addOptional(p, 'edges', [0 0 0 0]); % Edges of registration. Pixels outside the edges will not be registered
    addOptional(p, 'refsize', 500, @isnumeric);  % Set the number of frames from which we make the reference
    addOptional(p, 'refoffset', 500, @isnumeric);  % The offset in frames for the reference image, accounts for weirdness in the first few frames
    addOptional(p, 'ref_downsample_xy', [], @isnumeric);  % Pixels to downsample in xy, will be set to downsample_xy if empty
      
    % Rarely changed image-processing variables
    addOptional(p, 'medfilt2size', [2 2]); % Neighbor area for 2D median filter. Leave empty if no median filter
    addOptional(p, 'highpassnorm', true); % Use highpass filter and local normalization before registration.
    addOptional(p, 'hp_norm_sigmas', [8, 30], @isnumeric);  % sigmas for highpass filter and local normalization
    addOptional(p, 'binbeforehighpassnorm', false); % Bin image before applying spatial filters. Recommended for imaging
                                                    % systems with magenification (e.g., GRIN doublets)                                                    % magnification
  
    % Efficiency variables
    addOptional(p, 'savewarp', true); % Save warp parameteres
    addOptional(p, 'reuseref', false); % Use pre-calculated ref if can find it (must be the right size).
    
    
    % Unpack if needed
    if iscell(varargin) && size(varargin,1) * size(varargin,2) == 1
        varargin = varargin{:};
    end

    parse(p, varargin{:});
    p = p.Results;

%% Input cleanup
if isempty(p.target), p.target = p.runs(1); end
    
%% IO paths
% sbx paths
sbxpaths = cell(length(p.runs),1);
for i = 1 : length(p.runs)
    sbxpaths{i} = sbxPath(mouse, date, p.runs(i),'xyreg', 'server', p.server, 'pmt', p.pmt);
end

% ref paths
[rpath, fname] = fileparts(sbxpaths{p.runs == p.target});        
regpath = [fileparts(rpath) '\reg_demonsreg'];
refname = [regpath '\' fname '_demonsreg.tif'];

% output paths
outpaths = cell(length(p.runs),1);
for i = 1 : length(p.runs)
    outpaths{i} = sprintf('%s_demonsreg-%i.sbxreg', sbxpaths{i}(1:end-7), p.pmt);
end

%% Target
% get ref

if exist(refname, 'file') && p.reuseref
    % Re-use a ref file
    fprintf('Using existing reference image...\n')
    ref = readtiff(refname);
else
    % Calculate ref fresh
    fprintf('Calculating reference image...')
    tic;
    ref = sbxAlignTargetCore(sbxpaths{p.runs == p.target}, p.pmt, p.ref_downsample_xy,...
        p.refoffset, p.refsize, p.edges);
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    writetiff(ref, refname, class(ref));
end

%% Registration
for i = 1:length(sbxpaths)
    % Get the path and info file
    movpath = sbxpaths{i};
    info = sbxInfo(movpath);
    nframes = info.max_idx + 1;

    % Parpool
    openParallel;
    
    % Parallel processing
    nchunks = ceil(nframes/p.chunksize);
    
    % Get the meshgrid to pass to all the workers
    [xx, yy] = meshgrid(1 : info.sz(2), 1 : info.sz(1));
    
    % registration paths
    [rpath, fname] = fileparts(movpath);
    regpath = [fileparts(rpath) '\reg_demonsreg'];
    
    % Initialize registered data (using the grid file xx to get size)
    fprintf('Initializing array...')
    tic;
    data_reg = uint16(zeros(size(xx)));
    data_reg = repmat(data_reg, [1, 1, nframes]);
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    % Initialize a cell for registered data
    data_reg_cell = cell(nchunks, 1);
    
    fprintf('Parallel registration...')
    tic;
    % Parallel processing
    parfor c = 1 : nchunks
        % savepath for warp files
        savepath = [regpath '\Transforms\' fname '_Dtransform_' sprintf('%02d', c)];
        
        % Process data
        data_reg_cell{c} = ...
            uint16(sbxDemonsRegOneshotCore(movpath, savepath, (c-1)*p.chunksize+1, ...
            p.chunksize, ref, xx, yy, 'pmt', p.pmt, 'ref_downsample_xy', p.ref_downsample_xy,...
            'hp_norm_sigmas', p.hp_norm_sigmas, 'savewarp', p.savewarp, ...
            'medfilt2size', p.medfilt2size, 'binbeforehighpassnorm', p.binbeforehighpassnorm,...
            'highpassnorm', p.highpassnorm, 'edges', p.edges));
    end
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    % Reconstruct image stack
    fprintf('Reconstruct image stack...')
    tic;
    for c = 1 : nchunks
        % time indices
        i_start = (c-1) * p.chunksize + 1;
        i_end = min(c * p.chunksize, nframes);
        
        % Reconstruct
        data_reg(:,:,i_start : i_end) = data_reg_cell{c};
        
    end
    % Free memory
    clear data_reg_cell
    t = toc;
    fprintf(' Done. Elapsed time: %i seconds.\n', round(t));
    
    fprintf('Saving final stack...')
    % Save
    if p.saveastiff
        writetiff(data_reg,[outpaths{i}(1:end-7), '.tif'], 'double');  
    else
        sbxWrite(outpaths{i}, data_reg, info, p.force, true);
    end
end


end