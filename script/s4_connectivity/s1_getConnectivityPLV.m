clc
clear
close all

% Sets the paths.
config.path.meg   = '../../data/segments/';
config.path.beam  = '../../data/sources/beamformers/';
config.path.conn  = '../../data/connectivity/plv_1DW/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = false;

% Defines the template to use.
config.template   = 'template_AAL';

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Beta' 'Gamma' };

% Defines the filter to use (or false to use the band-specific filter).
config.filter     = 'Broadband';

% Defines the parameters.
config.downsample = false;
config.keeptrials = false;
config.threshold  = 1;
config.weight     = true;


% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
ft_path
ft_defaults

% Disables the FT feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';


% Loads the template for the label information.
template = load ( config.template, 'grid', 'atlas' );

% Gets the number of sources and areas.
nareas   = numel ( template.atlas.name );
nsources = numel ( template.grid.area ( template.grid.area ~= 0 ) );


% Creates and output folder, if needed.
if ~exist ( config.path.conn, 'dir' ), mkdir ( config.path.conn ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.beam, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    % Loads the MEG data and the beamformer information.
    filename   = files ( file ).name;
    [ ~, basename ] = fileparts ( filename );
    
    beamdata   = load ( sprintf ( '%s%s', config.path.beam,  basename ), 'subject', 'task', 'stage', 'channel', 'whitener', 'lambda' );
    
    if exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s.mat', config.path.conn, beamdata.subject, beamdata.task, beamdata.stage, beamdata.channel, beamdata.whitener, beamdata.lambda ( 1: end -1 ), template.atlas.atlas ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'', whitening ''%s'', regularization %s (already calculated).\n', beamdata.subject, beamdata.task, beamdata.channel, beamdata.whitener, beamdata.lambda );
        continue
    end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', channel group ''%s'', whitening ''%s'', regularization %s.\n', beamdata.subject, beamdata.task, beamdata.channel, beamdata.whitener, beamdata.lambda );
    
    fprintf ( 1, '  Loading data.\n' );
    
    sensfile   = dir ( sprintf ( '%s%s_%s%s*.mat', config.path.meg, beamdata.subject, beamdata.task, beamdata.stage ) );
    sensfile   = sensfile (1).name;
    
    sensdata   = myft_load ( sprintf ( '%s%s', config.path.meg,  sensfile ) );
    beamdata   = load ( sprintf ( '%s%s', config.path.beam, basename ) );
    
    % Gets the number of padding samples of the data.
    padding    = round ( sensdata.trialinfo.trialpad (1) * sensdata.trialdata.fsample );
    
    
    % Gets the beamformer filter for all the bands, if requested.
    if config.filter
        
        fprintf ( 1, '  Calculating the number of directions for each source.\n' );
        
        % Extracts the selected beamformer.
        bffilter     = strcmpi ( { beamdata.band.name }, config.filter );
        source       = beamdata.band ( bffilter ).sources;
        
        % Gets the requested orientations for each source.
        source       = my_selectOri ( source, config.threshold );
        
        % Gets the original map.
        rawmap       = myfun_getFilter ( source, template );
    end
    
    
    % Keeps only the requested bands.
    bandinfos  = beamdata.band;
    bandinfos  = bandinfos ( ismember ( { bandinfos.name }, config.bands ) );
    
    % Reserves memory for the band data.
    bands        = cell ( numel ( bandinfos ), 1 );
    
    % Goes through each band.
    for bindex = 1: numel ( bandinfos )
        
        % Gets the band information.
        band         = bandinfos ( bindex );

        fprintf ( 1, '  Working in ''%s'' band (%0.1f to %0.1f Hz).\n', band.name, band.edges );
        
        
        % Gets the beamformer filter for the band, if requested.
        if ~config.filter
            
            fprintf ( 1, '    Calculating the number of directions for each source.\n' );
            
            % Gets the beamformer filter.
            source       = bandinfos ( bindex ).sources;
            
            % Gets the requested orientations for each source.
            source       = my_selectOri ( source, config.threshold );
            
            % Gets the original map.
            rawmap       = myfun_getFilter ( source, template );
        end
        
        
        fprintf ( 1, '    Filtering the data.\n' );
        
        % Gets the required sensors and filters in the current band.
        myscript_filterSens
        
        % Gets the data dimensions.
        nchans       = size ( banddata, 1 );
        nsamples     = size ( banddata, 2 );
        ntrials      = size ( banddata, 3 );
        
        
        fprintf ( 1, '    Calculating the representative time series for each area.\n' );
        
        % Gets the maps to the principal components.
        pcamaps      = myfun_getPCAMap ( rawmap, banddata );
        
        % Gets the map to the centroid.
        cmap         = myfun_getCentMap ( rawmap, template.atlas );
        
        % Gets the map to the average time series.
        mmap         = myfun_getMeanMap ( rawmap );
        
        % Concatenates all the maps.
        maps         = cat ( 1, rawmap, pcamaps, cmap, mmap );
        
        % If no weighted, removes the weights from the maps.
        if ~config.weight
            for mindex = 1: numel ( maps )
                maps ( mindex ).weight (:) = 1;
            end
        end
        
        
        
        fprintf ( 1, '    Calculating the all-with-all PLV.\n' );
        
        % Reserves memory for the sensor-level complex covariance matrix.
        sensccov     = complex ( zeros ( nchans, nchans, ntrials ) );
        
        % Goes through each trial.
        for tindex = 1: ntrials
            
            % Calculates the sensor-level complex covariance matrix.
            sensccov ( :, :, tindex ) = banddata ( :, :, tindex ) * banddata ( :, :, tindex )' / nsamples;
        end
        
        % Goes through each mapping.
        for mindex = 1: numel ( maps )
            
            % Gets the current map.
            map        = maps ( mindex );
            
            
            fprintf ( 1, '      Working with mapping %s.\n', map.label );
            
            % Gets the size of the problem.
            nsources   = size ( map.filter, 1 );
            
            % Reserves memory for the source-level covariance matrix.
            cplv       = complex ( zeros ( nsources, nsources, ntrials, 'single' ) );
            
            % Goes through each trial.
            for tindex = 1: ntrials
                
                % Gets the sources time-series for the current trial.
                sourcedata   = map.filter * banddata ( :, :, tindex );
                
                % Gets the vector of phases for each time series.
                sourcevector = sourcedata ./ abs ( sourcedata );
                
                % Calculates the PLV as a matrix multiplication.
                cplv ( :, :, tindex )   = sourcevector * sourcevector' / nsamples;
            end
            
            % Stores the all-with all complex Hibert coherence.
            maps ( mindex ).cplv  = cplv;
        end
        
        
        
        fprintf ( 1, '    Calculating the inter-area PLV coherence.\n' );
        
        % Goes through each mapping.
        for mindex = 1: numel ( maps )
            
            % Gets the current map.
            map          = maps ( mindex );
            
            
            fprintf ( 1, '      Working with mapping %s.\n', map.label );
            
            % Combines the per-source FC into per-area FC.
            map          = myfun_combineFC ( map );
            
            % Averages along trials, if required.
            if ~config.keeptrials
                map.plv_rms   = nanmean ( map.plv_rms, 3 );
                map.ciplv_rms = nanmean ( map.ciplv_rms, 3 );
                map.plv_avg   = nanmean ( map.plv_avg, 3 );
                map.ciplv_avg = nanmean ( map.ciplv_avg, 3 );
            end
            
            % Stores the Hibert coherence.
            maps ( mindex ).plv_avg   = map.plv_avg;
            maps ( mindex ).plv_rms   = map.plv_rms;
            maps ( mindex ).ciplv_avg = map.ciplv_avg;
            maps ( mindex ).ciplv_rms = map.ciplv_rms;
        end
        
        
        % Modifies the band structure.
        band             = rmfield ( band, 'sources' );
        band             = rmfield ( band, 'label' );
        band.label       = template.atlas.name;
        band.nick        = template.atlas.nick;
        
        % Goes through each map.
        for mindex = 1: numel ( maps )
            
            % Gets the current map.
            map              = maps ( mindex );
            
            % Adds the Hilber Coherence to the band structure.
            band.( sprintf ( 'plv_avg_%s',   map.label ) ) = map.plv_avg;
            band.( sprintf ( 'plv_rms_%s',   map.label ) ) = map.plv_rms;
            band.( sprintf ( 'ciplv_avg_%s', map.label ) ) = map.ciplv_avg;
            band.( sprintf ( 'ciplv_rms_%s', map.label ) ) = map.ciplv_rms;
        end
        
        % Stores the modified band structure.
        bands { bindex } = band;
    end
    
    % Concatenates all the bands.
    bands        = cat ( 1, bands {:} );
    
    
    fprintf ( 1, '  Saving the connectivity data.\n' );
    
    % Prepares the output.
    conndata          = [];
    conndata.subject  = beamdata.subject;
    conndata.task     = beamdata.task;
    conndata.stage    = beamdata.stage;
    conndata.channel  = beamdata.channel;
    conndata.whitener = beamdata.whitener;
    conndata.lambda   = beamdata.lambda;
    conndata.atlas    = template.atlas.atlas;
    conndata.band     = bands;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s_w%s_r%s_%s', config.path.conn, conndata.subject, conndata.task, conndata.stage, conndata.channel, conndata.whitener, conndata.lambda ( 1: end -1 ), conndata.atlas ), '-struct', 'conndata' );
end
