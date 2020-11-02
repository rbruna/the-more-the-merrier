clc
clear
close all

% Sets the paths.
config.path.meg  = '../../data/segments/';
config.path.lead = '../../data/sources/leadfield/';
config.path.filt = '../../data/sources/beamformers_3D/';
config.path.patt = '*.mat';

% Action when the task has already been processed.
config.overwrite = false;

% Defines the channels to use.
config.channel   = 'EEG';

% Defines the regularization parameter.
config.lambda    = '25%';

% Level of whitening: 0 (no whitening), 1 (scaling), 2 (PCA).
config.whitener  = 1;

% Defines the bands.
config.bands (1).name  = 'Theta';
config.bands (1).order = 1800;
config.bands (1).edges = [  4.0  8.0 ];
config.bands (1).fs    = 1000;

config.bands (2).name  = 'Alpha';
config.bands (2).order = 1800;
config.bands (2).edges = [  8.0 12.0 ];
config.bands (2).fs    = 1000;

config.bands (3).name  = 'Beta1';
config.bands (3).order = 1800;
config.bands (3).edges = [ 12.0 20.0 ];
config.bands (3).fs    = 1000;

config.bands (4).name  = 'Beta2';
config.bands (4).order = 1800;
config.bands (4).edges = [ 20.0 30.0 ];
config.bands (4).fs    = 1000;

config.bands (5).name  = 'Beta';
config.bands (5).order = 1800;
config.bands (5).edges = [ 12.0 30.0 ];
config.bands (5).fs    = 1000;

config.bands (6).name  = 'Gamma';
config.bands (6).order = 1800;
config.bands (6).edges = [ 30.0 45.0 ];
config.bands (6).fs    = 1000;

config.bands (7).name  = 'Broadband';
config.bands (7).order = 1800;
config.bands (7).edges = [  2.0 45.0 ];
config.bands (7).fs    = 1000;


% Creates and output folder, if needed.
if ~exist ( config.path.filt, 'dir' ), mkdir ( config.path.filt ); end

% Saves the original path.
pathvar = path;

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
ft_path
ft_defaults

% Disables the FT feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';

% Adds the FreeSrufer toolbox to the path.
ft_hastoolbox ( 'spm8', 1, 1 );


% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.meg, config.path.patt ) );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Preloads the data.
    megdata         = load ( sprintf ( '%s%s', config.path.meg,  files ( findex ).name ), 'subject', 'task', 'stage', 'channel' );
    
    % Creates the output structure and fills it with the metadata.
    output          = [];
    output.subject  = megdata.subject;
    output.task     = megdata.task;
    output.stage    = megdata.stage;
    output.channel  = config.channel;
    output.whitener = 'None';
    output.lambda   = config.lambda;
    
    % Defines the whitener type.
    switch config.whitener
        case 0, output.whitener = 'None';
        case 1, output.whitener = 'Scaling';
        case 2, output.whitener = 'PCA';
        case 3, output.whitener = 'MNE';
    end
    
    if ~exist ( sprintf ( '%s%s.mat', config.path.lead, megdata.subject ), 'file' )
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (no lead field).\n', output.subject, output.task, output.channel );
        continue
    end
    
    if ~config.overwrite && exist ( sprintf ( '%s%s_%s%s_%s_w%s_r%s.mat', config.path.filt, output.subject, output.task, output.stage, output.channel, output.whitener, output.lambda ( 1: end - 1 ) ), 'file' )
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (already calculated).\n', output.subject, output.task, output.channel );
        continue
    end
    
    fprintf ( 1, 'Calculating beamformers for subject ''%s'', task ''%s'', channel group ''%s''.\n', output.subject, output.task, output.channel );
    
    
    % Loads the data and the leadfield.
    megdata         = load ( sprintf ( '%s%s', config.path.meg,  files ( findex ).name ) );
    leaddata        = load ( sprintf ( '%s%s.mat', config.path.lead, megdata.subject ) );
    
    
    % Extracts the data padding.
    padding        = megdata.trialinfo.trialpad (1);
    
    % Goes through each band.
    for band =  1: numel ( config.bands )
        
        % Gets the current band information.
        bandinfo             = config.bands ( band );
        
        % Gets the leadfield definition.
        grid                 = leaddata.grid;
        
        fprintf ( 1, '  Working in ''%s'' band (%0.1f to %0.1f Hz).\n', bandinfo.name, bandinfo.edges );
        fprintf ( 1, '    Filtering the data.\n' );
        
        % Gets the data.
        banddata             = megdata.trialdata;
        
        % Initializes the list of channels.
        channel              = banddata.label;
        
        % Removes the bad channels, if any.
        channel              = setdiff ( channel, megdata.chaninfo.bad );
        
        % Removes the channels not present in the lead field.
        channel              = intersect ( channel, grid.label, 'stable' );
        
        % Keeps only the channels selected to be used.
        channel              = ft_channelselection ( config.channel, channel );
        
        cfg                  = [];
        cfg.channel          = channel;
        
        banddata             = ft_selectdata ( cfg, banddata );
        
        
        % Filters the data in the selected band and removes the padding.
        fir                  = fir1 ( bandinfo.order, bandinfo.edges / ( bandinfo.fs / 2 ) );
        banddata             = ft_myfiltfilt ( fir, 1, banddata );
        
        cfg                  = [];
        cfg.toilim (1)       = banddata.time {1} (1)     + padding;
        cfg.toilim (2)       = banddata.time {1} ( end ) - padding;
        cfg.feedback         = 'no';
        
        banddata             = ft_redefinetrial ( cfg, banddata );
        
        
        fprintf ( 1, '    Calculating the covariance matrix.\n' );
        
        % Calculates the covariance matrix.
        cfg                  = [];
        cfg.channel          = grid.label;
        cfg.covariance       = 'yes';
        cfg.covariancewindow = 'all';
        cfg.keeptrials       = 'no';
        cfg.removemean       = 'yes';
        cfg.feedback         = 'no';
        
        timelock             = ft_timelockanalysis ( cfg, banddata );
        
        
        fprintf ( 1, '    Creating the beamformer spatial filter.\n' );
        
        % Generates the data whitener, if requested.
        whitener             = my_whitener ( timelock, config.whitener, true );
        
        % Calculates the spatial filter.
        cfg                  = [];
        cfg.grid             = grid;
        cfg.keepori          = 'yes';
        cfg.keepleadfield    = 'no';
        cfg.method           = 'lcmv';
        
        cfg.lcmv.lambda      = config.lambda;
        cfg.lcmv.projectmom  = 'no';
        cfg.lcmv.keepmom     = 'no';
        cfg.lcmv.keepcov     = 'yes';
        cfg.lcmv.keepfilter  = 'yes';
        cfg.lcmv.feedback    = 'no';
        cfg.lcmv.subspace    = whitener;
        
        sources              = my_sourceanalysis ( cfg, timelock );
        
        
        % Sets the output
        bandinfo.label       = sources.label;
        bandinfo.sources     = sources;
        
        output.band ( band ) = bandinfo;
    end
    
    fprintf ( 1, '  Saving the data.\n' );
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s%s_%s_w%s_r%s', config.path.filt, output.subject, output.task, output.stage, output.channel, output.whitener, output.lambda ( 1: end - 1 ) ), '-struct', 'output' );
end
