clc
clear
close all

% Paths IN/OUT.
config.path.meg   = '../../data/segments/';
config.path.pow   = '../../data/spectra_dpss_05/';
config.path.patt  = '*.mat';

% Frequency band to calculate.
config.band       = [  2 45 ];

% Taper to use.
% config.taper      = 'hann';
config.taper      = 'dpss';
config.smoothing  = 0.5;

% Sets the action when the task has already been processed.
config.overwrite  = false;


% Saves the original path.
pathvar = path;

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds, if needed, the FieldTrip folder to the path.
ft_path
ft_defaults

% Disables the FT feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'spm8', 1, 1 );


% Generates the output folder, if needed.
if ~exist ( config.path.pow, 'dir' ), mkdir ( config.path.pow ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.meg, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    
    
    % Pre-loads the data.
    megdata   = load ( sprintf ( '%s%s', config.path.meg, files ( file ).name ), 'subject', 'task', 'stage', 'channel' );
    
    % Gets the text of the message.
    msgtext   = sprintf ( 'subject ''%s'', task ''%s''', megdata.subject, megdata.task );
    if ~isempty ( megdata.stage )
        msgtext   = sprintf ( '%s, stage ''%s''', msgtext, megdata.stage );
    end
    msgtext   = sprintf ( '%s, channel group ''%s''', msgtext, megdata.channel );
    
    if exist ( sprintf ( '%s%s_%s%s_%s.mat', config.path.pow, megdata.subject, megdata.task, megdata.stage, megdata.channel ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring %s (Already calculated).\n', msgtext );
        continue
    end
    
    fprintf ( 1, 'Working with %s.\n', msgtext );
    
    
    fprintf ( 1, '  Loading the data.\n' );
    
    % Loads the data.
    megdata   = load ( sprintf ( '%s%s', config.path.meg, files ( file ).name ) );
    trialdata = megdata.trialdata;
    
    
    fprintf ( 1, '  Calculating power spectrum.\n' );
    
    % Extracts the channel data.
    rawdata      = cat ( 3, trialdata.trial {:} );
    
    % Removes the padding.
    padding      = round ( megdata.trialinfo.trialpad (1) * trialdata.fsample );
    rawdata      = rawdata ( :, padding + 1: end - padding, : );
    
    % Gets the data dimensions.
    nchan        = size ( rawdata, 1 );
    nsample      = size ( rawdata, 2 );
    ntrial       = size ( rawdata, 3 );
    
    
    % Demeans the data.
    rawdata      = bsxfun ( @minus, rawdata, mean ( rawdata, 2 ) );
    
    
    % Generates the Slepian tapers, if requested.
    if strcmp ( config.taper, 'dpss' )
        taper        = dpss ( nsample, nsample * ( config.smoothing / trialdata.fsample ) )';
        taper        = taper ( 1: end - 1, : );
        
    % Otherwise generates a normalized taper.
    else
        taper        = window ( config.taper, nsample )';
        taper        = taper / norm ( taper );
    end
    
    % Gets the data dimensions.
    ntaper       = size ( taper, 1 );
    
    
    % Generates the frequencies vector.
    freqs        = ( 0: nsample - 1 ) / ( nsample / trialdata.fsample );
    
    % Selects the band of frequencies, if provided.
    if isfield ( config, 'band' ) && ~isempty ( config.band )
        freqindex    = freqs >= config.band (1) & freqs <= config.band (2);
    else
        freqindex    = freqs <= trialdata.fsample / 2;
    end
    freqs        = freqs ( freqindex );
    
    % Gets the data dimensions.
    nfreq        = numel ( freqs );
    
    
    % Initializes the frequency matrix.
    powdata      = zeros ( nchan, nfreq, ntrial, ntaper );
    
    % Goes through each taper.
    for tindex = 1: ntaper
        
        % Applies the taper.
        tapdata      = bsxfun ( @times, rawdata, taper ( tindex, : ) );
        
        % Calculates the Fourier transform of the data.
        tapfreq      = fft ( tapdata, [], 2 ) / sqrt ( nsample );
        
        % Keeps only the selected frequencies.
        tapfreq      = tapfreq ( :, freqindex, :, : );
        
        % Calculates the power.
        powdata ( :, :, :, tindex ) = abs ( tapfreq ) .^ 2;
    end
    
    % Averages the power along trials and tapers.
    powdata      = mean ( powdata, 4 );
    powdata      = mean ( powdata, 3 );
    
    % Corrects the amplitude with the imaginary part.
    corrpow      = freqs > 0 & freqs < trialdata.fsample / 2;
    powdata ( :, corrpow ) = 2 * powdata ( :, corrpow );
    
    
    % Stores the spectrum in FieldTrip form.
    freqdata           = [];
    freqdata.label     = trialdata.label;
    freqdata.powspctrm = powdata;
    freqdata.freq      = freqs;
    freqdata.dimord    = 'chan_freq';
    
    
    fprintf ( 1, '  Saving the spectrum.\n' );
    
    % Saves the data.
    powdata          = [];
    powdata.subject  = megdata.subject;
    powdata.task     = megdata.task;
    powdata.stage    = megdata.stage;
    powdata.channel  = megdata.channel;
    powdata.fileinfo = megdata.fileinfo;
    powdata.freqdata = freqdata;
    
    save ( '-v6', sprintf ( '%s%s_%s%s_%s', config.path.pow, powdata.subject, powdata.task, powdata.stage, powdata.channel ), '-struct', 'powdata' );
end

% Restores the original path.
path ( pathvar );
