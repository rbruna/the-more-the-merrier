clc
clear
close all

% Paths IN/OUT.
config.path.meg   = '../../data/segments/';
config.path.beam  = '../../data/sources/beamformers/';
config.path.pow   = '../../data/power/fourier/';
config.path.patt  = '*.mat';

% Action when the task have already been processed.
config.overwrite  = false;

% Name of the beamformer band to use.
config.beamformer = 'Broadband';

% Frequency band to calculate.
config.band       = [  2 45 ];

% Taper to use.
config.taper      = 'hann';


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


% Loads the template.
template = load ( 'template.mat' );

% Creates and output folder, if needed.
if ~exist ( config.path.pow, 'dir' ), mkdir ( config.path.pow ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.beam, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    % Loads the MEG data and the beamformer information.
    filename = files ( file ).name;
    [ ~, basename ] = fileparts ( filename );
    
    beamformer = load ( sprintf ( '%s%s', config.path.beam,  basename ), 'subject', 'task', 'stage', 'channel' );
    
    if exist ( sprintf ( '%s%s_%s_%s.mat', config.path.pow, beamformer.subject, beamformer.task, beamformer.channel ), 'file' ) && ~config.overwrite
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (already calculated).\n', beamformer.subject, beamformer.task, beamformer.channel );
        continue
    end
    
    % Lists the MEG files.
    megfile    = dir ( sprintf ( '%s%s_%s%s*.mat', config.path.meg, beamformer.subject, beamformer.task, beamformer.stage ) );
    
    if isempty ( megfile )
        fprintf ( 1, 'Ignoring subject ''%s'', task ''%s'', channel group ''%s'' (no MEG data).\n', beamformer.subject, beamformer.task, beamformer.channel );
        continue
    end
    
    fprintf ( 1, 'Working on subject ''%s'', task ''%s'', channel group ''%s''.\n', beamformer.subject, beamformer.task, beamformer.channel );
    
    fprintf ( 1, '  Loading data.\n' );
    
    megfile    = megfile (1).name;
    
    megdata    = myft_load ( sprintf ( '%s%s', config.path.meg,  megfile ) );
    beamformer = load ( sprintf ( '%s%s', config.path.beam, basename ) );
    
    % Gets the number of padding samples of the data.
    padding    = round ( megdata.trialinfo.trialpad (1) * megdata.trialdata.fsample );
    
    
    % Keeps only the requested band.
    bandinfo   = beamformer.band;
    bandinfo   = bandinfo ( strcmp ( { bandinfo.name }, config.beamformer ) );
    
    % Gets the beamformer information.
    source        = bandinfo.sources;
    filter        = cat ( 1, source.filter {:} );
    
    
    % Selects only the channels present in the beamformer filter.
    cfg          = [];
    cfg.channel  = source.label;
    cfg.trials   = 'all';
    
    trialdata    = ft_selectdata ( cfg, megdata.trialdata );
    
    % Gets the order of the channels in the beamformer filter.
    chanorder    = my_matchstr ( trialdata.label, source.label );
    
    % Extracts and sorts the channel data.
    rawdata      = cat ( 3, trialdata.trial {:} );
    rawdata      = rawdata ( chanorder, :, : );
    
    % Removes the padding.
    rawdata      = rawdata ( :, padding + 1: end - padding, : );
    
    
    fprintf ( 1, '  Calculating the source-space spectra.\n' );
    
    % Gets the data dimensions.
    nsource      = size ( filter, 1 );
    nchan        = size ( rawdata, 1 );
    nsample      = size ( rawdata, 2 );
    ntrial       = size ( rawdata, 3 );
    
    % Calculates a normalized taper.
    taper        = window ( config.taper, nsample )';
    taper        = taper / norm ( taper );
    
    % Applies the taper.
    rawdata      = bsxfun ( @times, rawdata, taper );
    
    
    % Calculates the Fourier transform of the data.
    freqdata     = fft ( rawdata, [], 2 ) / sqrt ( nsample );
    freqs        = ( 0: nsample - 1 ) / ( nsample / trialdata.fsample );
    
    % Keeps only the first half of the spectrum.
    freqindex    = freqs <= trialdata.fsample / 2;
    freqdata     = freqdata ( :, freqindex, : );
    freqs        = freqs ( freqindex );
    
    % Selects the band of frequencies, if provided.
    if isfield ( config, 'band' ) && ~isempty ( config.band )
        freqindex    = freqs >= config.band (1) & freqs <= config.band (2);
    elseif isfield ( bandinfo, 'edges' )
        freqindex    = freqs >= bandinfo.edges (1) & freqs <= bandinfo.edges (2);
    else
        freqindex    = true ( size ( freqs ) );
    end
    
    % Keeps only the selected frequencies.
    freqdata     = freqdata ( :, freqindex, : );
    freqs        = freqs ( freqindex );
    nfreq        = numel ( freqs );
    
    
    % Goes to source-space.
    freqdata     = filter * freqdata ( :, : );
    freqdata     = reshape ( freqdata, nsource, nfreq, ntrial );
    
    % Calculates the power spectral density.
    powdata      = mean ( abs ( freqdata ) .^ 2, 3 );
    
    % Corrects the calculated power with the imaginary part.
    corrpow      = freqs > 0 & freqs < trialdata.fsample / 2;
    powdata ( :, corrpow ) = 2 * powdata ( :, corrpow );
    
    
    % Creates the sources structure from the template.
    source           = template.grid;
    
    % Stores the spectrum in FieldTrip form.
    pow              = nan ( numel ( source.inside ), nfreq, 'single' );
    pow ( source.inside, : ) = powdata;
    
    source.pow       = pow;
    source.normpow   = bsxfun ( @rdivide, pow, sum ( pow, 2 ) );
    source.freq      = freqs;
    source.dimord    = 'pos_freq';
    
    bandinfo.sources = source;
    
    
    fprintf ( 1, '  Saving the power data.\n' );
    
    % Prepares the output.
    output          = [];
    output.subject  = beamformer.subject;
    output.task     = beamformer.task;
    output.channel  = beamformer.channel;
    output.whitener = beamformer.whitener;
    output.lambda   = beamformer.lambda;
    output.taper    = config.taper;
    output.band     = bandinfo;
    
    % Saves the data.
    save ( '-v6', sprintf ( '%s%s_%s_%s_w%s_r%s_%s', config.path.pow, output.subject, output.task, output.channel, output.whitener, output.lambda ( 1: end - 1 ), output.taper ), '-struct', 'output' );
end

% Restores the original path.
path ( pathvar );
