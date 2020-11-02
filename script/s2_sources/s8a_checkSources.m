clc
clear
close all

% Sets the path.
config.path.filt     = '../../data/sources/beamformers/';
config.path.figs     = '../../figs/beamformers/';
config.path.patt     = '*.mat';

% Selects which versions of the figure to save.
config.savefig       = false;


% Saves the original path.
pathvar = path;

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
ft_path
ft_defaults

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'spm8', 1, 1 );

% Disables the FT feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';



% Gets the template.
template = load ( 'template' );

% Generates the output folder, if needed.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end

% Gets the files list.
files = dir ( sprintf ( '%s%s', config.path.filt, config.path.patt ) );

% Goes through all the files.
for findex = 1: numel ( files )
    
    % Loads the beam former data.
    filtdata      = load ( sprintf ( '%s%s', config.path.filt, files ( findex ).name ) );
    
    % Calculates the relative alpha power.
    alphaband     = filtdata.band ( strcmpi ( 'alpha', { filtdata.band.name } ) );
    allband       = filtdata.band ( strcmpi ( 'broadband',  { filtdata.band.name } ) );
    relpow        = alphaband.sources.pow ./ allband.sources.pow;
    
    grid          = template.grid;
    grid.pow      = relpow;
    
    % Interporlates the sources to the template space.
    cfg           = [];
    cfg.parameter = 'pow';
    
    mrigrid       = ft_sourceinterpolate ( cfg, grid, template.mri );
    
    % Plots the relative alpha power in the sources.
    cfg           = [];
    cfg.funparameter = 'pow';
    cfg.location  = [ +0.00, -0.06, +0.00 ];
    
    ft_sourceplot ( cfg, mrigrid );
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_w%s_r%s.png', config.path.figs, filtdata.subject, filtdata.task, filtdata.stage, filtdata.channel, filtdata.whitener, filtdata.lambda ( 1: end - 1 ) ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_w%s_r%s.fig', config.path.figs, filtdata.subject, filtdata.task, filtdata.stage, filtdata.channel, filtdata.whitener, filtdata.lambda ( 1: end - 1 ) ) )
    end
    
    close all
end

% Restores the original path.
path ( pathvar );
