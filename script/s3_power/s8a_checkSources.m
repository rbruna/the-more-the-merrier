clc
clear
close all

% Sets the path.
config.path.filt     = '../../data/power/fourier/';
config.path.figs     = '../../figs/power/';
config.path.patt     = '*.mat';

config.savefig       = false;


% Saves the original path.
pathvar = path;

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );

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
    
    % Loads the MRI data and extracts the masks.
    filtdata      = load ( sprintf ( '%s%s', config.path.filt, files ( findex ).name ) );
    
    alphapow      = sum ( filtdata.band.sources.pow ( :, filtdata.band.sources.freq >= 8 & filtdata.band.sources.freq < 12 ), 2 );
    allpow        = sum ( filtdata.band.sources.pow, 2 );
    
    relpow        = alphapow ./ allpow;
    
    grid          = template.grid;
    grid.pow      = relpow;
    
    cfg           = [];
    cfg.parameter = 'pow';
    
    mrigrid       = ft_sourceinterpolate ( cfg, grid, template.mri );
    
    cfg           = [];
    cfg.funparameter = 'pow';
    cfg.location  = [ +0.00, -0.06, +0.00 ];
    
    ft_sourceplot ( cfg, mrigrid );
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s_%s_w%s_r%s_%s', config.path.figs, filtdata.subject, filtdata.task, filtdata.channel, filtdata.whitener, filtdata.lambda ( 1: end - 1 ), filtdata.taper ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s_%s_w%s_r%s_%s.fig', config.path.figs, filtdata.subject, filtdata.task, filtdata.channel, filtdata.whitener, filtdata.lambda ( 1: end - 1 ), filtdata.taper ) )
    end
    
    close all
end

% Restores the original path.
path ( pathvar );
