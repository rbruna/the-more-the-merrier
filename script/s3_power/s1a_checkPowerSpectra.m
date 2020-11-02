clc
clear
close all

% Paths IN/OUT.
config.path.pow   = '../../data/spectra_hann/';
config.path.figs  = '../../figs/spectra_hann/';
config.path.patt  = '*.mat';

config.savefig    = false;


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
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end

% Gets the list of subjects.
files    = dir ( sprintf ( '%s%s', config.path.pow, config.path.patt ) );

% Goes through each subject.
for file = 1: numel ( files )
    
    powdata   = load ( sprintf ( '%s%s', config.path.pow, files ( file ).name ) );
    freqdata  = powdata.freqdata;
    
    
    % Plots only the magnetometers.
    cfg         = [];
    cfg.channel = 'MEGMAG';
    
    magdata     = ft_selectdata ( cfg, freqdata );
    
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
    
    ft_multiplotER ( cfg, magdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGMAG_multi.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGMAG_multi.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
    cfg.xlim    = [ 7 12 ];
    
    ft_topoplotER ( cfg, magdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGMAG_topo.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGMAG_topo.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
    cfg.channel = { 'MEG1631', 'MEG1641', 'MEG1731', 'MEG1741', 'MEG1831', 'MEG1841', 'MEG1911', 'MEG1921', 'MEG1931', 'MEG1941', 'MEG2011', 'MEG2021', 'MEG2031', 'MEG2041', 'MEG2111', 'MEG2121', 'MEG2131', 'MEG2141', 'MEG2231', 'MEG2241', 'MEG2311', 'MEG2321', 'MEG2331', 'MEG2341', 'MEG2441', 'MEG2511', 'MEG2541' };
    
    ft_singleplotER ( cfg, magdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGMAG_O.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGMAG_O.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    
    % Plots only the gradiometers.
    cfg         = [];
    cfg.channel = 'MEGGRAD';
    
    graddata    = ft_selectdata ( cfg, freqdata );
    
    graddata    = ft_combineplanar ( [], graddata );
    
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    
    ft_multiplotER ( cfg, graddata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGGRAD_multi.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGGRAD_multi.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    cfg.xlim    = [ 7 12 ];
    
    ft_topoplotER ( cfg, graddata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGGRAD_topo.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGGRAD_topo.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    cfg.channel = { 'MEG1632+1633', 'MEG1642+1643', 'MEG1732+1733', 'MEG1742+1743', 'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923', 'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2022+2023', 'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113', 'MEG2122+2123', 'MEG2132+2133', 'MEG2142+2143', 'MEG2232+2233', 'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2332+2333', 'MEG2342+2343', 'MEG2442+2443', 'MEG2512+2513', 'MEG2542+2543' };
    
    ft_singleplotER ( cfg, graddata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_MEGGRAD_O.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_MEGGRAD_O.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    
    % Plots only the EEG.
    cfg         = [];
    cfg.channel = 'EEG';
    
    eegdata     = ft_selectdata ( cfg, freqdata );
    
    cfg         = [];
    cfg.layout  = 'neuromagEEG060.mat';
    
    ft_multiplotER ( cfg, eegdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_EEG_multi.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_EEG_multi.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromagEEG060.mat';
    cfg.xlim    = [ 7 12 ];
    
    ft_topoplotER ( cfg, eegdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_EEG_topo.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_EEG_topo.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
    
    cfg         = [];
    cfg.layout  = 'neuromagEEG060.mat';
    cfg.channel = { 'EEG038', 'EEG044', 'EEG045', 'EEG046', 'EEG047', 'EEG048', 'EEG049', 'EEG050', 'EEG051', 'EEG052', 'EEG053', 'EEG054', 'EEG055', 'EEG056', 'EEG057', 'EEG058', 'EEG059', 'EEG060' };
    
    ft_singleplotER ( cfg, eegdata );
    set ( gcf, 'PaperOrientation', 'portrait' )
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s_%s%s_%s_EEG_O.png', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s_%s%s_%s_EEG_O.fig', config.path.figs, powdata.subject, powdata.task, powdata.stage, powdata.channel ) )
    end
    
    close all
end

% Restores the original path.
path ( pathvar );
