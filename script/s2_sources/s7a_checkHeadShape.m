clc
clear
close all

% Sets the path.
config.path.lead     = '../../data/sources/leadfield/';
config.path.figs     = '../../figs/headshape/';
config.path.patt     = '*.mat';

% Selects which versions of the figure to save.
config.savefig       = false;
config.savegif       = true;


% Saves the original path.
pathvar = path;

% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );

% Adds, if needed, the FieldTrip folder to the path.
ft_path
ft_defaults

% Disables the FT feedback.
global ft_default;
ft_default.showcallinfo = 'no';
ft_default.checkconfig  = 'silent';

% Adds the FT toolboxes that will be required.
ft_hastoolbox ( 'spm8', 1, 1 );
ft_hastoolbox ( 'iso2mesh', 1, 1 );
ft_hastoolbox ( 'openmeeg', 1, 1 );


% Loads the template.
template = load ( 'template.mat' );
dipoleu  = template.grid.inside & template.grid.pos ( :, 3 ) >= 0;
dipoled  = template.grid.inside & template.grid.pos ( :, 3 ) <  0;
dipoler  = template.grid.inside & template.grid.pos ( :, 1 ) >= 0;
dipolel  = template.grid.inside & template.grid.pos ( :, 1 ) <  0;

% Generates the output folder, if needed.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ); end

% Gets the files list.
files = dir ( sprintf ( '%s%s', config.path.lead, config.path.patt ) );

% Goes through all the files.
for file = 1: numel ( files )
    
    % Loads the MRI data and extracts the masks.
    leaddata      = load ( sprintf ( '%s%s', config.path.lead, files ( file ).name ) );
    
    % If BEM checks the surfaces.
    if numel ( leaddata.mesh.bnd ) == 3
        
        % Checks the triangle files and the geometry file.
        % Will find intersections and self-intersections.
        om_save_tri ( 'brain.tri', leaddata.mesh.bnd (1).pos, leaddata.mesh.bnd (1).tri );
        om_save_tri ( 'skull.tri', leaddata.mesh.bnd (2).pos, leaddata.mesh.bnd (2).tri );
        om_save_tri ( 'scalp.tri', leaddata.mesh.bnd (3).pos, leaddata.mesh.bnd (3).tri );
        om_write_geom ( 'geom.geom', { 'brain.tri', 'skull.tri', 'scalp.tri' } );
        
        system ( 'om_mesh_info -i brain.tri' );
        system ( 'om_mesh_info -i skull.tri' );
        system ( 'om_mesh_info -i scalp.tri' );
        system ( 'om_check_geom -g geom.geom' );
        
        % Deletes the files.
        delete ( 'geom.geom', 'brain.tri', 'skull.tri', 'scalp.tri' );
        
        % Checks that the meshes are closed.
        % A closed mesh has an Euler characteristic of 2.
        fprintf ( 1, 'Subject %s.\n', leaddata.subject );
        fprintf ( 1, 'The Euler characteristic of the first mesh is %i.\n',  mesheuler ( leaddata.mesh.bnd (1).tri ) );
        fprintf ( 1, 'The Euler characteristic of the second mesh is %i.\n', mesheuler ( leaddata.mesh.bnd (2).tri ) );
        fprintf ( 1, 'The Euler characteristic of the third mesh is %i.\n',  mesheuler ( leaddata.mesh.bnd (3).tri ) );
    end
    
    
    % Converts the meshes and the grid to millimeters.
    grid          = leaddata.grid;
    mesh          = leaddata.mesh;
    headshape     = leaddata.headshape;
    
    % Plots the meshes.
    for mindex = 1: numel ( mesh.tissue )
        switch mesh.tissue { mindex }
            case 'brain', meshcolor = 'brain';
            case 'skull', meshcolor = 'white';
            case 'scalp', meshcolor = 'skin';
            otherwise,    meshcolor = 'white';
        end
        
        ft_plot_mesh  ( mesh.bnd ( mindex ), 'facecolor', meshcolor, 'edgecolor', 'none', 'facealpha', .5 );
    end
    hold on
    
    
    % Gets the head shape points, the fiducials and the HPI coils.
    hpiindex  = strncmp ( headshape.label, 'hpi_', 4 );
    hpipos    = headshape.pnt (  hpiindex, : );
    hspos     = headshape.pnt;
    fidpos    = headshape.fid.pnt;
    
    % Plots the head shape.
    plot3 ( hspos  ( :, 1 ), hspos  ( :, 2 ), hspos  ( :, 3 ), '.b' )
    plot3 ( fidpos ( :, 1 ), fidpos ( :, 2 ), fidpos ( :, 3 ), '.k', 'MarkerSize', 20 )
    plot3 ( hpipos ( :, 1 ), hpipos ( :, 2 ), hpipos ( :, 3 ), '.r', 'MarkerSize', 20 )
    
    % Lights the scene.
    set ( gcf, 'Name', leaddata.subject );
    view ( [   90,   0 ] ), camlight
    view ( [ -150,   0 ] ), camlight
    lighting gouraud
    drawnow
    
    % Saves the figure.
    print ( '-dpng', sprintf ( '%s%s.png', config.path.figs, leaddata.subject ) )
    
    if config.savefig
        savefig ( sprintf ( '%s%s.fig', config.path.figs, leaddata.subject ) )
    end
    if config.savegif
        my_savegif ( sprintf ( '%s%s.gif', config.path.figs, leaddata.subject ) )
    end
    
    close all
    clc
end

% Restores the original path.
path ( pathvar );
