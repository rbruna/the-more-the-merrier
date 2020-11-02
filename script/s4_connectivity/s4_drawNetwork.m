clc
clear
close all

% Sets the paths.
config.path.corr  = '../../stats/correlations/';
config.path.figs  = '../../figs/violins/';
config.path.patt  = '*.mat';

% Defines the label of the correlations.
config.label      = '99';


template     = load ( 'template_AAL.mat' );
atlas        = template.atlas;
atlas.name   = atlas.name  ( [ 1: 40 43: 70 79: 90 ] );
atlas.nick   = atlas.nick  ( [ 1: 40 43: 70 79: 90 ] );
atlas.pos    = atlas.pos   ( [ 1: 40 43: 70 79: 90 ], : );
atlas.order  = atlas.order ( ismember ( atlas.order, [ 1: 40 43: 70 79: 90 ] ) );
[ ~, order ] = ismember ( atlas.order, sort ( atlas.order ) );
atlas.name   = atlas.name  ( order );
atlas.nick   = atlas.nick  ( order );
atlas.pos    = atlas.pos   ( order, : );
atlas.order  = order;


% Adds the 'functions' folder to the path.
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );
addpath ( sprintf ( '%s/functions/', pwd ) );


% Loads both sets of correlations results.
data              = load ( sprintf ( '%splv_%s', config.path.corr, config.label ) );

data_ciplv   = data.ciplv_inter ( :, 2, 2 );
narea        = size ( data.ciplv_intra, 1 );
mindex       = triu ( true ( narea ), 1 );
matrix_ciplv = zeros ( narea );
matrix_ciplv ( mindex ) = data_ciplv;
pthres_ciplv = sort ( data_ciplv, 'descend' );
pthres_ciplv = pthres_ciplv ( ceil ( numel ( pthres_ciplv ) / 50 ) );
nthres_ciplv = sort ( data_ciplv, 'ascend' );
nthres_ciplv = nthres_ciplv ( ceil ( numel ( nthres_ciplv ) / 50 ) );
matrix_ciplv = matrix_ciplv + matrix_ciplv';
matrix_ciplv = matrix_ciplv ( atlas.order, atlas.order );

data_plv     = data.plv_inter ( :, 2, 2 );
narea        = size ( data.plv_intra, 1 );
mindex       = triu ( true ( narea ), 1 );
matrix_plv   = zeros ( narea );
matrix_plv   ( mindex ) = data_plv;
pthres_plv   = sort ( data_plv, 'descend' );
pthres_plv   = pthres_plv ( ceil ( numel ( pthres_plv ) / 50 ) );
nthres_plv   = sort ( data_plv, 'ascend' );
nthres_plv   = nthres_plv ( ceil ( numel ( nthres_plv ) / 50 ) );
matrix_plv   = matrix_plv + matrix_plv';
matrix_plv   = matrix_plv ( atlas.order, atlas.order );


cfg        = [];
cfg.label  = atlas.nick;
cfg.nodes  = atlas.pos * 1000;
cfg.width  = 20;
cfg.groups = '';

my_drawNetworkFull ( ( matrix_plv >= pthres_plv ) - ( matrix_plv <= nthres_plv ), cfg )

print ( '-dpng', '-r300', 'Distribution of higher and lower correlations for PLV' )

my_drawNetworkFull ( ( matrix_ciplv >= pthres_ciplv ) - ( matrix_ciplv <= nthres_ciplv ), cfg )

print ( '-dpng', '-r300', 'Distribution of higher and lower correlations for ciPLV' )
