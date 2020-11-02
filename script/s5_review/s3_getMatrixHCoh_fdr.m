clc
clear
close all

% Sets the paths.
config.path.conn  = '../../data/connectivity/hcoh_99/';
config.path.patt  = '*.mat';

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma' };
config.bandname   = { '\theta' '\alpha' 'Low \beta' 'High \beta' '\gamma' };
config.metrics    = { 'Average' 'RMS' 'First PCA' 'Centroid' 'Average\ntime series' };

% Defines the cotical ROIs.
config.roi        = [ 1: 40 43: 70 79: 90 ];

addpath ( sprintf ( '%s/functions/', pwd ) );
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );


% Gets the list of subjects.
files = dir ( sprintf ( '%s%s', config.path.conn, config.path.patt ) );

% Initializes the data cell array.
data  = cell ( numel ( files ), 1 );

% Goes through each subject.
for findex = 1: numel ( files )
    
    % Loads the data.
    datum = load ( sprintf ( '%s%s', config.path.conn, files ( findex ).name ) );
    
    % Gets the desired bands.
    order = my_matchstr ( { datum.band.name }, config.bands );
    band  = datum.band ( order );
    
    % Fills the metadata.
    datum.band  = { band.name };
    datum.label = band (1).label;
    datum.nick  = band (1).nick;
    datum.appr  = { 'Mean PLV' 'RMS' 'PCA' 'Centroid' 'Average time series' };
    
    % Gets the connectivity for each label.
    for bindex = 1: numel ( band )
        datum.plv   ( :, :, 1, bindex ) = band ( bindex ).hcoh_avg_all;
        datum.plv   ( :, :, 2, bindex ) = band ( bindex ).hcoh_rms_all;
        datum.plv   ( :, :, 3, bindex ) = band ( bindex ).hcoh_avg_pc;
        datum.plv   ( :, :, 4, bindex ) = band ( bindex ).hcoh_avg_centroid;
        datum.plv   ( :, :, 5, bindex ) = band ( bindex ).hcoh_avg_mean;
        datum.ciplv ( :, :, 1, bindex ) = band ( bindex ).cihcoh_avg_all;
        datum.ciplv ( :, :, 2, bindex ) = band ( bindex ).cihcoh_rms_all;
        datum.ciplv ( :, :, 3, bindex ) = band ( bindex ).cihcoh_avg_pc;
        datum.ciplv ( :, :, 4, bindex ) = band ( bindex ).cihcoh_avg_centroid;
        datum.ciplv ( :, :, 5, bindex ) = band ( bindex ).cihcoh_avg_mean;
    end
    
    datum.label = datum.label ( config.roi );
    datum.nick  = datum.nick  ( config.roi );
    datum.plv   = datum.plv   ( config.roi, config.roi, :, : );
    datum.ciplv = datum.ciplv ( config.roi, config.roi, :, : );
    
    data { findex } = datum;
end

% Concatenates all the datas.
data     = cat ( 1, data {:} );
data ( strcmp ( { data.subject }, 'mq-06' ) ) = [];

% Lists the subjects and channel groups.
subjs    = unique ( { data.subject } );
chans    = unique ( { data.channel } );


% Gets the subject and channel type indexes.
sind     = my_matchstr ( subjs, { data.subject } );
cind     = my_matchstr ( chans, { data.channel } );

% Reserves memory for the channel types subdata.
cdata    = cell ( numel ( chans ), 1 );

% Goes through each channel type.
for cindex = 1: numel ( chans )
    
    % Gets only the selected subdata.
    sdata    = data ( cind == cindex );
    ssind    = sind  ( cind == cindex );
    
    % Stores the subdata.
    cdata { cindex } = sdata;
end

% Concatenates all the subdatas.
cdata    = cat ( 2, cdata {:} );


% Gets the size of the data.
ndata    = size ( cdata, 1 );
narea    = numel ( cdata (1).label );
nband    = numel ( cdata (1).band );
napp     = numel ( cdata (1).appr );

% Calculates the correlation from the normalized data.
plveeg   = cat ( 5, cdata ( :, 1 ).plv );
plvmeg   = cat ( 5, cdata ( :, 2 ).plv );
ciplveeg = cat ( 5, cdata ( :, 1 ).ciplv );
ciplvmeg = cat ( 5, cdata ( :, 2 ).ciplv );

plveeg   = zscore ( plveeg, [], 5 );
plvmeg   = zscore ( plvmeg, [], 5 );
r_plv    = sum ( plveeg .* plvmeg, 5 ) / ( ndata - 1 );

ciplveeg = zscore ( ciplveeg, [], 5 );
ciplvmeg = zscore ( ciplvmeg, [], 5 );
r_ciplv  = sum ( ciplveeg .* ciplvmeg, 5 ) / ( ndata - 1 );

% Calculates the associated p-value.
t_plv    = r_plv .* sqrt ( ( ndata - 2 ) ./ ( 1 - r_plv .^ 2 ) );
p_plv    = tcdf ( -t_plv, ndata - 2 );

t_ciplv  = r_ciplv .* sqrt ( ( ndata - 2 ) ./ ( 1 - r_ciplv .^ 2 ) );
p_ciplv  = tcdf ( -t_ciplv, ndata - 2 );


% Marks the upper triangular and the diagonal.
iinter   = find ( triu ( true ( narea ), 1 ) );
iintra   = find ( eye ( narea ) == 1 );

% Rewrites the data as columns.
r_plv    = reshape ( r_plv,   [], napp, nband );
p_plv    = reshape ( p_plv,   [], napp, nband );
r_ciplv  = reshape ( r_ciplv, [], napp, nband );
p_ciplv  = reshape ( p_ciplv, [], napp, nband );


round ( 100 * squeeze ( sum ( p_plv   ( iinter, :, : ) <= 0.05 ) )' / numel ( iinter ) )
% round ( 100 * squeeze ( sum ( p_plv   ( iintra, :, : ) <= 0.05 ) ) / numel ( iintra ) )
round ( 100 * squeeze ( sum ( p_ciplv ( iinter, :, : ) <= 0.05 ) )' / numel ( iinter ) )
% round ( 100 * squeeze ( sum ( p_ciplv ( iintra, :, : ) <= 0.05 ) ) / numel ( iintra ) )


q_plv   = sort ( p_plv   ( iinter, :, : ), 1 ) <= 0.05 .* ( 1: numel ( iinter ) )' / numel ( iinter );
q_ciplv = sort ( p_ciplv ( iinter, :, : ), 1 ) <= 0.05 .* ( 1: numel ( iinter ) )' / numel ( iinter );
q_plv   = flipud ( q_plv );
q_ciplv = flipud ( q_ciplv );
q_plv   = cumsum ( q_plv, 1 );
q_ciplv = cumsum ( q_ciplv, 1 );
q_plv   = flipud ( q_plv ) > 0;
q_ciplv = flipud ( q_ciplv ) > 0;

round ( 100 * squeeze ( sum ( q_plv   ) )' / numel ( iinter ) )
round ( 100 * squeeze ( sum ( q_ciplv ) )' / numel ( iinter ) )
