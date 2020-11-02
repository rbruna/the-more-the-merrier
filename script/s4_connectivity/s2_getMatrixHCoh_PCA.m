clc
clear
close all

% Sets the paths.
config.path.conn  = '../../data/connectivity/hcoh_99/';
config.path.figs  = '../../figs/violins/';
config.path.patt  = '*.mat';

% Defines a label for the figures.
config.label      = '99';

% Defines the bands to use.
config.bands      = { 'Theta' 'Alpha' 'Beta1', 'Beta2' 'Gamma' };
config.bandname   = { '\theta' '\alpha' 'Low \beta' 'High \beta' '\gamma' };
config.metrics    = { 'RMS' 'PCA\n99.9%%' 'PCA\n99%%' 'PCA\n95%%' 'First PCA' };

% Defines the cotical ROIs.
config.roi        = [ 1: 40 43: 70 79: 90 ];

addpath ( sprintf ( '%s/functions/', pwd ) );
addpath ( sprintf ( '%s/functions/', fileparts ( pwd ) ) );


% Creates the output folder, if required.
if ~exist ( config.path.figs, 'dir' ), mkdir ( config.path.figs ), end


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
        datum.plv   ( :, :, 1, bindex ) = band ( bindex ).hcoh_rms_all;
        datum.plv   ( :, :, 2, bindex ) = band ( bindex ).hcoh_rms_pc999;
        datum.plv   ( :, :, 3, bindex ) = band ( bindex ).hcoh_rms_pc99;
        datum.plv   ( :, :, 4, bindex ) = band ( bindex ).hcoh_rms_pc95;
        datum.plv   ( :, :, 5, bindex ) = band ( bindex ).hcoh_rms_pc;
        datum.ciplv ( :, :, 1, bindex ) = band ( bindex ).cihcoh_rms_all;
        datum.ciplv ( :, :, 2, bindex ) = band ( bindex ).cihcoh_rms_pc999;
        datum.ciplv ( :, :, 3, bindex ) = band ( bindex ).cihcoh_rms_pc99;
        datum.ciplv ( :, :, 4, bindex ) = band ( bindex ).cihcoh_rms_pc95;
        datum.ciplv ( :, :, 5, bindex ) = band ( bindex ).cihcoh_rms_pc;
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
    sdata  = data ( cind == cindex );
    ssind  = sind  ( cind == cindex );
    
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
inter    = find ( triu ( true ( narea ), 1 ) );
intra    = find ( eye ( narea ) == 1 );

% Rewrites the data as columns.
r_plv    = reshape ( r_plv,   [], napp, nband );
p_plv    = reshape ( p_plv,   [], napp, nband );
r_ciplv  = reshape ( r_ciplv, [], napp, nband );
p_ciplv  = reshape ( p_ciplv, [], napp, nband );

% Calculates the significance threshold.
t_thres  = tinv ( 0.95, ndata - 2 );
r_int    = 0: 0.001: 1;
t_int    = r_int .* sqrt ( ( ndata - 2 ) ./ ( 1 - r_int .^ 2 ) );
r_thres  = r_int ( abs ( t_int - t_thres ) == min ( abs ( t_int - t_thres ) ) );


height = [  5.00  5.00 12.00  0.00 ];
height (4) = 2.7 * nband - 0.10 + 0.80;

figure ( 'Units', 'centimeters' )
drawnow
set ( gcf, 'Position', height )

for bindex = 1: nband
    
    % Calculates the vertical offset.
    offset = ( nband - bindex ) * [ 0.00 2.70 0.00 0.00 ] + [ 0.00 0.80 0.00 0.00 ];
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 0.45 2.00 7.75 0.50 ] + offset )
    
    text ( 0, 0, sprintf ( 'Inter-area - %s band', config.bandname { bindex } ), 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'FontSize', 8 )
    xlim ( [ +0 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 8.80 2.00 3.10 0.50 ] + offset )
    
    text ( 0, 0, sprintf ( 'Intra-area - %s band', config.bandname { bindex } ), 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'FontSize', 8 )
    xlim ( [ +0 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 0.45 0.20 7.75 2.00 ] + offset )
    
    plot ( [ 0 1 ] * ( napp + 1 ), [ 1 1 ] * r_thres, 'LineStyle', '--', 'Color', [ 0.7 0.7 0.7 ] )
    my_violinplot ( num2cell ( r_plv ( inter, :, bindex ), 1 ), [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    set ( gca, 'YTick', [ -1 0 +1 ] );
    set ( gca, 'XTick', [] );
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 8.80 0.20 3.10 2.00 ] + offset )
    
    plot ( [ 0 1 ] * ( napp + 1 ), [ 1 1 ] * r_thres, 'LineStyle', '--', 'Color', [ 0.7 0.7 0.7 ] )
    my_violinplot ( num2cell ( r_plv ( intra, 1: 2, bindex ), 1 ), [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    set ( gca, 'YTick', [ -1 0 +1 ] );
    set ( gca, 'XTick', [] );
end

% Creates a new axis for the name of the metrics.
axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8, 'Visible', 'off' )
drawnow
set ( gca, 'Position', [ 0.45 0.00 7.75 1.00 ] )

% Goes through each metric.
for mindex = 1: size ( r_plv, 2 )
    
    % Creates the label.
    text ( mindex, 0, sprintf ( config.metrics { mindex } ), 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', 'FontSize', 8 )
end

% Sets the limits.
ylim ( [ -1 +1 ] )
xlim ( [ 0.5 5.5 ] )


% Creates a new axis for the name of the metrics.
axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8, 'Visible', 'off' )
drawnow
set ( gca, 'Position', [ 8.80 0.00 3.10 1.00 ] )

% Goes through each metric.
for mindex = 1: 2
    
    % Creates the label.
    text ( mindex, 0, sprintf ( config.metrics { mindex } ), 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', 'FontSize', 8 )
end

% Sets the limits.
ylim ( [ -1 +1 ] )
xlim ( [ 0.5 2.5 ] )

drawnow
print ( '-dpng', '-r300', sprintf ( '%sViolins - HCoh %s PCA.png', config.path.figs, config.label ) )
% close


figure ( 'Units', 'centimeters' )
drawnow
set ( gcf, 'Position', height )

for bindex = 1: nband
    
    % Calculates the vertical offset.
    offset = ( nband - bindex ) * [ 0.00 2.70 0.00 0.00 ] + [ 0.00 0.80 0.00 0.00 ];
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 0.45 2.00 7.75 0.50 ] + offset )
    
    text ( 0, 0, sprintf ( 'Inter-area - %s band', config.bandname { bindex } ), 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'FontSize', 8 )
    xlim ( [ +0 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 8.80 2.00 3.10 0.50 ] + offset )
    
    text ( 0, 0, sprintf ( 'Intra-area - %s band', config.bandname { bindex } ), 'HorizontalAlign', 'left', 'VerticalAlign', 'bottom', 'FontSize', 8 )
    xlim ( [ +0 +1 ] )
    ylim ( [ -1 +1 ] )
    axis off
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 0.45 0.20 7.75 2.00 ] + offset )
    
    plot ( [ 0 1 ] * ( napp + 1 ), [ 1 1 ] * r_thres, 'LineStyle', '--', 'Color', [ 0.7 0.7 0.7 ] )
    my_violinplot ( num2cell ( r_ciplv ( inter, :, bindex ), 1 ), [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    set ( gca, 'YTick', [ -1 0 +1 ] );
    set ( gca, 'XTick', [] );
    
    axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8 )
    drawnow
    set ( gca, 'Position', [ 8.80 0.20 3.10 2.00 ] + offset )
    
    plot ( [ 0 1 ] * ( napp + 1 ), [ 1 1 ] * r_thres, 'LineStyle', '--', 'Color', [ 0.7 0.7 0.7 ] )
    my_violinplot ( num2cell ( r_ciplv ( intra, 1: 2, bindex ), 1 ), [ -1 +1 ] )
    ylim ( [ -1 +1 ] )
    set ( gca, 'YTick', [ -1 0 +1 ] );
    set ( gca, 'XTick', [] );
end

% Creates a new axis for the name of the metrics.
axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8, 'Visible', 'off' )
drawnow
set ( gca, 'Position', [ 0.45 0.00 7.75 1.00 ] )

% Goes through each metric.
for mindex = 1: size ( r_plv, 2 )
    
    % Creates the label.
    text ( mindex, 0, sprintf ( config.metrics { mindex } ), 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', 'FontSize', 8 )
end

% Sets the limits.
ylim ( [ -1 +1 ] )
xlim ( [ 0.5 5.5 ] )


% Creates a new axis for the name of the metrics.
axes ( 'Units', 'centimeters', 'NextPlot', 'add', 'FontSize', 8, 'Visible', 'off' )
drawnow
set ( gca, 'Position', [ 8.80 0.00 3.10 1.00 ] )

% Goes through each metric.
for mindex = 1: 2
    
    % Creates the label.
    text ( mindex, 0, sprintf ( config.metrics { mindex } ), 'HorizontalAlign', 'center', 'VerticalAlign', 'middle', 'FontSize', 8 )
end

% Sets the limits.
ylim ( [ -1 +1 ] )
xlim ( [ 0.5 2.5 ] )

drawnow
print ( '-dpng', '-r300', sprintf ( '%sViolins - ciHCoh %s PCA.png', config.path.figs, config.label ) )
% close

squeeze ( sum ( p_plv   ( inter, :, : ) <= 0.05 ) ) / numel ( inter )
squeeze ( sum ( p_plv   ( intra, :, : ) <= 0.05 ) ) / numel ( intra )
squeeze ( sum ( p_ciplv ( inter, :, : ) <= 0.05 ) ) / numel ( inter )
squeeze ( sum ( p_ciplv ( intra, :, : ) <= 0.05 ) ) / numel ( intra )
