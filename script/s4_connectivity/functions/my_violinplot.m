%% data debe ser un cell de n arrays donde n es cada una de las condiciones
% por ejemplo para plotear volumenes de dcl y control. data={vol_dcl,vol_cnt}
%label=cell con nombres de cada grupo en el mismo orden que data
function my_violinplot ( data, lims, varargin )

if nargin < 2
    lims = [ -inf +inf ];
end

cindex = strcmpi ( varargin, 'Color' ) + 1;
if isempty ( cindex )
    colors = get ( gca, 'ColorOrder' );
else
    colors = varargin { cindex };
end

set ( gca, 'NextPlot', 'add' )

for dindex = 1: numel ( data )
    
    % Gets the color for this set of data.
    color = colors ( rem ( dindex - 1, size ( colors, 1 ) ) + 1, : );
    
    % If all NaN skips the data.
    if all ( isnan ( data { dindex } ) )
        continue
    end
    
    % Calculates the probability density function.
%     [fv,voli] = ksdensity(data{i},'BandWidth',.1,'Support',lims);
%     [ fv, voli ] = ksdensity ( data { dindex }, 'BandWidth', .1, 'Support', lims );
    bw = ( ceil ( 100 / numel ( data { dindex } ) ) ) / 20;
    [ fv, voli ] = ksdensity ( data { dindex }, 'BandWidth', bw, 'Support', lims );
    
    % Converts the values in a circle.
    fv = cat ( 2, fv, -fliplr ( fv ) );
    voli = cat ( 2, voli, fliplr ( voli ) );
    
    % Normalizes the PDF.
    fv = fv / ( max ( fv ) ) / 2.5;
    
    % Calculates the meadian and the quartiles.
    Q1 = prctile ( data { dindex }, 25 );
    Q2 = prctile ( data { dindex }, 50 );
    Q3 = prctile ( data { dindex }, 75 );
    
    % Plots the left side of the violin.
    fill ( fv + dindex, voli, 'b', 'FaceColor', color, 'FaceAlpha', 0.3, 'EdgeColor', 'none' );
    plot ( fv + dindex, voli, 'Color', color );
    
    % Plots the quartiles.
    plot( [ 0 0 ] + dindex, [ Q1, Q3 ],'-','Color', color );
    plot( [ 0 0 ] + dindex, [ Q1, Q3 ],'o','MarkerSize',2,'MarkerFaceColor', color,'MarkerEdgeColor', color );
    plot( 0 + dindex, Q2,'o','MarkerSize',4,'MarkerFaceColor', color,'MarkerEdgeColor',color );
end

xlim ( [ 1 dindex ] + [ -0.5 +0.5 ] )
