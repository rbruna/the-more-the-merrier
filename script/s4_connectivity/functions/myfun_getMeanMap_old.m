function map = myfun_getMeanMap ( filter, sourcearea )

% Gets the size of the problem.
nchans       = size ( filter, 2 );
nareas       = max ( sourcearea );

% Reserves memory for the filter.
meanfilter   = zeros ( nareas, nchans );

% Calculates the average of each area.
for aindex = 1: nareas
    meanfilter ( aindex, : ) = mean ( filter ( sourcearea == aindex, : ) );
end

% Generates the map for the average time series.
map         = [];
map.label   = 'mean';
map.filter  = meanfilter;
map.area    = ( 1: nareas )';
map.weight  = ones ( nareas, 1 );
map.diags   = diag ( nan ( nareas, 1 ) );
