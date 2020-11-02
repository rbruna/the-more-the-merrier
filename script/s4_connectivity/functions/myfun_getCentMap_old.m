function map = myfun_getCentMap ( filter, sourcepos, atlas )

% Gets the number of areas.
nareas      = numel ( atlas.name );

% Calculates the source nearest to the centroid of each area.
sourcedist  = sqrt ( sum ( bsxfun ( @minus, reshape ( sourcepos, [], 1, 3 ), reshape ( atlas.pos, 1, [], 3 ) ) .^ 2, 3 ));
[ ~, cent ] = min ( sourcedist, [], 1 );

% Gets the centroid beamformer filter.
filter      = filter ( cent, : );

% Generates the map for the centroid.
map         = [];
map.label   = 'centroid';
map.filter  = filter;
map.area    = ( 1: nareas )';
map.weight  = ones ( nareas, 1 );
map.diags   = diag ( nan ( nareas, 1 ) );
