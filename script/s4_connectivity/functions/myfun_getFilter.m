function map = myfun_getFilter ( source, template )

% Keeps only the sources with a valid area identifier.
areas        = template.grid.inside & template.grid.area > 0;
filters      = source.filter ( areas );
sourcearea   = template.grid.area ( areas );
sourcepos    = template.grid.pos  ( areas, : );

fprintf ( 1, '    Sources with one direction: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 1 ) );
fprintf ( 1, '    Sources with two directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 2 ) );
fprintf ( 1, '    Sources with three directions: %i.\n', sum ( cellfun ( @(x) size ( x, 1 ), filters ) == 3 ) );


% Lists the number of orientations for source position.
nori         = cellfun ( @(x) size ( x, 1 ), filters );

% Sets the source information in cell form.
nori         = num2cell ( nori );
sourcearea   = num2cell ( sourcearea );
sourcepos    = num2cell ( sourcepos, 2 );

% Repeats the area and position for each source.
sourcearea   = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcearea, nori, 'UniformOutput', false );
sourcepos    = cellfun ( @(x,y) repmat ( x, y, 1 ), sourcepos,  nori, 'UniformOutput', false );

% Goes back to vector form.
sourcearea   = cat ( 1, sourcearea {:} );
sourcepos    = cat ( 1, sourcepos  {:} );

% Gets the weight for each source.
sourcepow    = cellfun ( @diag, source.cov ( areas ), 'UniformOutput', false );
sourcepow    = cat ( 1, sourcepow {:} );

% Writes the beamformer filter in matrix form.
filter       = cat ( 1, filters {:} );


% Calculates the extended diagonal of the data.
diags        = cellfun ( @nan, nori, 'UniformOutput', false );
diags        = blkdiag ( diags {:} );

% Gets the original map.
map          = [];
map.label    = 'all';
map.filter   = filter;
map.pos      = sourcepos;
map.area     = sourcearea;
map.weight   = sourcepow;
map.diags    = diags;
