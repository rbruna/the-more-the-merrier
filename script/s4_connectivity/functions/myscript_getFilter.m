
% Gets the requested orientations for each source.
source       = my_selectOri ( source, config.threshold );

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

% Gets the weight for each source, if requested.
if config.weight
    sourcepow    = cellfun ( @diag, source.cov ( areas ), 'UniformOutput', false );
    sourcepow    = cat ( 1, sourcepow {:} );
else
    sourcepow    = ones ( size ( sourcearea ) );
end

% Writes the beamformer filter in matrix form.
filter       = cat ( 1, filters {:} );
