function maps = myfun_getPCAMap ( rawmap, banddata )

% Gets the original filter and the area labels.
filter     = rawmap.filter;
sourcearea = rawmap.area;

% Gets the size of the problem.
nsamples   = size ( banddata, 2 );
ntrials    = size ( banddata, 3 );
nareas     = max ( sourcearea );

% Defines the PCAs to extract.
pcainfos   = struct ( 'label', { 'pc', 'pc95' 'pc99' 'pc999' }, 'thres', { 1 0.95 0.99 0.999 } );

% Initializes the maps cell array.
maps       = cell ( numel ( pcainfos ), 1 );

% Goes through each PCA map.
for mindex = 1: numel ( pcainfos )
    
    % Initializes the map.
    map.label  = pcainfos ( mindex ).label;
    map.filter = cell ( nareas, 1 );
    map.pos    = nan;
    map.weight = cell ( nareas, 1 );
    map.area   = cell ( nareas, 1 );
    map.diags  = [];
    
    % Stores the map.
    maps { mindex } = map;
end


% Calculates the sensor-level covariance matrix.
sensccov   = real ( banddata ( :, : ) ) * real ( banddata ( :, : ) )' / ( nsamples * ntrials );

% Goes through each area.
for aindex = 1: nareas
    
    % Gets the source covariance.
    sourcecov     = filter ( sourcearea == aindex, : ) * sensccov * filter ( sourcearea == aindex, : )';
    
    % Applies PCA.
    [ pcacoeffs, pcapow ] = svd ( sourcecov );
    pcapow        = diag ( pcapow );
    pcaratio      = pcapow / sum ( pcapow );
    
    % Goes through each PCA subset.
    for mindex = 1: numel ( pcainfos )
        
        % Gets the corrent subset.
        pcainfo    = pcainfos ( mindex );
        map        = maps  { mindex };
        
        % Gets the requested components.
        if pcainfo.thres < 1
            pccomps    = 1: find ( cumsum ( pcaratio ) >= pcainfo.thres, 1, 'first' );
        else
            pccomps    = 1: pcainfo.thres;
        end
        
        % Gets the coefficients for the requested components.
        pccoeffs   = pcacoeffs ( :, pccomps );
        
        % Gets the filters, area labels and weights.
        pcfilter   = pccoeffs' * filter ( sourcearea == aindex, : );
        pcarea     = aindex * ones ( numel ( pccomps ), 1 );
        pcpow      = pcapow ( pccomps );
        
        
        % Stores the map information.
        map.filter { aindex } = pcfilter;
        map.weight { aindex } = pcpow;
        map.area   { aindex } = pcarea;
        
        maps { mindex } =  map;
    end
end

% Goes through each map.
for mindex = 1: numel ( pcainfos )
    
    % Gets the corrent subset.
    map        = maps  { mindex };
    
    % Concatenates the filters, area mappings and weights.
    map.filter = cat ( 1, map.filter {:} );
    map.weight = cat ( 1, map.weight {:} );
    map.area   = cat ( 1, map.area   {:} );
    
    % Adds the diagonal information.
    map.diags  = diag ( nan ( numel ( map.area ), 1 ) );
    
    % Stores the subset.
    maps { mindex } = map;
end

% Concatenates all the maps.
maps       = cat ( 1, maps {:} );
