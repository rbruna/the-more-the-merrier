
% Defines the PCAs to extract.
pcainfos      = struct ( 'label', { 'pc', 'pc95' 'pc99' 'pc999' }, 'thres', { 1 0.95 0.99 0.999 } );

% Initializes the maps cell array.
pcamaps       = cell ( numel ( pcainfos ), 1 );

% Goes through each PCA map.
for mindex = 1: numel ( pcainfos )
    
    % Initializes the map.
    pcamap.label  = pcainfos ( mindex ).label;
    pcamap.filter = cell ( nareas, 1 );
    pcamap.weight = cell ( nareas, 1 );
    pcamap.area   = cell ( nareas, 1 );
    pcamap.diags  = [];
    
    % Stores the map.
    pcamaps { mindex } = pcamap;
end


% Calculates the sensor-level covariance matrix.
sensccov      = real ( banddata ( :, : ) ) * real ( banddata ( :, : ) )' / ( nsamples * ntrials );

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
        pcainfo       = pcainfos ( mindex );
        pcamap        = pcamaps  { mindex };
        
        % Gets the requested components.
        if pcainfo.thres < 1
            pccomps       = 1: find ( cumsum ( pcaratio ) >= pcainfo.thres, 1, 'first' );
        else
            pccomps       = 1: pcainfo.thres;
        end
        
        % Gets the coefficients for the requested components.
        pccoeffs      = pcacoeffs ( :, pccomps );
        
        % Gets the filters and area labels.
        pcfilter      = pccoeffs' * filter ( sourcearea == aindex, : );
        pcarea        = aindex * ones ( numel ( pccomps ), 1 );
        
        
        % Gets the weight for each component, if requested.
        if config.weight
            pcpow         = pcpow ( pccomps );
        else
            pcpow         = ones ( numel ( pccomps ), 1 );
        end
        
        
        % Stores the map information.
        pcamap.filter { aindex } = pcfilter;
        pcamap.weight { aindex } = pcpow;
        pcamap.area   { aindex } = pcarea;
        
        pcamaps { mindex } =  pcamap;
    end
end

% Goes through each map.
for mindex = 1: numel ( pcainfos )
    
    % Gets the corrent subset.
    pcamap        = pcamaps  { mindex };
    
    % Concatenates the filters, area mappings and weights.
    pcamap.filter = cat ( 1, pcamap.filter {:} );
    pcamap.weight = cat ( 1, pcamap.weight {:} );
    pcamap.area   = cat ( 1, pcamap.area   {:} );
    
    % Adds the diagonal information.
    pcamap.diags  = diag ( nan ( numel ( pcamap.area ), 1 ) );
    
    % Stores the subset.
    pcamaps { mindex } = pcamap;
end

% Concatenates all the maps.
pcamaps       = cat ( 1, pcamaps {:} );
