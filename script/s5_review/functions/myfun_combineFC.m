function map = myfun_combineFC ( map )

% Gets the map information.
area          = map.area;
weights       = map.weight;
diags         = map.diags;


% Gets the connectivity metric.
if isfield ( map, 'chcoh' )
    metric        = 'hcoh';
    cconn         = map.chcoh;
elseif isfield ( map, 'cplv' )
    metric        = 'plv';
    cconn         = map.cplv;
end

% Gets the size of the problem.
nareas        = max ( area );
ntrials       = size ( cconn, 3 );

% Sets the all-sources diagonal to NaN.
cconn         = bsxfun ( @plus, cconn, diags );

% Reserves memory for the inter-area Hilbert coherence matrices.
conn_avg      = zeros ( nareas, nareas, ntrials, 'single' );
ciconn_avg    = zeros ( nareas, nareas, ntrials, 'single' );

conn_rms      = zeros ( nareas, nareas, ntrials, 'single' );
ciconn_rms    = zeros ( nareas, nareas, ntrials, 'single' );

% Goes through each area 'i'.
for aindex1 = 1: nareas
    
    % Gets the sources in this area and their power.
    sindex1       = ( area == aindex1 );
    weight1       = weights ( sindex1 );
    
    % Gets the intra-area sub-matrix.
    cconn_ii      = cconn ( sindex1, sindex1, : );
    cconn_ii      = reshape ( cconn_ii, [], ntrials );
    
    
    % Gets the combined power for each source pair.
    mweight       = sqrt ( weight1 * weight1' );
    mweight       = mweight (:);
    mweight       = repmat ( mweight, 1, ntrials );
    
    % Sets the required pairs to NaN.
    mweight ( isnan ( cconn_ii ) ) = nan;
    
    % Normalizes the weights.
    mweight       = bsxfun ( @times, mweight, sum ( isfinite ( mweight ), 1 ) ./ nansum ( mweight, 1 ) );
    
    
    % Calculates the average intra-area connectivity.
    conn_avg_ii   = nanmean ( abs ( cconn_ii ) .* mweight, 1 );
    iconn_avg_ii  = nanmean ( abs ( imag ( cconn_ii ) ) .* mweight, 1 );
    rconn_avg_ii  = nanmean ( abs ( real ( cconn_ii ) ) .* mweight, 1 );
    ciconn_avg_ii = iconn_avg_ii ./ sqrt ( 1 - rconn_avg_ii .^ 2 );
    
    % Calculates the RMS intra-area connectivity.
    conn_rms_ii   = sqrt ( nanmean ( abs ( cconn_ii ) .^ 2 .* mweight, 1 ) );
    iconn_rms_ii  = sqrt ( nanmean ( abs ( imag ( cconn_ii ) ) .^ 2 .* mweight, 1 ) );
    rconn_rms_ii  = sqrt ( nanmean ( abs ( real ( cconn_ii ) ) .^ 2 .* mweight, 1 ) );
    ciconn_rms_ii = iconn_rms_ii ./ sqrt ( 1 - rconn_rms_ii .^ 2 );
    
    % Stores the average PLV value.
    conn_avg   ( aindex1, aindex1, : ) = conn_avg_ii;
    ciconn_avg ( aindex1, aindex1, : ) = ciconn_avg_ii;
    
    % Stores the RMS PLV value.
    conn_rms   ( aindex1, aindex1, : ) = conn_rms_ii;
    ciconn_rms ( aindex1, aindex1, : ) = ciconn_rms_ii;
    
    % Goes through each area 'j'.
    for aindex2 = aindex1 + 1: nareas
        
        % Gets the sources in this area and their power.
        sindex2       = ( area == aindex2 );
        weight2       = weights ( sindex2 );
        
        % Gets the intra-area sub-matrix.
        cconn_ij      = cconn ( sindex1, sindex2, : );
        cconn_ij      = reshape ( cconn_ij, [], ntrials );
        
        
        % Gets the combined power for each source pair.
        mweight       = sqrt ( weight1 * weight2' );
        mweight       = mweight (:);
        mweight       = repmat ( mweight, 1, ntrials );
        
        % Sets the required pairs to NaN.
        mweight ( isnan ( cconn_ij ) ) = nan;
        
        % Normalizes the weights.
        mweight       = bsxfun ( @times, mweight, sum ( isfinite ( mweight ), 1 ) ./ nansum ( mweight, 1 ) );
        
        
        % Calculates the average inter-area connectivity.
        conn_avg_ij   = nanmean ( abs ( cconn_ij ) .* mweight, 1 );
        iconn_avg_ij  = nanmean ( abs ( imag ( cconn_ij ) ) .* mweight, 1 );
        rconn_avg_ij  = nanmean ( abs ( real ( cconn_ij ) ) .* mweight, 1 );
        ciconn_avg_ij = iconn_avg_ij ./ sqrt ( 1 - rconn_avg_ij .^ 2 );
        
        % Calculates the RMS inter-area connectivity.
        conn_rms_ij   = sqrt ( nanmean ( abs ( cconn_ij ) .^ 2 .* mweight, 1 ) );
        iconn_rms_ij  = sqrt ( nanmean ( abs ( imag ( cconn_ij ) ) .^ 2 .* mweight, 1 ) );
        rconn_rms_ij  = sqrt ( nanmean ( abs ( real ( cconn_ij ) ) .^ 2 .* mweight, 1 ) );
        ciconn_rms_ij = iconn_rms_ij ./ sqrt ( 1 - rconn_rms_ij .^ 2 );
        
        % Stores the average PLV value.
        conn_avg   ( aindex1, aindex2, : ) = conn_avg_ij;
        conn_avg   ( aindex2, aindex1, : ) = conn_avg_ij;
        ciconn_avg ( aindex1, aindex2, : ) = ciconn_avg_ij;
        ciconn_avg ( aindex2, aindex1, : ) = ciconn_avg_ij;
        
        % Stores the RMS PLV value.
        conn_rms   ( aindex1, aindex2, : ) = conn_rms_ij;
        conn_rms   ( aindex2, aindex1, : ) = conn_rms_ij;
        ciconn_rms ( aindex1, aindex2, : ) = ciconn_rms_ij;
        ciconn_rms ( aindex2, aindex1, : ) = ciconn_rms_ij;
    end
end

% Stores the per-area FC values.
map.( sprintf ( '%s_avg', metric ) ) = conn_avg;
map.( sprintf ( '%s_rms', metric ) ) = conn_rms;
map.( sprintf ( 'ci%s_avg', metric ) ) = ciconn_avg;
map.( sprintf ( 'ci%s_rms', metric ) ) = ciconn_rms;
