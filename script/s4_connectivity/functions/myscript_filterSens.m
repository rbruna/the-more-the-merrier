
% Selects only the channels present in the leadfield.
cfg          = [];
cfg.channel  = source.label;
cfg.trials   = 'all';

banddata     = ft_selectdata ( cfg, sensdata.trialdata );

% Gets the order of the channels in the beamformer filter.
chanorder    = my_matchstr ( banddata.label, source.label );

% Extracts and sorts the channel data.
banddata     = cat ( 3, banddata.trial {:} );
banddata     = banddata ( chanorder, :, : );

% Filters the data in the selected band.
fir          = fir1 ( band.order, band.edges / ( band.fs / 2 ) );
banddata     = permute ( banddata, [ 2 1 3 ] );
banddata     = my_filtfilt ( fir, 1, banddata, true );
banddata     = permute ( banddata, [ 2 1 3 ] );

% Removes the padding.
banddata     = banddata ( :, padding + 1: end - padding, : );

% Downsamples the data, if required.
if config.downsample
    banddata     = banddata ( :, floor ( config.downsample / 2 ): config.downsample: end, : );
end
