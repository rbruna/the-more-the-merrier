function layout = my_prepare_layout ( data )

% If no FIFF file, relies on FieldTrip.
if ~ft_senstype ( data, 'neuromag' )
    
    % Initializes the configuration structure.
    cfg    = [];
    
    % Determines the appropriate layout from the data.
    switch ft_senstype ( data )
        case { 'eeg1020' 'eeg1010' 'eeg1005' }
            cfg.layout = upper ( ft_senstype ( data ) );
        otherwise
            error ( 'Unknown layout.' );
    end
    
    layout = ft_prepare_layout ( cfg );
    return
end


cfg.skipscale       = 'yes';
cfg.skipcomnt       = 'yes';
cfg.overlap         = 'keep';

if isfield ( data, 'elec' )
    cfg.elec            = data.elec;
    cfg.grad            = [];
    eleclayout          = ft_prepare_layout ( cfg );
else
    eleclayout          = [];
end

cfg.elec            = [];
cfg.grad            = data.grad;
gradlayout          = ft_prepare_layout ( cfg );

fields = fieldnames ( gradlayout );

if ~isempty ( eleclayout )
    for field = 1: numel ( fields )
        layout.( fields { field } ) = cat ( 1, eleclayout.( fields { field } ), gradlayout.( fields { field } ) );
    end
else
    layout = gradlayout;
end

% Corrects the sensors without position.
layout.pos ( isnan ( layout.pos ) ) = 0;
