function my_maximize ( handles )

% Gets the handle objects.
handles = handle ( handles );

% Operates for each object.
for hindex = 1: numel ( handles )
    
    % Gets the Java frame.
    frame = get ( handles ( hindex ), 'JavaFrame' );
    
    % Maximizes the frame.
    frame.setMaximized ( true );
end