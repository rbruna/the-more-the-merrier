clc
clear
close all

% Starts with the ActiCap 64 layout.
load 'acticap-64ch-standard2'
lay = rmfield ( lay, 'cfg' );

% We need to add 4 electrodes:
% * Fpz between Fp1 and Fp2 (and slightly more frontal).
% * T9 to the right of T7.
% * T10 to the right of T8.
% * Iz dorsal to Oz.

% Gets the position of the 'reference' electrodes.
e_Fp1 = find ( strcmp ( lay.label, 'Fp1' ) );
e_Fp2 = find ( strcmp ( lay.label, 'Fp2' ) );
e_T7  = find ( strcmp ( lay.label, 'T7' ) );
e_T8  = find ( strcmp ( lay.label, 'T8' ) );
e_C5  = find ( strcmp ( lay.label, 'C5' ) );
e_C6  = find ( strcmp ( lay.label, 'C6' ) );
e_Oz  = find ( strcmp ( lay.label, 'Oz' ) );

p_Fp1 = lay.pos ( e_Fp1, : );
p_Fp2 = lay.pos ( e_Fp2, : );
p_T8  = lay.pos ( e_T8,  : );
p_T7  = lay.pos ( e_T7,  : );
p_C5  = lay.pos ( e_C5,  : );
p_C6  = lay.pos ( e_C6,  : );
p_Oz  = lay.pos ( e_Oz,  : );

% Generates the position of the new electrodes from the references.
p_Fpz = ( p_Fp1 + p_Fp2 ) / 2 + [ 0.000 0.010 ];
p_T10 = 2 * p_T8 - p_C6 - [ 0.005 0.000 ];
p_T9  = 2 * p_T7 - p_C5 + [ 0.005 0.000 ];
p_Iz  = p_Oz - [ 0.000 0.060 ];

% Adds the new electrodes.
lay.pos    = cat ( 1, lay.pos, p_Fpz );
lay.width  = cat ( 1, lay.width, lay.width (1) );
lay.height = cat ( 1, lay.height, lay.height (1) );
lay.label  = cat ( 1, lay.label, 'Fpz' );

lay.pos    = cat ( 1, lay.pos, p_T9 );
lay.width  = cat ( 1, lay.width, lay.width (1) );
lay.height = cat ( 1, lay.height, lay.height (1) );
lay.label  = cat ( 1, lay.label, 'T9' );

lay.pos    = cat ( 1, lay.pos, p_T10 );
lay.width  = cat ( 1, lay.width, lay.width (1) );
lay.height = cat ( 1, lay.height, lay.height (1) );
lay.label  = cat ( 1, lay.label, 'T10' );

lay.pos    = cat ( 1, lay.pos, p_Iz );
lay.width  = cat ( 1, lay.width, lay.width (1) );
lay.height = cat ( 1, lay.height, lay.height (1) );
lay.label  = cat ( 1, lay.label, 'Iz' );


% Gets the mapping to Elekta's naming system.
map = { ...
    'Fp1'     'EEG001';
    'Fpz'     'EEG002';
    'Fp2'     'EEG003';
    'AF7'     'EEG004';
    'AF3'     'EEG005';
    'AF4'     'EEG006';
    'AF8'     'EEG007';
    'F7'      'EEG008';
    'F5'      'EEG009';
    'F3'      'EEG010';
    'F1'      'EEG011';
    'Fz'      'EEG012';
    'F2'      'EEG013';
    'F4'      'EEG014';
    'F6'      'EEG015';
    'F8'      'EEG016';
    'FT9'     'EEG017';
    'FT7'     'EEG018';
    'FC5'     'EEG019';
    'FC1'     'EEG020';
    'FC2'     'EEG021';
    'FC6'     'EEG022';
    'FT8'     'EEG023';
    'FT10'    'EEG024';
    'T9'      'EEG025';
    'T7'      'EEG026';
    'C5'      'EEG027';
    'C3'      'EEG028';
    'C1'      'EEG029';
    'Cz'      'EEG030';
    'C2'      'EEG031';
    'C4'      'EEG032';
    'C6'      'EEG033';
    'T8'      'EEG034';
    'T10'     'EEG035';
    'TP9'     'EEG036';
    'TP7'     'EEG037';
    'CP3'     'EEG038';
    'CP1'     'EEG039';
    'CP2'     'EEG040';
    'CP4'     'EEG041';
    'TP8'     'EEG042';
    'TP10'    'EEG043';
    'P7'      'EEG044';
    'P5'      'EEG045';
    'P3'      'EEG046';
    'P1'      'EEG047';
    'Pz'      'EEG048';
    'P2'      'EEG049';
    'P4'      'EEG050';
    'P6'      'EEG051';
    'P8'      'EEG052';
    'PO7'     'EEG053';
    'PO3'     'EEG054';
    'PO4'     'EEG055';
    'PO8'     'EEG056';
    'O1'      'EEG057';
    'Oz'      'EEG058';
    'O2'      'EEG059';
    'Iz'      'EEG060' };

% Goes through each electrode.
for mindex = 1: size ( map, 1 )
    
    % Replaces the channel label.
    chindex = strcmp ( lay.label, map ( mindex, 1 ) );
    lay.label ( chindex ) = map ( mindex, 2 );
end

% Sorts the channels.
chindex = my_matchstr ( lay.label, map ( :, 2 ) );

lay.pos    = lay.pos    ( chindex, : );
lay.width  = lay.width  ( chindex );
lay.height = lay.height ( chindex );
lay.label  = lay.label  ( chindex );

% Saves the layout.
save -v6 neuromagEEG060 lay

% Plots the result.
ft_plot_lay ( lay  )
