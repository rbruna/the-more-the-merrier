function header = my_read_header ( filename )

% Tries to use the specific function.
if ft_filetype ( filename, 'neuromag_fif' )
    header = myfiff_read_header ( filename );
    
    % Tries to get the sensors information.
    [ grad, elec ]     = ft_mne2grad ( header.orig, false, [] );
    if ~isempty ( grad )
        header.grad        = grad;
    end
    if ~isempty ( elec )
        header.elec        = elec;
        header.elec.type   = 'eeg';
    end
    
elseif ft_filetype ( filename, 'egi_mff' )
    header = mymff_read_header ( filename );
    
elseif ft_filetype ( filename, 'ns_cnt' )
    header = mycnt_read_header ( filename );
    
% If no specific function relies on FieldTrip.
else
    header = ft_read_header ( filename );
end


% Checks the sensor definition.
if isfield ( header, 'grad' )
	header.grad        = ft_datatype_sens ( header.grad );
end
if isfield ( header, 'elec' )
	header.elec        = ft_datatype_sens ( header.elec );
end

% Extends the channel definitions.
if ~isfield ( header, 'chantype' )
    header.chantype    = ft_chantype ( header );
end
if ~isfield ( header, 'chanunit' )
    header.chanunit    = ft_chanunit ( header );
end

% Makes sure that all the vectors are column arrays.
if isfield ( header, 'label' )
    header.label       = header.label (:);
end
if isfield ( header, 'chantype' )
    header.chantype    = header.chantype (:);
end
if isfield ( header, 'chanunit' )
    header.chanunit    = header.chanunit (:);
end

% The metadata values must be in float precission.
header.Fs          = double ( header.Fs );
header.nSamples    = double ( header.nSamples );
header.nSamplesPre = double ( header.nSamplesPre );
header.nTrials     = double ( header.nTrials );
header.nChans      = double ( header.nChans );
