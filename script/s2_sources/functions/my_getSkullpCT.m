function skullmask = my_getSkullpCT ( mri )

% Gets a skull estimation from the pseudo-CT.

% Gets the scalp mask.
scalppos   = strcmpi ( mri.masklabel, 'scalp' );
scalpmask  = mri.mask ( :, :, :, scalppos );

% Gets a estimation of the brain.
brainmask  = mri.white + mri.gray + mri.csf > .5;


% Defines the thresholds for the fisrt and second passes.
threshold1 = 200;
threshold2 = 100;

% Gets the skull.
skullmask  = mri.pct > threshold1;
% skullmask  = smooth3 ( mri.pct, 'gaussian', 3 ) > threshold1;

% Removes the points outside the scalp.
skullmask  = skullmask .* single ( scalpmask );

% Forces the minimum skull thickness.
skullmask  = skullmask | dilateC ( dilateC ( brainmask ) );


% Gets the biggest connected element.
skullmask  = bwlabeln ( skullmask, 6 );
skullmask  = skullmask == mode ( skullmask ( skullmask ~= 0 ) );

% Fixes the skull.
skullmask  = dilateO2 ( skullmask );
skullmask  = imfill   ( skullmask, 'holes' );
skullmask  = erodeO2  ( skullmask );


% Goes back to the original image.
skullmask  = skullmask & mri.pct > threshold2;

% Forces the minimum skull thickness.
skullmask  = skullmask | dilateC ( dilateC ( brainmask ) );

% Keeps only the biggest connected element.
skullmask  = bwlabeln ( skullmask, 6 );
skullmask  = ( skullmask == mode ( skullmask ( skullmask ~= 0 ) ) );

% Fixes the skull.
skullmask  = dilateC  ( skullmask );
skullmask  = imfill   ( skullmask, 'holes' );
skullmask  = erodeC   ( skullmask );

return
% Dilates the image.
skullmask  = dilateO2 ( skullmask );

% Goes back to the original image.
skullmask  = skullmask & mri.pct > threshold2;

% Forces the minimum skull thickness.
skullmask  = skullmask | dilateC ( dilateC ( brainmask ) );

% Keeps only the biggest connected element.
skullmask  = bwlabeln ( skullmask, 6 );
skullmask  = ( skullmask == mode ( skullmask ( skullmask ~= 0 ) ) );

% Fixes the skull.
skullmask  = dilateC  ( skullmask );
skullmask  = imfill   ( skullmask, 'holes' );
skullmask  = erodeC   ( skullmask );
