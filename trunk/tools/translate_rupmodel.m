% ----------------------------------------------------------------------
% Read fields from HDF5 file and outputs text file named rupmodel.txt
%
%                 translate_rupmodel( ruptureFile )
%
%         Input: ruptureFile - HDF5 file describing the rupture
% ----------------------------------------------------------------------
   function translate_rupmodel( ruptureFile )

  % Filename for rupture model HDF5 file (NOW READ FROM ARGUMENT)
%filename = 'hs_r01_hypoF_vr92_tr15.h5';
%filename = 'hs_r01_hypoO_vr92_tr15.h5';
%filename = 'hs_r01_hypoH_vr92_tr15.h5';
%filename = 'oakland.h5';
%filename = 'alumrock.h5';

% Final slip (left-lateral, reverse, opening)
finalSlip = hdf5read(ruptureFile, ...
                     '/rupture_model/slip_summary/final_slip/HaywardCalaveras/');

% Rise time (t95)
t95 = hdf5read(ruptureFile, ...
               '/rupture_model/slip_summary/slip_time_fn/rise_time/HaywardCalaveras');

% Slip time (t0, t1)
% t0 = time when slip starts at location
% t1 = t0 + t95 [we will use t95 from above rather than use t95=t1-t0]
slipTime = hdf5read(ruptureFile, ...
               '/rupture_model/slip_summary/slip_time/HaywardCalaveras');

% Get orientation of subfault (strike, dip, rake)
% Rake values are equal to 0.0 in orientation; rake values should be
% calculated from slip vector.
orientation = hdf5read(ruptureFile, ...
                       '/rupture_model/slip_summary/slip_orientation/HaywardCalaveras');

% Coordinates of subfault centroids (longitude, latitude, elevation)
coordinates = hdf5read(ruptureFile, ...
                       '/rupture_model/auxiliary/coordinates/HaywardCalaveras');

% Area of subfaults
subfaultArea = hdf5read(ruptureFile, ...
                        '/rupture_model/auxiliary/subfault_area/HaywardCalaveras');


% ----------------------------------------------------------------------
% Read shear modulus from HDF5 file.
% ----------------------------------------------------------------------

% Filename for rupture model HDF5 file (NOW READ FROM ARGUMENT)
%filename = 'shearmodulus.h5';

% shear modulus at subfaults
%shearmod = hdf5read(shearmodFile, '/rupture_model/auxiliary/shear_modulus/HaywardCalaveras/');

% ----------------------------------------------------------------------
% Compute additional fields as necessary.
% ----------------------------------------------------------------------

% Compute rake angle
rakeAngle = 180.0/pi*atan2(finalSlip(2,:), finalSlip(1,:));

% Magnitude of slip
slipMag = sqrt(finalSlip(1,:).^2 + finalSlip(2,:).^2);

% the shearmod field contains incorrect values
%moment = shearmod' .* subfaultArea' .* slipMag;
moment = subfaultArea' .* slipMag;


% ----------------------------------------------------------------------
% Dump fields to an ASCII file.
% ----------------------------------------------------------------------

data = [coordinates(1,:); coordinates(2,:); coordinates(3,:); ...
       orientation(1,:); orientation(2,:); rakeAngle; ...
       moment; t95'; slipTime(1,:)];

fout = fopen('rupmodel.txt', 'w');

fprintf(fout, '# Rupture model for Hayward scenarios\n');
fprintf(fout, '#\n');
fprintf(fout, '# Columns are:\n');
fprintf(fout, '# (1) Longitude (degrees, WGS84)\n');
fprintf(fout, '# (2) Latitude (degrees, WGS84)\n');
fprintf(fout, '# (3) Elevation (m)\n');
fprintf(fout, '# (4) Strike (degrees)\n');
fprintf(fout, '# (5) Dip (degrees)\n');
fprintf(fout, '# (6) Rake (degrees)\n');
fprintf(fout, '# (7) Moment/shearmod = subFaultArea*slipMag (m^3)\n');
fprintf(fout, '# (8) Rise time, t95 (s)\n');
fprintf(fout, '# (9) Slip (start) time (s)\n');

fprintf(fout, '%10.5f %10.5f %8.1f %6.1f %6.1f %6.1f %12.6e %5.2f %7.3f\n', data);

fclose(fout);
