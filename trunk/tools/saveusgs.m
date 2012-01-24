%
% SAVEUSGS
%
%  Save receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios.
%
%              saveusgs(t, ux, uy, uz, fc, filename)
%
%       Input: 
%               t  - Time vector of the same length as the data vectors (1st column of data).
%               ux - X- or East-West direction data component (2nd column of data )
%               uy - Y- or North-South direction data component (3rd column of data )
%               uz - Z- or Up direction data component (4th column of data)
%               fc - corner frequency = band width [Hz]
%               filename - Name of receiever data file
%               
%       Note: [ux,uy,uz] is [East-West,North-South,Up] or [X,Y,Z] components
%       depending on how the data file was written by WPP. 
%
%
function saveusgs( t, ux, uy, uz, fc, fileName)
N = length(t);

% open file
fd = fopen(fileName, 'w');

% write header
fprintf(fd,'# Author: Anders Petersson\n');
fprintf(fd,'# Scenario: %s\n', 'TEST');
fprintf(fd,'# Date: %s\n', date());
fprintf(fd,'# Bandwidth (Hz): %15.7e\n', fc);
fprintf(fd,'# Station: XXYY\n');
fprintf(fd,'# Target location XX YY\n');
fprintf(fd,'# Actual location XX YY\n');
fprintf(fd,'# Distance from target to actual location [m] XY\n');
fprintf(fd,'# Column 1: Time [s]\n');
fprintf(fd,'# Column 2: East-west velocity [m/s]\n');
fprintf(fd,'# Column 3: North-south velocity [m/s]\n');
fprintf(fd,'# Column 4: Up-down velocity [m/s]\n');
  
% output all the data
for i=1:N
  fprintf(fd,'%15.7e %15.7e %15.7e %15.7e\n', t(i), ux(i), uy(i), uz(i));
end

% close file
fclose(fd);


