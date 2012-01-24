%
% SAVESHAKEMAP
%
%  Interpolate a value from the solution array u to location (lon, lat)
%
%              saveshakemap( smLon, smLat, smHVM, fileName, scenario )
%
%       Input: smLon - 1-D array of Longitude coordinates
%              smLat - 1-D array of Latitude coordinates
%              smHVM - 2-D array of values (max horizontal velocities)
%              fileName - Name of file to write
%
function saveshakemap( smLon, smLat, smHVM, fileName, scenario )
Nx=length(smLon);
Ny=length(smLat);
fd = fopen(fileName, 'w');
% Header
fprintf(fd,'# Author: Anders Petersson\n');
fprintf(fd,'# Scenario: %s\n', scenario);
fprintf(fd,'# Date: %s\n', date());
fprintf(fd,'# Bandwidth (Hz): 0.500000\n');
fprintf(fd,'# Column 1: WGS84 longitude (deg)\n');
fprintf(fd,'# Column 2: WGS84 latitude (deg)\n');
fprintf(fd,'# Column 3: Peak ground velocity (m/s)\n');
fprintf(fd,'# Column 4 (optional): Peak ground acceleration (m/s^2)\n');
fprintf(fd,'# Column 5 (optional): SA @ T=0.3 s (m/s^2)\n');
fprintf(fd,'# Column 6 (optional): SA @ T=1.0 s (m/s^2)\n');
fprintf(fd,'# Column 7 (optional): SA @ T=3.0 s (m/s^2)\n');
% output all the data
for j=Ny:-1:1
  for i=1:Nx
    fprintf(fd,'%15.7e %15.7e %15.7e\n', smLon(i), smLat(j), smHVM(i,j));
  end
end
% done!
fclose(fd);


