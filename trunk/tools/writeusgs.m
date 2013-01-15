%
% WRITEUSGS
%
%  Write time series data in format specified by USGS for the Hayward
%  fault earthquake scenarios.
%
%             writeusgs( filename, stname, ux, uy, uz, dt, t0, enz, vel, utc )
%
%       Input: filename - write to this file
%              stname   - Station name, (to be saved in the header of the file)
%              ux - First component of time series data
%              uy - Second component
%              uz - Third component
%              dt - Time step used for time series data
%              lat, lon - Latitude and longitude of event
%              t0 - time offset to some reference time
%              utc - Reference time coordinate, 7 integers, [year,month,day,hour,minute,second,millisecond]
%              enz - 1  [ux,uy,uz] is [East-West,North-South,Up]
%                    0  [ux,uy,uz] is [X, Y, Z ]
%              vel - Defaults to 1
%                    0 Displacements
%                    1 Velocities
%
function []= writeusgs( filename, stname, ux, uy, uz, dt, t0, enz, vel, utc )
%
% Default to xyz-velocities
if nargin < 9
   utc(1) = 2011;
   utc(2)=10;
   utc(3)=25;
   utc(4)=19;
   utc(5)=00;
   utc(6)=00;
   utc(7)=11;
end;
if nargin < 8
  vel = 1;
end;
if nargin < 8
  enz = 0;
end;
%
fd=fopen(filename,'w');
fprintf(fd,'# Author: Matlab/Octave\n');
fprintf(fd,'# Scenario: test\n');
%fprintf(fd,['# Date: ' date '\n']);
%fprintf(fd,'# Date: UTC 10/25/2011:19:00:00.011\n');
fprintf(fd,'# Date: UTC %02d/%02d/%04d:%02d:%02d:%02d.%03d\n',utc(2),utc(3),utc(1),utc(4),utc(5),utc(6),utc(7));
fprintf(fd,'# Bandwidth: x Hz\n');
fprintf(fd,['# Station: ' stname '\n']);
fprintf(fd, '# Target location (WGS84 longitude, latitude) (deg): xx yy\n');
fprintf(fd, '# Actual location (WGS84 longitude, latitude) (deg): xx yy\n');
fprintf(fd, '# Distance from target to actual location (m): 0\n');
fprintf(fd, '# nColumns: 4\n');
fprintf(fd, '# Column 1: Time (s)\n');
if enz == 1 && vel == 1
   fprintf(fd, '# Column 2: East-west velocity (m/s)\n');
   fprintf(fd, '# Column 3: North-south velocity (m/s)\n');
   fprintf(fd, '# Column 4: Up-down velocity (m/s)\n');
elseif enz == 0 && vel == 0
   fprintf(fd, '# Column 2: X displacement (m)\n');
   fprintf(fd, '# Column 3: Y displacement (m)\n');
   fprintf(fd, '# Column 4: Z displacement (m)\n');
elseif enz == 1 && vel == 0
   fprintf(fd, '# Column 2: East-west displacement (m)\n');
   fprintf(fd, '# Column 3: North-south displacement (m)\n');
   fprintf(fd, '# Column 4: Up-down displacement (m)\n');
else
   fprintf(fd, '# Column 2: X velocity (m/s)\n');
   fprintf(fd, '# Column 3: Y velocity (m/s)\n');
   fprintf(fd, '# Column 4: Z velocity (m/s)\n');
end;
n = length(ux);
t = t0;
for i = 1:n
   fprintf(fd, ' %.12e  %.12e  %.12e  %.12e \n',t,ux(i),uy(i),uz(i));
   t = t+ dt;
end;
fclose(fd);

