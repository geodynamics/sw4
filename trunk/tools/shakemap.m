%
% SHAKEMAP
%
%  Interpolate a value from the solution array u to location (lon, lat)
%
%              val = plotusgs(u, lon, lat)
%
%       Input: u - 2-D solution array
%              lon - longitude
%              lat - latitude
%       Output: val - (scalar) solution value at (lon,lat)
function [val]=shakemap( u, lon, lat )
mMetersPerDegree=111319.5;
[Ny Nx]=size(u);
% convert (lon,lat) to (x,y)

% big domain (120 by 200 km
gridLat = 38.8;
gridLon = -122.25;

% small domain 100 by 100 km
%gridLat = 38.0;
%gridLon = -121.8;

mGeoAz = 144;
  
deg2rad=pi/180;
phi = mGeoAz * deg2rad;

% compute x and y
x = mMetersPerDegree*(cos(phi)*(lat-gridLat) + cos(lat*deg2rad)*(lon-gridLon)*sin(phi));
y = mMetersPerDegree*(-sin(phi)*(lat-gridLat) + cos(lat*deg2rad)*(lon-gridLon)*cos(phi));

% grid size
h=100;
%h=200;

% grid indices
i = 1 + round(x/h);
j = 1 + round(y/h);

if (i<1 | i > Nx | j<1 | j > Ny)
  val = 0;
else
  val = u(j,i);
end

