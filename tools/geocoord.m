function [lat,lon]=geocoord(x, y, geoAz, latOrigin, lonOrigin)
% conversion factor between degrees and radians
  deg2rad = pi/180.0;
  metersPerDegree = 111319.5;
  phi = geoAz * deg2rad;
% Compute the latitude
  lat = latOrigin + (x*cos(phi) - y*sin(phi))/metersPerDegree;
% Compute the longitude
  lon = lonOrigin + (x*sin(phi) + y*cos(phi))/(metersPerDegree*cos(lat*deg2rad));
end
