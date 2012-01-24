%
% PLOTSHAKEMAP
%
%  Interpolate a value from the solution array u to location (lon, lat)
%
%              [smLon smLat smHVM]=plotshakemap( hvm, plotit )
%
%       Input: hvm - 2D array holding max horizontal velocities
%              plotit - 0: no plotting
%
%       Output: smLon - 1-D array of Longitude coordinates
%               smLat - 1-D array of Latitude coordinates
%               smHVM - 2-D array of interpolated values from input array hvm
%
function [smLon,smLat,smHVM]=plotshakemap( hvm, plotit )
smLon0 = -123.53;
smLat0 = 36.51333;
smLon = smLon0 + (0:180)/60;
smLat = smLat0 + (0:150)/60;
smHVM=zeros(181,151);
for i=1:181
 for j=1:151
  smHVM(i,j)=shakemap(hvm,smLon(i),smLat(j));
 end
end
if (plotit)
  contour(smLon, smLat, smHVM'); colorbar
end
