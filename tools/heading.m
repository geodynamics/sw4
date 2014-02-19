% 
% (initial) heading from (lon1,lat1) towards (lon2,lat2)
%
% [azdeg, az] = heading(lon1,lat1,lon2,lat2)
%
% Input locations in degrees
%
% Output heading in radians (az) and degrees (azdeg)
%
function [azdeg,az] = heading(lon1,lat1,lon2,lat2)
dlonr=(lon2-lon1)*pi/180;
la1r=lat1*pi/180;
l01r=lon1*pi/180;
la2r=lat2*pi/180;
lo2r=lon2*pi/180;
y=sin(dlonr*cos(la2r));
x=cos(la1r)*sin(la2r)-sin(la1r)*cos(la2r)*cos(dlonr);
az = atan2(y,x);
azdeg=az*180/pi;
if (azdeg < 0)
  azdeg += 360;
endif


