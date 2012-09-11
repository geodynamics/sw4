%
% SACTOUSGS
%
%   Read three SAC files and combine them to one file in USGS format.
%
%        sactousgs( fname, sac1, sac2, sac3, sacdir )
%
%          Input: fname - File name of output USGS file
%                 sac1, sac2, sac3 - Three sac files representing recording of three components
%                                    at one single location
%                 sacdir - Directory of sac1,sac2,sac3, defaults to current directory
%          Output: The USGS file is written to 'fname'
%
function sactousgs( fname, sac1, sac2, sac3, sacdir )

if nargin <5
   sacdir = './';
end;
n=length(sacdir);
if sacdir(n) ~= '/'
   sacdir(n+1) = '/';
end;
convfactor=pi/180;

[u1,dt1,lat1,lon1,b1,e1,npts1,year1,jday1,hour1,min1,sec1,msec1,cmpaz1,cmpinc1,idep1,stnm1]=readsac([sacdir sac1]);
[u2,dt2,lat2,lon2,b2,e2,npts2,year2,jday2,hour2,min2,sec2,msec2,cmpaz2,cmpinc2,idep2,stnm2]=readsac([sacdir sac2]);
[u3,dt3,lat3,lon3,b3,e3,npts3,year3,jday3,hour3,min3,sec3,msec3,cmpaz3,cmpinc3,idep3,stnm3]=readsac([sacdir sac3]);

% Verify that origins of files agree
eflag = 0;
if (stnm1 ~= stnm2) | (stnm1 ~= stnm3) | (stnm2 ~= stnm3)
   disp(['Error sac stations names do not agree ' stnm1 ' ' stnm2 ' ' stnm3]);
   eflag = 1;
end;
if (lat1 ~= lat2) | (lat1 ~= lat3) | (lat2 ~= lat3)
   disp(['Error sac stations latitudes do not agree ' num2str(lat1) ' ' num2str(lat2) ' ' num2str(lat3)]);
   eflag = 1;
end;
if (lon1 ~= lon2) | (lon1 ~= lon3) | (lon2 ~= lon3)
   disp(['Error sac stations longitudes do not agree ' num2str(lon1) ' ' num2str(lon2) ' ' num2str(lon3)]);
   eflag = 1;
end;
if (dt1 ~= dt2) | (dt1~=dt3) | (dt2 ~= dt3)
   disp(['Error sac stations time steps do not agree ' num2str(dt1) ' ' num2str(dt2) ' ' num2str(dt3)]);
   eflag = 1;
end;			
if (npts1 ~= npts2) | (npts1~=npts3) | (npts2 ~= npts3)
   disp(['Error sac stations number of points do not agree ' num2str(npts1) ' ' num2str(npts2) ' ' num2str(npts3)]);
   eflag = 1;
end;			
if (b1 ~= b2) | (b1~=b3) | (b2 ~= b3)
   disp(['Error sac stations time offsets do not agree ' num2str(b1) ' ' num2str(b2) ' ' num2str(b3)]);
   eflag = 1;
end;			
if (year1 ~= year2) | (year1~=year3) | (year2 ~= year3)
   disp(['Error sac stations years do not agree ' num2str(year1) ' ' num2str(year2) ' ' num2str(year3)]);
   eflag = 1;
end;			
if (jday1 ~= jday2) | (jday1~=jday3) | (jday2 ~= jday3)
   disp(['Error sac stations jdays do not agree ' num2str(jday1) ' ' num2str(jday2) ' ' num2str(jday3)]);
   eflag = 1;
end;			
if (hour1 ~= hour2) | (hour1~=hour3) | (hour2 ~= hour3)
   disp(['Error sac stations hours do not agree ' num2str(hour1) ' ' num2str(hour2) ' ' num2str(hour3)]);
   eflag = 1;
end;			
if (min1 ~= min2) | (min1~=min3) | (min2 ~= min3)
   disp(['Error sac stations minutes do not agree ' num2str(min1) ' ' num2str(min2) ' ' num2str(min3)]);
   eflag = 1;
end;			
if (sec1 ~= sec2) | (sec1~= sec3) | (sec2 ~= sec3)
   disp(['Error sac stations seconds do not agree ' num2str(sec1) ' ' num2str(sec2) ' ' num2str(sec3)]);
   eflag = 1;
end;			
if (msec1 ~= msec2) | (msec1~=msec3) | (msec2 ~= msec3)
   disp(['Error sac station time offsets do not agree ' num2str(msec1) ' ' num2str(msec2) ' ' num2str(msec3)]);
   eflag = 1;
end;			
if eflag == 0
   cmpaz1  = convfactor*cmpaz1;
   cmpaz2  = convfactor*cmpaz2;
   cmpaz3  = convfactor*cmpaz3;
   cmpinc1 = convfactor*cmpinc1;
   cmpinc2 = convfactor*cmpinc2;
   cmpinc3 = convfactor*cmpinc3;
   if idep1 == 6  
      vel = 0;
   elseif idep1 == 7
      vel = 1;
   else
      disp(['Idep not recognized, using velocities']);
      vel = 1;
   end;
   [day,mon]= convertjday( jday1, year1 );
   utc(1) = year1;
   utc(2) = mon;
   utc(3) = day;
   utc(4) = hour1;
   utc(5) = min1;
   utc(6) = sec1;
   utc(7) = msec1;
   un = sin(cmpinc1)*cos(cmpaz1)*u1+sin(cmpinc2)*cos(cmpaz2)*u2 + sin(cmpinc3)*cos(cmpaz3)*u3;
   ue = sin(cmpinc1)*sin(cmpaz1)*u1+sin(cmpinc2)*sin(cmpaz2)*u2 + sin(cmpinc3)*sin(cmpaz3)*u3;
   uu = cos(cmpinc1)*u1 + cos(cmpinc2)*u2 + cos(cmpinc3)*u3;
   writeusgs( fname, stnm1, ue, un, uu, dt1, b1, 1, vel, utc );
else
   disp('Return witout creating USGS file');
end;

