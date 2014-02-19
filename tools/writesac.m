%
% WRITESAC
%
%   Write SAC receiever data.
%
%     writesac( u, dt, lat, lon, t0, fname, format )
%
%          Input: 
%                 u        - The data component on SAC file
%                 dt       - Uniform time step for u
%                 lat, lon - Latitude and longitude position of
%                            the SAC station.
%                 t0       - start time
%                 fname  - Name of SAC file
%                 format - Little endian ('l') or big endian ('b')
%                          byte order for binary data. Default is 'l'.
%
%
   function [] = writesac(u, dt, lat, lon, t0, fname, format )
  if nargin < 7 
     format = 'l';
  end;

  evlat = 0;
  evlon = 0;
  evdepth = 0;
% required fields in the header
  npts = length(u);
  nvhdr = 6;
  b=t0;
  e=t0+npts*dt;
  iftype = 1;
  leven = 1;
  delta = dt;
% dummy values
  fdum = -12345.0;
  idum = -12345;
  cdum='-12345..';
  iftype=1;
  fid = fopen(fname,'wb',format);
  if fid < 0 
    disp( ['Error: could not open file ' fname] );
  else
% records 1-70 are float32
%
% DELTA (offset 0)
     fwrite( fid,dt,'float32');    

%     fseek(fid,4*4,0);
     for i=1:4
       fwrite(fid,fdum,'float32');
     end

% B (5)
     fwrite( fid,t0,'float32');    
% E (6)
     fwrite( fid,t0+(npts-1)*dt,'float32');    

%     fseek(fid,25*4,0);
     for i=1:24
       fwrite(fid,fdum,'float32');
     end

% station lat (31)
     fwrite(fid,lat,'float32');
% station lon (32)
     fwrite(fid,lon,'float32');

%     fseek(fid,2*4,0);
     for i=1:2
       fwrite(fid,fdum,'float32');
     end

% event lat (35)
     fwrite(fid,evlat,'float32');
% event lat (36)
     fwrite(fid,evlon,'float32');

%     fseek(fid,4,0);
     fwrite(fid,fdum,'float32');

% event depth (38)
     fwrite(fid,evdepth,'float32');

%     fseek(fid,4*40,0);
     for i=1:31
       fwrite(fid,fdum,'float32');
     end
% (70)
% records 71-110 are int32
     for i=1:6
       fwrite(fid,idum,'int32');
     end
% nvhdr (76)
     fwrite(fid,nvhdr,'int32');
% (77)
     fwrite(fid,idum,'int32');
% (78)
     fwrite(fid,idum,'int32');
% npts (79)
     fwrite(fid,npts,'int32');

%     fseek(fid,78*4,0);

     for i=1:5
       fwrite(fid,idum,'int32');
     end

% iftype (85)
     fwrite(fid,iftype,'int32');
% (86)
     for i=1:19
       fwrite(fid,idum,'int32');
     end
% LEVEN (105)
     fwrite(fid,leven,'int32');
% (106)
     for i=1:4
       fwrite(fid,idum,'int32');
     end
% (110)
% offset 110-157 are characters
% the size of 'uchar' is machine dependent, but we must have 32 bit records
     for i=1:48
       fwrite(fid,cdum,'int32');
     end
%     for i=1:48
%       fwrite(fid,idum,'int32');
%     end

% time series (158)
     nu = fwrite(fid,u,'float32')

     fclose(fid);
  end
