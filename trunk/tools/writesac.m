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
  npts = length(u)
  fdum = -12345;
  idum = -12345;
  cdum='blah';
  iftype=1;
  fid = fopen(fname,'wb',format);
  if fid < 0 
    disp( ['Error: could not open file ' fname] );
  else
% records 1-70 are float32
% delta
     fwrite( fid,dt,'float32');    

%     fseek(fid,4*4,0);
     for i=1:4
       fwrite(fid,fdum,'float32');
     end

% start time B
     fwrite( fid,t0,'float32');    
% end time E
     fwrite( fid,t0+npts*dt,'float32');    

%     fseek(fid,25*4,0);
     for i=1:24
       fwrite(fid,fdum,'float32');
     end

% station lat (32)
     fwrite(fid,lat,'float32');
% station lon (33)
     fwrite(fid,lon,'float32');

%     fseek(fid,2*4,0);
     for i=1:2
       fwrite(fid,fdum,'float32');
     end

% event lat (36)
     fwrite(fid,evlat,'float32');
% event lat (37)
     fwrite(fid,evlon,'float32');

%     fseek(fid,4,0);
     fwrite(fid,fdum,'float32');

% event depth (39)
     fwrite(fid,evdepth,'float32');

%     fseek(fid,4*40,0);
     for i=1:31
       fwrite(fid,fdum,'float32');
     end

% records 71-110 are int32
     for i=1:9
       fwrite(fid,idum,'int32');
     end

% npts (80)
     fwrite(fid,npts,'int32');

%     fseek(fid,78*4,0);

     for i=1:5
       fwrite(fid,idum,'int32');
     end

% iftype (86)
     fwrite(fid,iftype,'int32');


     for i=1:24
       fwrite(fid,idum,'int32');
     end

% records 111-158 are characters
     for i=1:48
       fwrite(fid,cdum,'uchar');
     end

% time series
     nu = fwrite(fid,u,'float32')

     fclose(fid);
  end
