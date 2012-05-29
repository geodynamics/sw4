%
% READSAC
%
%   Read SAC receiever data.
%
%     [u,dt,lat,lon,t0] = readsac( fname, format )
%
%          Input: fname  - Name of SAC file
%                 format - Little endian ('l') or big endian ('b')
%                          byte order for binary data. Default is 'l'.
%
%          Output: u        - The data component on SAC file
%                  dt       - Uniform time step for u
%                  lat, lon - Latitude and longitude position of
%                             the SAC station.
%
   function [u,dt,lat,lon,t0] = readsac( fname, format )
  if nargin < 2 
     format = 'l';
  end;

  fid = fopen(fname,'r',format);
  if fid < 0 
    disp( ['Error: could not open file ' fname] );
  else
     dt = fread( fid,1,'float32');    
     fseek(fid,4*4,0);
     t0 = fread( fid,1,'float32');    
     t1 = fread( fid,1,'float32');    
     fseek(fid,24*4,0);
     lat = fread(fid,1,'float32');
     lon = fread(fid,1,'float32');
     fseek(fid,2*4,0);
     evlat = fread(fid,1,'float32');
     evlon = fread(fid,1,'float32');
     fseek(fid,4,0);
     evdepth = fread(fid,1,'float32');
     fseek(fid,4*31,0);
% integers from offset 70
     fseek(fid,4*6,0);
     nvhdr = fread(fid,1,'int32');
     fseek(fid,4*2,0);
     npts=fread(fid,1,'int32');
     fseek(fid,78*4,0);
% output required header data
disp(['Begin time (B) = ' num2str(t0) ' End time (E) = ' num2str(t1) ' Station lat lon = ' num2str(lat) ' ' num2str(lon) ' nvhdr = ' num2str(nvhdr) ' npts = ' num2str(npts)]);
% read time series
     u=fread(fid,npts,'float32');
     fclose(fid);
%    dt=fread(fid,1,'float32')
%    fseek(fid,78*4,0);
%    npts=fread(fid,1,'int')
%    fseek(fid,78*4,0);
%    u=fread(fid,npts,'float32');
%    fclose(fid);
  end
