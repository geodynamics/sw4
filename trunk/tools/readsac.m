%
% READSAC
%
%   Read SAC receiever data.
%
%     [u, dt, lat, lon, b, e, npts, year, jday, hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] = readsac( fname, format )
%
%          Input: fname  - Name of SAC file
%                 format - Little endian ('l') or big endian ('b')
%                          byte order for binary data. Default is 'l'.
%
%          Output: u        - The data component on SAC file
%                  dt       - Uniform time step for u
%                  stalat, stalon - Latitude and longitude position of
%                             the SAC station.
%                  b        - begin time relative reference datum
%                  e        - end time relative reference datum
%                  npts     - Number of elements in u
%                  year       - reference datum
%                  Julian day - reference datum
%                  hour       - reference datum
%                  minute     - reference datum
%                  second     - reference datum
%                  micro sec  - reference datum
%                  cmpaz      - Azimuth angle of component (degrees)
%                  cmpinc     - Inclination angle of compondent (degrees)
%                  idep       - Code of component stored (6-displacement, 7-velocity)
%                  stnam      - Name of the station
%
 function  [u, dt, lat, lon, b, e, npts, year, jday, hour, min, sec, msec, cmpaz, cmpinc, idep, stnam ] = readsac( fname, format )
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
     fseek(fid,4*18,0);
     cmpaz  = fread(fid,1,'float32');
     cmpinc = fread(fid,1,'float32');
     fseek(fid,4*11,0);
%     fseek(fid,4*31,0);
% integers from offset 70
     year = fread(fid,1,'int32');
     jday = fread(fid,1,'int32');
     hour = fread(fid,1,'int32');
     min = fread(fid,1,'int32');
     sec = fread(fid,1,'int32');
     msec = fread(fid,1,'int32');
     nvhdr = fread(fid,1,'int32');
     fseek(fid,4*2,0);
     npts=fread(fid,1,'int32');
     fseek(fid,6*4,0);
     idep = fread(fid,1,'int32');
     fseek(fid,23*4,0);
     stnam = fread(fid,8,'char');
     stnam = stnam';
     fseek(fid,37*4,0);
% output reference time stamp
disp(['Year = ' num2str(year) ' Julian Day = ' num2str(jday) ' Hour = ' num2str(hour) ' Min = ' num2str(min) ' Sec = ' num2str(sec) ' Micro Sec = ' num2str(msec) ]);
% output required header data
disp(['Begin time (B) = ' num2str(t0) ' End time (E) = ' num2str(t1) ' Station lat lon = ' num2str(lat) ' ' num2str(lon) ' nvhdr = ' num2str(nvhdr) ' npts = ' num2str(npts)]);
% read time series
     u=fread(fid,npts,'float32');
% copy begin and end times
     b = t0;
     e = t1;
     disp(['cmpaz  = ' num2str(cmpaz)]);
     disp(['cmpinc = ' num2str(cmpinc)]);
     disp(['idep    = ' num2str(idep)]);
     disp(['stnam   =  ' stnam ]);

     fclose(fid);
  end
