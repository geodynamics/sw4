%
% RfileHeader
%
%    Read the header information from an rfile
function [a]=rfileheader( fname )

  machineformat='ieee-le';

  fd=fopen(fname,'r',machineformat);
  if fd ~= -1 
% Read header
    magic    = fread(fd,1,'int');
    printf('magic = %d\n', magic);

    prec    = fread(fd,1,'int');
    att    = fread(fd,1,'int');
    az = fread(fd,1,'double');
    lon0 = fread(fd,1,'double');
    lat0 = fread(fd,1,'double');
    mlen = fread(fd,1,'int');
   mercstr = fread(fd,[1 mlen],'uchar');
   nb = fread(fd,1,'int');
% Display header
%   if verbose == 1
    printf('magic = %d, prec = %d, att = %d, az = %g, lon0 = %g, lat0 = %g\n', magic, prec, att, az, lon0, lat0);
    printf('mlen = %d, nb = %d\n', mlen, nb);
    printf('mercstr = %s\n', mercstr);
    printf('\n');
%   end;
   for p=1:nb
      hh(p) = fread(fd,1,'double');
      hv(p) = fread(fd,1,'double');
      z0(p) = fread(fd,1,'double');
      nc(p) = fread(fd,1,'int');
      ni(p) = fread(fd,1,'int');
      nj(p) = fread(fd,1,'int');
      nk(p) = fread(fd,1,'int');
%      if verbose == 1
%      disp(['    patch nr ' num2str(p) ' has h = ' num2str(h(p)) ' zmin = ' num2str(zmin(p))]);
      printf('patch = %d, hh=%g, hv=%g, z0=%g, nc=%d, ni=%d, nj=%d, nk=%d\n', p, hh(p), hv(p), z0(p), nc(p), ni(p), nj(p), nk(p));
%      end;
   end;
    fclose(fd);
else
   disp(['Error: could not open file ' imfile ]);
end;
