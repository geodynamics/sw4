%
% READMAT
%
%         [arr,x,y,z,pstr]=readmat( imfile, pnr, verbose )
%
%                Input: imfile  - name of rfile
%                       pnr     - patch number 
%                       verbose - Set to 1 to display file header information.
%                Output: arr     - Data for patch 'pnr', as a 4d array.
%                        x, y, z - The spatial coordinates of the image as 1-D arrays.
%                        pstr    - String specifying the Proj4 projection command
function [arr,x,y,z,pstr]=readmat( imfile, pnr, verbose )
if nargin < 3
   verbose = 0;
end;
% read block 1 by default
if nargin < 2
   pnr = 1;
end;

fd=fopen(imfile,'r');
if fd ~= -1 
% Read header
   magic   =fread(fd,1,'int');
   prec    =fread(fd,1,'int');
   att     =fread(fd,1,'int');
   alpha   =fread(fd,1,'double');
   xlon0   =fread(fd,1,'double');
   ylat0   =fread(fd,1,'double');
   len     =fread(fd,1,'int');
   pstr    =fread(fd,[1 len],'uchar');
% parsing the string as an array of uchar
   pstr = num2str(pstr,'%c');
   nblocks  =fread(fd,1,'int');
% Display header
   if verbose == 1
      disp(['Header: magic  = ' num2str(magic) ' prec= ' num2str(prec) ' att= ' num2str(att)]);   
      disp(['        azimuth= ' num2str(alpha,14) ' lon0= ' num2str(xlon0,14) ' lat0= ' num2str(ylat0,14)] );
      disp(['        len = ' num2str(len) ' proj-string= "' num2str(pstr(1:len)), '"']);
      disp(['        nblocks= ' num2str(nblocks) ]);
   end;
   for p=1:nblocks
      hh(p) = fread(fd,1,'double'); % horizontal grid size
      hv(p) = fread(fd,1,'double'); % vertical grid size
      z0(p) = fread(fd,1,'double'); % base z-level (not used for topo)
      nc(p) = fread(fd,1,'int'); % number of components in this block
      ni(p) = fread(fd,1,'int');
      nj(p) = fread(fd,1,'int');
      nk(p) = fread(fd,1,'int');
      if verbose == 1
         disp(['  Block number ' num2str(p)]);
	 disp(['      hh = ' num2str(hh(p)) ' hv = ' num2str(hv(p)) ' z0 = ' num2str(z0(p))]);
	 disp(['      nc = ' num2str(nc(p)) ' ni = ' num2str(ni(p)) ' nj = ' num2str(nj(p)) ' nk = ' num2str(nk(p))]);
      end;
   end;

% Read data
  readz = 0;

% skip over previous blocks
  if pnr <= nblocks
     for p=1:pnr-1
	fseek(fd,ni(p)*nj(p)*nk(p)*nc(p)*prec,'cof');
     end;
  end;

  if prec == 4
    arr = fread(fd, ni(pnr)*nj(pnr)*nk(pnr)*nc(pnr),'float');
  else
    arr = fread(fd, ni(pnr)*nj(pnr)*nk(pnr)*nc(pnr),'double');
end;
%  arr = reshape(arr,nc(pnr),ni(pnr),nj(pnr),nk(pnr));
% "C"-order
  arr = reshape(arr,nc(pnr),nk(pnr),nj(pnr),ni(pnr));

% make grid arrays
  x = ((0:ni(pnr)-1)*hh(pnr))';
  y = ((0:nj(pnr)-1)*hh(pnr))';
  if nk(pnr)>1
    z = z0(pnr) + ((0:nk(pnr)-1)*hv(pnr))';
  else
    z = z0(pnr);
  end;

else
   disp(['Error: could not open file ' imfile ]);
end;
fclose(fd);
