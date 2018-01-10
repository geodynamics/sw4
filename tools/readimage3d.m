%
% READIMAGE3D
%
%    Read volumeimage produced by sw4 on format with the *.sw4img extension.
%    This script will not read the older WPP image format.
%
%         [im,x,y,z,t,timestring]=readimage3d( imfile, pnr, verbose )
%
%                Input: imfile  - Image file name
%                       pnr     - Patch number, if more than one grid is used.
%                       verbose - Set to 1 to display file header information.
%                Output: im      - The image, as a 2D array.
%                        x, y, z - The spatial coordinates of the image, one of these is a scalar.
%                        t       - Simulation time at which the image was output.
%                        timestring - String holding creation date
function [im,x,y,z,t,timestring]=readimage3d( imfile, pnr, verbose )
if nargin < 3
   verbose= 0;
end;
if nargin < 2
   pnr = 1;
end;

fd=fopen(imfile,'r');
if fd ~= -1 
% Read header
   prec    =fread(fd,1,'int');
   npatches=fread(fd,1,'int');
   t       =fread(fd,1,'double');
   plane   =fread(fd,1,'int');
   coord   =fread(fd,1,'double');
   mode    =fread(fd,1,'int');
   gridinfo=fread(fd,1,'int');
   timecreated=fread(fd,[1 25],'uchar');
   timestring=num2str(timecreated,'%c');
   mstr=getimage3dmode(mode);
% Display header
   if verbose == 1
      disp(['Found:  prec  = ' num2str(prec)  ' t= ' num2str(t) ' plane= ' num2str(plane)]);   
      disp(['        coord = ' num2str(coord) ' mode= ' mstr ' npatches= ' num2str(npatches) ]);
      disp(['        gridinfo = ' num2str(gridinfo) ]); 
     disp(['    file created ' timecreated(1:24)]);
   end;
   for p=1:npatches
      h(p) = fread(fd,1,'double');
      zmin(p) = fread(fd,1,'double');
      ib(p) = fread(fd,1,'int');
      ni(p) = fread(fd,1,'int');
      jb(p) = fread(fd,1,'int');
      nj(p) = fread(fd,1,'int');
      kb(p) = fread(fd,1,'int');
      nk(p) = fread(fd,1,'int');
      if verbose == 1
         disp(['    patch nr ' num2str(p) ' has h = ' num2str(h(p)) ' zmin = ' num2str(zmin(p))]);
      end;
   end;
% Read data
  readz = 0;
  if pnr <= npatches
     for p=1:pnr-1
	fseek(fd,(ni(p)-ib(p)+1)*(nj(p)-jb(p)+1)*(nk(p)-kb(p)+1)*prec,'cof');
     end;
     if prec == 4
        im0 = fread(fd,(ni(pnr)-ib(pnr)+1)*(nj(pnr)-jb(pnr)+1)*(nk(pnr)-kb(pnr)+1),'float');
     else
        im0 = fread(fd,(ni(pnr)-ib(pnr)+1)*(nj(pnr)-jb(pnr)+1)*(nk(pnr)-kb(pnr)+1),'double');
     end;
     im = reshape(im0,ni(pnr)-ib(pnr)+1,nj(pnr)-jb(pnr)+1,nk(pnr)-kb(pnr)+1);
% read grid z-coordinates
     if (pnr == npatches) && (gridinfo == 1 )
        if prec == 4
           z0 = fread(fd,(ni(pnr)-ib(pnr)+1)*(nj(pnr)-jb(pnr)+1)*(nk(pnr)-kb(pnr)+1),'float');
        else
           z0 = fread(fd,(ni(pnr)-ib(pnr)+1)*(nj(pnr)-jb(pnr)+1)*(nk(pnr)-kb(pnr)+1),'double');
        end;
	z = reshape(z0,ni(pnr)-ib(pnr)+1,nj(pnr)-jb(pnr)+1,nk(pnr)-kb(pnr)+1);
        readz = 1;
     end;
     x = (ib(pnr)-1:ni(pnr)-1)*h(pnr);
     y = (jb(pnr)-1:nj(pnr)-1)*h(pnr);
     if readz == 0
       z = (kb(pnr)-1:nk(pnr)-1)*h(pnr);
     end;
     fclose(fd);
  else
     disp(['Error: number of patches on file ' num2str(npatches) ' is smaller than input pnr = ' num2str(pnr)]);
  end;
else
   disp(['Error: could not open file ' imfile ]);
end;
