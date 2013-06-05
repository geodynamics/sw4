%
% READIMAGEPATCH3D
%
%    Read image produced by sw4 on format with the *.3Dimg extension.
%    This script will not read the older WPP image3d format.
%
%    NOTE: currently, no curvilinear grid coordinates are assigned to z
%
%         [im,x,y,z,t]=readimagepatch3d( imfile, pnr, verbose )
%
%                Input: imfile  - Image file name
%                       pnr     - Patch number, if more than one grid is used.
%                       verbose - Set to 1 to display file header information.
%                Output: im      - The image, as a 3D array.
%                        x, y, z - The spatial coordinates of the image, as 1D arrays
%                        t       - Simulation time at which the 3D image was output.

function [im3d,x,y,z,time]=readimagepatch3d( file, pnr, verbose )

if nargin < 3
   verbose=0;
end;

if nargin < 2
   pnr = 1;
end;

fd=fopen(file,'r');

% Read header
   prec    =fread(fd,1,'int');
   npatches=fread(fd,1,'int');
   time       =fread(fd,1,'double');
   plane   =fread(fd,1,'int');
   coord   =fread(fd,1,'double');
   mode    =fread(fd,1,'int');
   gridinfo=fread(fd,1,'int');
   timecreated=fread(fd,[1 25],'uchar');
   timestring=num2str(timecreated,'%c');
   mstr=getimage3dmode(mode);
% Display header
   if verbose == 1
      disp(['Found:  prec  = ' num2str(prec)  ' time = ' num2str(time) ' plane = ' num2str(plane)]);   
      disp(['        coord = ' num2str(coord) ' mode = ' mstr ' npatches = ' num2str(npatches) ]);
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

if pnr < 1 | pnr > npatches 
   disp(['Error in input patch nr']);
   disp(['File contains ' num2str(npatches) ' patches']);
else
   offset = 0;
   for p=1:pnr-1
     nvals = ni(p)*nj(p)*nk(p);
     offset  = offset + nvals;
   end;
   if( prec == 4 )
     fseek(fd,offset*4,'cof');
   else
     fseek(fd,offset*8,'cof');
   end;
% now read the patch that we are interested in
   npts = ni(pnr)*nj(pnr)*nk(pnr);
   if prec == 4
      im3 = fread(fd,npts,'float');
   else
      im3 = fread(fd,npts,'double');
   end;
   im3d = reshape(im3, ni(pnr), nj(pnr), nk(pnr));
   n1 = ni(pnr);
   n2 = nj(pnr);
   n3 = nk(pnr);
   x = [0:n1-1]*h(pnr);
   y = [0:n2-1]*h(pnr);
   z = [0:n3-1]*h(pnr)+zmin(pnr);
end;
