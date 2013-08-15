%
% READIMAGE
%
%    Read image produced by sw4 on format with the *.sw4img extension.
%    This script will not read the older WPP image format.
%
%         [im,x,y,z,plane,t,timestring]=readimage( imfile, pnr, verbose )
%
%                Input: imfile  - Image file name
%                       pnr     - Patch number, if more than one grid is used.
%                       verbose - Set to 1 to display file header information.
%                Output: im      - The image, as a 2D array.
%                        x, y, z - The spatial coordinates of the image, one of these is a scalar.
%                        plane   - 0: x=const, 1: y=const, 2: z=const
%                        t       - Simulation time at which the image was output.
%                        timestring - String holding creation date
function [im,x,y,z,plane,t,timestring]=readimage( imfile, pnr, verbose )
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
   mstr=getimagemodestr(mode);
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
      if verbose == 1
         disp(['    patch nr ' num2str(p) ' has h = ' num2str(h(p)) ' zmin = ' num2str(zmin(p))]);
      end;
   end;
% Read data
  readz = 0;
  if pnr <= npatches
     for p=1:pnr-1
	fseek(fd,(ni(p)-ib(p)+1)*(nj(p)-jb(p)+1)*prec,'cof');
     end;
     if prec == 4
        im0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'float');
     else
        im0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'double');
     end;
% read grid z-coordinates
     if (pnr == npatches) && (gridinfo == 1 )
%        disp(['Reading z-coordinates...']);   
        if prec == 4
           z0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'float');
        else
           z0 = fread(fd,[ni(pnr)-ib(pnr)+1 nj(pnr)-jb(pnr)+1],'double');
        end;
% transpose z0 so it gets the same dimensions as im
        z = z0';
        readz = 1;
     end;
% transpose im0 and return result in im
     im = im0';
     fclose(fd);
     if plane == 0 
        x = coord;
%        disp(['Making the grid: plane  = ' num2str(plane) ' readz = ' num2str(readz)]);   
        if readz == 0 
           y = h(pnr)*((ib(pnr):ni(pnr))-1);
           z = zmin(pnr) + h(pnr)*((jb(pnr):nj(pnr))-1);
        else
          n1 = size(im,1);
          n2 = size(im,2);
% make y two-dimensional (same size as z and im)
          y = zeros(n1,n2);
          for i=1:n1
            y(i,:) = h(pnr)*((ib(pnr):ni(pnr))-1);
          end
        end;
     elseif plane == 1 
        y = coord;
        if readz == 0 
           x = h(pnr)*((ib(pnr):ni(pnr))-1);
           z = zmin(pnr) + h(pnr)*((jb(pnr):nj(pnr))-1);
        else
          n1 = size(im,1);
          n2 = size(im,2);
% make x two-dimensional (same size as z and im)
          x = zeros(n1,n2);
          for i=1:n1
            x(i,:) = h(pnr)*((ib(pnr):ni(pnr))-1);
          end
        end;
     elseif plane == 2
        x = h(pnr)*((ib(pnr):ni(pnr))-1);
        y = h(pnr)*((jb(pnr):nj(pnr))-1);
        z = coord;
     end;
  else
     disp(['Error: number of patches on file ' num2str(npatches) ' is smaller than input pnr = ' num2str(pnr)]);
  end;
else
   disp(['Error: could not open file ' imfile ]);
end;
