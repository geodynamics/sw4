%
% READMAT
%
%         [topo,x,y]=readmat( imfile, verbose )
%
%                Input: imfile  - Image file name
%                       verbose - Set to 1 to display file header information.
%                Output: im      - The image, as a 2D array.
%                        x, y, z - The spatial coordinates of the image, one of these is a scalar.
%                        plane   - 0: x=const, 1: y=const, 2: z=const
%                        t       - Simulation time at which the image was output.
%                        timestring - String holding creation date
function [pstr,topo,x,y]=readmat( imfile, verbose )
if nargin < 2
   verbose= 0;
end;

fd=fopen(imfile,'r');
if fd ~= -1 
% Read header
   magic   =fread(fd,1,'int');
   prec    =fread(fd,1,'int');
   alpha   =fread(fd,1,'double');
   len     =fread(fd,1,'int');
   pstr    =fread(fd,[1 len],'uchar');
   pstr = num2str(pstr,'%c');
   nblocks  =fread(fd,1,'int');
% Display header
   if verbose == 1
      disp(['Header: magic  = ' num2str(magic)  ' prec= ' num2str(prec) ' alpha= ' num2str(alpha)]);   
      disp(['        len = ' num2str(len)]);
% not sure how to parse the string
      disp([' pstr=' num2str(pstr(1:len))]);
      disp([' nblocks= ' num2str(nblocks) ]);
   end;
   for p=1:nblocks
      h(p) = 800.0; % hard coded for now
      Ni(p) =fread(fd,1,'int');
      Nj(p) =fread(fd,1,'int');
      Nk(p) =fread(fd,1,'int');
      if verbose == 1
         disp(['    block nr ' num2str(p) ' has h = ' num2str(h(p)) ' Ni = ' num2str(Ni(p)) ' Nj = ' num2str(Nj(p)) ' Nk = ' num2str(Nk(p))]);
      end;
   end;
% Read data
else
   disp(['Error: could not open file ' imfile ]);
end;
fclose(fd);
