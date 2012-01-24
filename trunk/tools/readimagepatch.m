%
% READIMAGEPATCH
%
%  Reads one grid patch on an image file
%
%      [im,x,y]=readimagepatch( fil, inr )
%
%        Input: fil  - Name of image file.
%               inr  - Patch number on file, 1..maxpatch. Default value is 1.
%        Output: im  - The image patch
%                x,y - Corresponding x- and y-coordinates.
function [im,x,y]=readimagepatch( fil, inr )

if nargin < 2
  inr = 1;
end
fd=fopen(fil,'r');
pr=fread(fd,1,'int');
ni=fread(fd,1,'int');
if inr > ni
   disp( 'Error image nr too large');
else
   for i=1:ni
      h(i)  = fread(fd,1,'double');
      ib(i) = fread(fd,1,'int');
      ie(i) = fread(fd,1,'int');
      jb(i) = fread(fd,1,'int');
      je(i) = fread(fd,1,'int');
   end;
   for i=1:inr-1
      fseek(fd,(ie(i)-ib(i)+1)*(je(i)-jb(i)+1)*pr,0);
   end;
   if pr == 4 
      im0 = fread(fd,[ie(inr)-ib(inr)+1 je(inr)-jb(inr)+1],'float');
   else
      im0 = fread(fd,[ie(inr)-ib(inr)+1 je(inr)-jb(inr)+1],'double');
   end;
   x  = ((ib(inr):ie(inr))-1)*h(inr);
   y  = ((jb(inr):je(inr))-1)*h(inr);
   fclose(fd);
% transpose im0 and return result in im
   im = im0';
end;
