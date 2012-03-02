%
%  IMAGEINFO
%
%  Get some info about an image.
%
%       Syntax:
%               [cvals,np,mx,mn]=imageinfo( fil, nc )
%
%       Input:  fil - Name of image file
%               nc  - Number of countour levels (default=10).
%       Output: cvals - Contour level vector to be used with function plotimage
%               np    - Number of patches on file.
%               mx    - Maximum data value over all patches.
%               mn    - Minimum data value over all patches.
%
function [cvals,np,mx,mn]=imageinfo( fil, nc )

if nargin < 2
   nc = 10;
end;

fd = fopen(fil,'r');
pr = fread(fd,1,'int');
np = fread(fd,1,'int');
fclose(fd);

for b=1:np
   [im,x,y] = readimagepatch(fil,b);
   if b == 1 
      mx = max(max(im));
      mn = min(min(im));
   else
      mx = max([mx max(max(im))]);
      mn = min([mn min(max(im))]);
   end;
end;
cvals = linspace(mn,mx,nc);
