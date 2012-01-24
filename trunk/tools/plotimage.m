%
% PLOTIMAGE
%
%     plotimage( fil, cvals )
%
%   Plots the image on file 'fil' with contourf, using the contour levels cvals.
%   The boundary of each grid patch is outlined in black. A vector cvals can be 
%   obtained from function imageinfo.
%
function plotimage( fil, cvals )
fd=fopen(fil,'r');
pr=fread(fd,1,'int');
nb=fread(fd,1,'int');
fclose(fd);

for b=1:nb
   [im,x,y] = readimagepatch(fil,b);
   contour(x,y,im,cvals);
   if b==1
      hold on;
   end;
   [n m]=size(im);
   plot(x, y(1)*ones(size(x)),'k')
   plot(x, y(n)*ones(size(x)),'k')
   plot(x(1)*ones(size(y)), y,'k')
   plot(x(m)*ones(size(y)), y,'k')
end;
hold off;
