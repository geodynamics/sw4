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
t=fread(fd,1,'double');
plane=fread(fd,1,'int');
coord   =fread(fd,1,'double');
mode    =fread(fd,1,'int');
gridinfo=fread(fd,1,'int');
fclose(fd);

for b=1:nb
   [im,x,y,z] = readimage(fil,b);
   if plane==0
     contour(y,z,im,cvals);
   elseif plane==1
     contour(x,z,im,cvals);
   elseif plane==2
     contour(x,y,im,cvals);
   end
   if b==1
      hold on;
   end;
   [n m]=size(im);
% make a frame
  if plane==0
    if (b==nb) && (gridinfo == 1)
      plot(y(1,:), z(1,:),'k')
      plot(y(n,:), z(n,:),'k')
      plot(y(:,1), z(:,1),'k')
      plot(y(:,m), z(:,m),'k')
    else
      plot(y, z(1)*ones(size(y)),'k')
      plot(y, z(n)*ones(size(y)),'k')
      plot(y(1)*ones(size(z)), z,'k')
      plot(y(m)*ones(size(z)), z,'k')
    end;
  elseif plane==1
    if (b==nb) && (gridinfo == 1)
      plot(x(1,:), z(1,:),'k')
      plot(x(n,:), z(n,:),'k')
      plot(x(:,1), z(:,1),'k')
      plot(x(:,m), z(:,m),'k')
    else
      plot(x, z(1)*ones(size(x)),'k')
      plot(x, z(n)*ones(size(x)),'k')
      plot(x(1)*ones(size(z)), z,'k')
      plot(x(m)*ones(size(z)), z,'k')
    end;
  elseif plane==2
    plot(x, y(1)*ones(size(x)),'k')
    plot(x, y(n)*ones(size(x)),'k')
    plot(x(1)*ones(size(y)), y,'k')
    plot(x(m)*ones(size(y)), y,'k')
  end;
end;
hold off;
