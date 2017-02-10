%
% PLOTIMAGE
%
%     plotimage( fil, machineformat, cvals )
%
%   Plots the image on file 'fil' with contour, using the contour levels cvals.
%   If cvals is ommitted, 21 countour levels will be obtained through the imageinfo fcn.
%
%   The boundary of each grid patch is outlined in black. A vector cvals can be 
%   obtained from function imageinfo.
%
%   Input:
%         fil:                       Name of image file
%         machineformat (optional):  Passed to fopen to read big endian, little endian, etc
%         cvals (optional):          Vector of countour levels to plot
%
   function plotimage( fil, machineformat, cvals )
if nargin < 2
   machineformat='native';
end;
if nargin < 3
   nc=21;
   cvals = imageinfo(fil,nc,0,machineformat);
end;

fd=fopen(fil,'r',machineformat);
pr=fread(fd,1,'int');
nb=fread(fd,1,'int');
t=fread(fd,1,'double');
plane=fread(fd,1,'int');
coord   =fread(fd,1,'double');
mode    =fread(fd,1,'int');
gridinfo=fread(fd,1,'int');
fclose(fd);

x1min=1e9;
x2min=1e9;
x1max=-1e9;
x2max=-1e9;

for b=1:nb
	[im,x,y,z] = readimage(fil,b,0,machineformat);
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
      x1mi=min(min(y));
      x1ma=max(max(y));
      x2mi=min(min(z));
      x2ma=max(max(z));
    else
      plot(y, z(1)*ones(size(y)),'k')
      plot(y, z(n)*ones(size(y)),'k')
      plot(y(1)*ones(size(z)), z,'k')
      plot(y(m)*ones(size(z)), z,'k')
      x1mi=min(y);
      x1ma=max(y);
      x2mi=min(z);
      x2ma=max(z);
    end;
  elseif plane==1
    if (b==nb) && (gridinfo == 1)
      plot(x(1,:), z(1,:),'k')
      plot(x(n,:), z(n,:),'k')
      plot(x(:,1), z(:,1),'k')
      plot(x(:,m), z(:,m),'k')
      x1mi=min(min(x));
      x1ma=max(max(x));
      x2mi=min(min(z));
      x2ma=max(max(z));
    else
      plot(x, z(1)*ones(size(x)),'k')
      plot(x, z(n)*ones(size(x)),'k')
      plot(x(1)*ones(size(z)), z,'k')
      plot(x(m)*ones(size(z)), z,'k')
      x1mi=min(x);
      x1ma=max(x);
      x2mi=min(z);
      x2ma=max(z);
    end;
  elseif plane==2
    plot(x, y(1)*ones(size(x)),'k')
    plot(x, y(n)*ones(size(x)),'k')
    plot(x(1)*ones(size(y)), y,'k')
    plot(x(m)*ones(size(y)), y,'k')
    x1mi=min(x);
    x1ma=max(x);
    x2mi=min(y);
    x2ma=max(y);
  end;
% update global min/max
  x1min = min(x1min,x1mi);
  x1max = max(x1max,x1ma);
  x2min = min(x2min,x2mi);
  x2max = max(x2max,x2ma);
end;
axis([x1min x1max x2min x2max]);
hold off;
axis ij; % flip z-axis to point downwards
