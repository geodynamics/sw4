%
% READALLREFINEMENTS
%
%  Reads one grid patch on an image file
%
%      [im,x,y]=readAllRefinements( fil, nlev )
%
%        Input: fil  - Name of image file.
%        Output: im  - The last image patch on the file
%                x,y - Corresponding x- and y-coordinates.
function [im,x,y]=readAllRefinements( fil, nlev, pmin, pmax )

% optional arguments
userbounds=1;
if nargin < 3
  userbounds=0;
end
if nargin < 2
  nlev=20;
end

fd=fopen(fil,'r');
pr=fread(fd,1,'int');
ni=fread(fd,1,'int');

clear figure;
clear im;
clear zmin;

for i=1:ni
     h(i)  = fread(fd,1,'double');
     ib(i) = fread(fd,1,'int');
     ie(i) = fread(fd,1,'int');
     jb(i) = fread(fd,1,'int');
     je(i) = fread(fd,1,'int');
end;

zmin=zeros(ni);
zmin(ni) = 0.0;

for i=(ni-1):-1:1
    zmin(i) = ((je(i+1)-jb(i+1)))*h(i+1)+zmin(i+1);
end

% first read to find ranges
umin = 1e10;
umax = -1e10;
for i=1:ni
   if pr == 4 
      im0 = fread(fd,[ie(i)-ib(i)+1 je(i)-jb(i)+1],'float');
   else
      im0 = fread(fd,[ie(i)-ib(i)+1 je(i)-jb(i)+1],'double');
   end;
   umax1 = max(max(im0));
   umin1 = min(min(im0));
   umax = max(umax,umax1);
   umin = min(umin,umin1);
end
fclose(fd);

disp(sprintf('Global umin=%e, umax=%e', umin, umax))
if userbounds>0
  umin = pmin;
  umax = pmax;
  disp(sprintf('Using specified bounds: umin=%e, umax=%e', umin, umax))
end

% now plot it
fd=fopen(fil,'r');
pr=fread(fd,1,'int');
ni=fread(fd,1,'int');

clear figure;
clear im;
clear zmin;

for i=1:ni
     h(i)  = fread(fd,1,'double');
     ib(i) = fread(fd,1,'int');
     ie(i) = fread(fd,1,'int');
     jb(i) = fread(fd,1,'int');
     je(i) = fread(fd,1,'int');
end;

zmin=zeros(ni);
zmin(ni) = 0.0;

for i=(ni-1):-1:1
    zmin(i) = ((je(i+1)-jb(i+1)))*h(i+1)+zmin(i+1);
end

xmin=1e10;
xmax=-1e10;
ymin=1e10;
ymax=-1e10;


for i=1:ni
   if pr == 4 
      im0 = fread(fd,[ie(i)-ib(i)+1 je(i)-jb(i)+1],'float');
   else
      im0 = fread(fd,[ie(i)-ib(i)+1 je(i)-jb(i)+1],'double');
   end;
   x  = ((ib(i):ie(i))-1)*h(i);
   y  = ((jb(i):je(i))-1)*h(i)+zmin(i);
   
   xmin = min(xmin,x(ib(i)));
   ymin = min(ymin,y(jb(i)));
   xmax = max(xmax,x(ie(i)));
   ymax = max(ymax,y(je(i)));

% transpose im0 and return result in im
   im = im0';

   contour(x,-y,im, linspace(umin,umax,nlev));
   caxis([umin umax]);
   hold on;  

end

axis([xmin, xmax, -ymax, -ymin]);
%colorbar;
axis equal;
hold off;

fclose(fd);
