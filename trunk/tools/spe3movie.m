% SPE3MOVIE
%
% read all image files following a pattern, make countour plots and ...
%
function spe3movie(az,el,img)

pattern='spe3.cycle=%04i.x=5625.mag.sw4img';
patt2='spe3.cycle=%04i.z=0.mag.sw4img';
dt = 13/6904;
h=15;

mmax=0;

% levels
cmax=2.e-4;
lev=linspace(0,cmax,51);
lim=[0 cmax];

% read topo
topo=readimage('spe3.cycle=0000.z=0.topo.sw4img',1);

xmin=60*h;
xmax=9000-60*h;
ymin=60*h;
ymax=7000-60*h;

k0=97;
%for k=img
for k=3970:10:6000

% make a file name for the image data
  fn1 = sprintf(pattern, k);
  fn2 = sprintf(patt2, k);
  disp(['Filename #' num2str(k), ': ', fn2])
% read horizontal data
  [m3 x3 y3 z3]=readimage(fn2,1);
%  plot
  mp = max(max(m3));
  if (mp>mmax) mmax=mp; end;

  figure(1);
  clf;

% solution on topographic surface
  surf(y3,x3,topo',m3'); shading flat; 
  view(az,el);
  camzoom(1.7);
  caxis(lim);
  axis([0 7000 0 9000 900 2200] ); axis equal;
  axis off;
  hold on;
% sg boundary
  [n1 n2]=size(topo);
  plot3(ymin*ones(1,length(x3)),x3,topo(60,:),'k-','Linewidth',2);
  plot3(ymax*ones(1,length(x3)),x3,topo(n1-60,:),'k-','Linewidth',2);
  plot3(y3,xmin*ones(1,length(y3)),topo(:,60),'k-','Linewidth',2);
  plot3(y3,xmax*ones(1,length(y3)),topo(:,n2-60),'k-','Linewidth',2);
% interior grid lines size(topo)=[468, 601]
  for i=60:60:408
    plot3((i-1)*h*ones(1,length(x3)),x3,topo(i,:),'k-','Linewidth',1);
  end;
  for j=60:60:541
    plot3(y3,(j-1)*h*ones(1,length(y3)),topo(:,j),'k-','Linewidth',1);
  end;

  ts = sprintf('t = %5.2f; max = %e',dt*k,mp);
%  title(ts);
%  set(gca,'Fontsize',16);

% wait for input
%  input("Press return to continue...");

% make output file name (put images in subdirectory)
  ifn = sprintf('frames2/img%03i.tiff', k0);

% save image as png file
  print(ifn,'-dtiff')

  disp(['saved image file: ' ifn])
  
%  t0=t(1);
%  dt=t(2)-t(1);
%  writeusgs( pn1, bn1, ux, uy, uz, dt, t0, 0, 0);
  k0 = k0+1;
end;
disp(['Max mag:' num2str(mmax)]);
end
