% SPE3MOVIE
%
% read all image files following a pattern, make countour plots and ...
%
function spe3movie(img)

pattern='spe3.cycle=%04i.x=5625.mag.sw4img';
patt2='spe3.cycle=%04i.z=0.mag.sw4img';
dt = 13/6904;

mmax=0;

% levels
cmax=1.89e-4;
lev=linspace(0,cmax,51);
lim=[0 cmax];

k0=0;
%for k=3000:50:4000
for k=img

% make a file name for the image data
  fn1 = sprintf(pattern, k);
  fn2 = sprintf(patt2, k);
  printf("Filename #%i, '%s', '%s'\n", k, fn1, fn2)
  % read image (patch 1, 2)
  [m1 x1 y1 z1]=readimage(fn1,1);
  [m2 x2 y2 z2]=readimage(fn1,2);
% read horizontal data
  [m3 x3 y3 z3]=readimage(fn2,1);
%  plot
  maxp1 = max(max(m1));
  maxp2 = max(max(m2));
  mp = max(maxp1,maxp2);
  if (mp>mmax) mmax=mp; endif;

% vertical cross-section
  figure(1);
  clf;
  contour(y1,-z1,m1,lev); 
  hold on;
  contour(y2,-z2,m2,lev); 
  caxis(lim);
  axis([0 7000 -1700 2200]); axis equal
% sg boundary
%  plot([5000 5000],[-17000 0],'k-','Linewidth',1);
%  plot([25000 25000],[-17000 0],'k-','Linewidth',1);
%  plot([0 30000],[-12000 -12000],'k-','Linewidth',1);

  ts = sprintf('t = %5.2f; max = %e',dt*k,mp);
  title(ts);
  set(gca,'Fontsize',16);

% horizontal cross-section
  figure(2);
  clf;
  contour(x3,y3,m3,lev); 
  hold on;
  axis([0 9000 0 7000]); axis equal
% sg boundary
%  plot([5000 5000],[-17000 0],'k-','Linewidth',1);
%  plot([25000 25000],[-17000 0],'k-','Linewidth',1);
%  plot([0 30000],[-12000 -12000],'k-','Linewidth',1);

  ts = sprintf('t = %5.2f; max = %e',dt*k,mp);
  title(ts);
  set(gca,'Fontsize',16);

% wait for input
  input("Press return to continue...");

% make output file name (put images in subdirectory)
  ifn = sprintf('frames/img%03i.png', k0);

% save image as png file
%  print(ifn,'-dpng')

  printf("read file '%s', image file '%s'\n", fn1, ifn)
  
%  t0=t(1);
%  dt=t(2)-t(1);
%  writeusgs( pn1, bn1, ux, uy, uz, dt, t0, 0, 0);
  k0 = k0+1;
endfor;
printf('Max mag: %e\n', mmax);
end
