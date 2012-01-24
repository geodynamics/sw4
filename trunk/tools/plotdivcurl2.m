function [im,x,y]=plotdivcurl2( basename )
clf
gy2=readimagepatch('imageY.cycle=000.x=4975.grid',2);
gz2=readimagepatch('imageZ.cycle=000.x=4975.grid',2);

div=readimagepatch(sprintf('%s.div', basename) ,2);
curl=readimagepatch(sprintf('%s.curl', basename) ,2);

ymin = min(min(gy2));
ymax = max(max(gy2));
dy=0.02*(ymax-ymin);
ymin = ymin - dy;
ymax= ymax + dy;

zmin = min(min(gz2));
zmax = max(max(gz2));
dz=0.02*(zmax-zmin);
zmin = zmin - dz;
zmax= zmax + dz;

cmin=0;
cmax=6e-4;
dmax=1.85e-4;
nlev=10;


subplot(1,2,1);
contourf(gy2,-gz2,div,linspace(-dmax,dmax,2*nlev+1));
axis([ymin, ymax, -zmax,-zmin]);
caxis([-dmax,dmax]);
colorbar;
set(gca,'FontSize',16);
xlabel('Y');ylabel('-Z')
divmin=min(min(div));
divmax=max(max(div));
title(sprintf('%10.2e < div < %10.2e', divmin, divmax));
hold on;
plot(gy2(1,:),-gz2(1,:),'k');
hold off;

subplot(1,2,2);
contourf(gy2,-gz2,curl,linspace(cmin,cmax,nlev+1));
axis([ymin, ymax, -zmax,-zmin]);
caxis([cmin,cmax]);
colorbar;
set(gca,'FontSize',16);
xlabel('Y');ylabel('-Z')
curlmax=max(max(curl));
title(sprintf('|curl| < %10.2e', curlmax));
hold on;
plot(gy2(1,:),-gz2(1,:),'k');
hold off;
