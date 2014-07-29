%
% PLOTTS
%
%              plotts(basename, bk, ek, comp, voff)
%
%       Input: 
%              basename - Name of receiever data file
%              bk: int
%              ek: int
%              comp: component number in usgs file
%              voff: vertical offset
%               
function plotts(basename, bk, ek, comp, voff)

% assume comp=5 for now (pressure)
  comp=5;

% first find global min and max in time series
pmax = -9e7;
pmin = 9e7;

for k=bk:ek
% make string for filename
  fname=sprintf('%s%i.txt', basename, k);

  [t r u1 u2 u3 p]=readusgs(fname);
  lmax = max(p);
  lmin = min(p);
  if (lmax > pmax) 
    pmax = lmax;
  end;
  if (lmin < pmin) 
    pmin = lmin;
  end;

end;
pmax;
pmin;
grange=pmax-pmin;

% now read the time series again and plot them with consistent scales
clf;
hold on;
for k=bk:ek
% make string for filename
  fname=sprintf('%s%i.txt', basename, k);

  [t r u1 u2 u3 p]=readusgs(fname);
  lmax = max(p);
  lmin = min(p);
  lrange = lmax-lmin;

  if (lmax > pmax) 
    pmax = lmax;
  end;
  if (lmin < pmin) 
    pmin = lmin;
  end;

  scale = lrange/grange;
  scale = 1.;
  scale = grange/lrange;

  rscale = 10000;
  k0 = 5;
% reduced time: k-k0 is range on a 10 km scale
  toff = abs(k-k0)*10/0.4;

  plot(t - toff,rscale*((k-k0)*voff + p.*scale),'k','linewidth',2);

end;
hold off;

% tmp
return

if (erase ~= 0)
  clf;
end
% east component
subplot(3,1,1)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,ux,colorstring);
%set(h,'LineWidth',2.0)
set(gca,'FontSize',20)
axis tight;

% north component
subplot(3,1,2)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,uy,colorstring);
%set(h,'LineWidth',2.0)
set(gca,'FontSize',20)
axis tight;

% up component
subplot(3,1,3)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,uz,colorstring);
%set(h,'LineWidth',2.0)
set(gca,'FontSize',20)
axis tight;
