%
% PLOTUSGS2
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios and plot it in 3 subwindows
%
%              plotusgs2(filename, evtlon, evtlat, colorstring, erasefirst, tmin, tmax)
%
%       Input: filename - Name of receiever data file
%              evtlon: Source (epicenter) longitude
%              evtlat: Source (epicenter) latitude
%              colorstring: string passed to plot, like 'r' for red lines
%              erasefirst: 0 does a 'hold on' for the current plot, otherwise erases the current figure
%              tmin: starting time (for plot)
%              tmax: ending time (for plot)
%               
function plotusgs2( filename, evtlon, evtlat, colorstring, erase, tmin, tmax )

tshift = 0;

[t ux uy uz]=readusgs(filename, evtlon, evtlat, 1);

if (erase ~= 0)
  clf;
end
% east component
subplot(3,1,1)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,ux,colorstring);
set(h,'LineWidth',2.0)
set(gca,'FontSize',16)
title("Radial(top), transverse(middle), vertical(bottom)");
umax = max(ux);
umin = min(ux);
axis([tmin tmax umin umax]);
%axis tight;

% north component
subplot(3,1,2)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,uy,colorstring);
set(h,'LineWidth',2.0)
set(gca,'FontSize',16)
umax = max(uy);
umin = min(uy);
axis([tmin tmax umin umax]);
%axis tight;

% up component
subplot(3,1,3)
if (erase == 0)
  hold on;
end
h=plot(t+tshift,uz,colorstring);
set(h,'LineWidth',2.0)
set(gca,'FontSize',16)
umax = max(uz);
umin = min(uz);
axis([tmin tmax umin umax]);
%axis tight;
