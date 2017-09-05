%
% PLOTSAC
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios and plot it in 3 subwindows
%
%              plotsac(basename, colorstring, erasefirst, timeshift)
%
%       Input: basename - Name of receiever data file, basename.e, basename.n, basename.u
%              colorstring: string passed to plot, like 'r' for red lines
%              erasefirst: 0 does a 'hold on' for the current plot, otherwise erases the current figure
%              timeshift:  change independent variable to be t+timeshift
%               
function plotsac( basename, colorstring, erase, tshift )

if nargin < 4
   tshift = 0;
end;

if nargin < 3
  erase = 1;
end;

if nargin < 2
  colorstring='b';
end;

[ux dt, lat, lon]=readsac(sprintf("%s.%s", basename, "e"));
[uy]=readsac(sprintf("%s.%s", basename, "n"));
[uz]=readsac(sprintf("%s.%s", basename, "u"));

nt = length(ux);
t = dt*[0:nt-1];

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
