%
% PLOTUSGS
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios and plot it in 3 subwindows
%
%              plotusgs(filename, colorstring, erasefirst, timeshift)
%
%       Input: filename - Name of receiever data file
%              colorstring: string passed to plot, like 'r' for red lines
%              erasefirst: 0 does a 'hold on' for the current plot, otherwise erases the current figure
%              timeshift:  change independent variable to be t+timeshift
%               
function plotusgs( filename, colorstring, erase, tshift )
   lw=2.0;
if nargin < 4
   tshift = 0;
end;

if nargin < 3
  erase = 1;
end;

if nargin < 2
  colorstring='b';
end;

[t ux uy uz]=readusgs(filename);

## if (erase ~= 0)
##   clf;
## end
% east component
#subplot(3,1,1)
figure(1)
if (erase == 0)
  hold on;
else
  clf;
end
h=plot(t+tshift,ux,colorstring);
if lw >0
  set(h,'LineWidth',lw)
end;
set(gca,'FontSize',20)
axis tight;

% north component
#subplot(3,1,2)
figure(2)
if (erase == 0)
  hold on;
else
  clf;
end
h=plot(t+tshift,uy,colorstring);
if lw > 0
  set(h,'LineWidth',lw)
end;
set(gca,'FontSize',20)
axis tight;

% up component
#subplot(3,1,3)
figure(3)
if (erase == 0)
  hold on;
else
  clf;
end
h=plot(t+tshift,uz,colorstring);
if lw > 0
  set(h,'LineWidth',lw)
end;
set(gca,'FontSize',20)
axis tight;
