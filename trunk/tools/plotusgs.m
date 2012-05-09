%
% PLOTUSGS
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios and plot it in 3 subwindows
%
%              plotusgs(filename, colorstring, erasefirst)
%
%       Input: filename - Name of receiever data file
%              colorstring: string passed to plot, like 'r' for red lines
%              erasefirst: 0 does a 'hold on' for the current plot, otherwise erases the current figure
%               
function plotusgs( filename, colorstring, erase )

[t ux uy uz]=readusgs(filename);

if (erase ~= 0)
  clf;
end
% east component
subplot(3,1,1)
if (erase == 0)
  hold on;
end
h=plot(t,ux,colorstring);
set(h,'LineWidth',1.2)
axis tight;

% north component
subplot(3,1,2)
if (erase == 0)
  hold on;
end
h=plot(t,uy,colorstring);
set(h,'LineWidth',1.2)
axis tight;

% up component
subplot(3,1,3)
if (erase == 0)
  hold on;
end
h=plot(t,uz,colorstring);
set(h,'LineWidth',1.2)
axis tight;
