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
function plotusgs( filename, colorstring, erase, tshift, winL, winR )
lw=2.0;
if nargin < 6
   winL= 0;
   winR=-1;
end;
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
size(ux);
size(t);

n = length(t);
winfcn = ones(n,1);
if winL < winR
   a = pi/(winR-winL);
   b = -a*(winR+winL)/2;
   winfcn = zeros(n,1);
   for i=1:n
      if winL < t(i) && t(i) < winR 
	 winfcn(i) = cos(a*t(i)+b)^10;
      end;
   end;
end;
if (erase ~= 0)
  clf;
end
% east component
subplot(3,1,1)
%figure(1)
if (erase == 0)
  hold on;
end
%else
%  clf;
%end
h=plot(t+tshift,winfcn.*ux,colorstring);
if lw >0
   set(h,'LineWidth',lw);
end;
set(gca,'FontSize',20)
axis tight;

% north component
subplot(3,1,2)
%figure(2)
if (erase == 0)
  hold on;
end;
%else
%  clf;
%end
h=plot(t+tshift,winfcn.*uy,colorstring);
if lw >0
   set(h,'LineWidth',lw);
end;
set(gca,'FontSize',20)
axis tight;

% up component
subplot(3,1,3)
%figure(3)
if (erase == 0)
  hold on;
end;
%else
%  clf;
%end
h=plot(t+tshift,winfcn.*uz,colorstring);
if lw>0
   set(h,'LineWidth',lw)
end;
set(gca,'FontSize',20)
axis tight;

