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
   function plotusgs( fname, colorstring, erase )

fd=fopen(fname,'r');
for i=1:12,fgetl(fd);end;
q=fscanf(fd,'%f');
t=q(1:4:end,1);
ux=q(2:4:end,1);
uy=q(3:4:end,1);
uz=q(4:4:end,1);
if (erase ~= 0)
  clf;
end
% east component
subplot(3,1,1)
if (erase == 0)
  hold on;
end
plot(t,ux,colorstring);
axis tight;

% north component
subplot(3,1,2)
if (erase == 0)
  hold on;
end
plot(t,uy,colorstring);
axis tight;

% up component
subplot(3,1,3)
if (erase == 0)
  hold on;
end
plot(t,uz,colorstring);
axis tight;
