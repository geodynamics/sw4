%
%  WPPFILTER
%
%        Input:
%        Output:
   function [t,radf,tranf,upf]=wppfilter( filename, er, fc, style, doplot, hld )
if nargin < 6
  hld = 0;
end
if nargin < 5
  doplot = 0;
end
if nargin < 4
  style = 'k-';
end;

[t xc yc zc]=readusgs(filename);
% time step is assumed constant
dt = t(2)-t(1);

% rotate input to radial, transverse and vertical components

% radial-component
  rad=er(1)*xc+er(2)*yc;

% transverse-component
tran=-er(2)*xc+er(1)*yc;

% vertical-component
up = -zc;

% optionally filter
if fc > 0
  [b a]=mybutter2(2*dt*fc);
  radf = myfiltfilt(b,a,rad);
  tranf = myfiltfilt(b,a,tran);
  upf = myfiltfilt(b,a,up);
else
  radf=rad;
  tranf=tran;
  upf=up;
end

if doplot == 1
   subplot(3,1,1);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,radf,style);
   axis([t(1) t(end) min(radf) max(radf)]);

   subplot(3,1,2);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,tranf,style);
   axis([t(1) t(end) min(tranf) max(tranf)]);

   subplot(3,1,3);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,upf,style);
   axis([t(1) t(end) min(upf) max(upf)]);
end;

