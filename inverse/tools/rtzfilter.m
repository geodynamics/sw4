%
%  RTZFILTER:
%         Read and optionally plot a three component seismogram in SAC format
%  USAGE:
%         [t,radf,tranf,upf]=rtzfilter( basename, scale, fc, style, doplot, hld )
%
%  ARGUMENTS: style, doplot and hld are optional arguments
%        Input:
%              basename: SAC file names are assumed to be basename.r, basename.t, basename.z
%              scale:    Apply scale factor to all three components
%              fc:       Corner frequency for low-pass Butterworth filter. No filtering if fc<=0
%              style:    Plot style, e.g., 'b-' for blue solid lines
%              doplot:   1 for plotting 
%              hld:      1 for holding the current plot (for overlaying several plots)
%        Output:
%              t:        Vector of time-levels
%              radf:     Vector of radial component
%              tranf:    Vector of transverse component
%              upf:      Vector of vertical component 
%
   function [t,radf,tranf,upf]=rtzfilter( basename, scale, fc, style, doplot, hld )
if nargin < 6
  hld = 0;
end
if nargin < 5
  doplot = 0;
end
if nargin < 4
  style = 'k-';
end;

% radial-component
fname=sprintf('%s.r',basename);
[rad dt lat lon t0]=readsac(fname);

% transverse-component
fname=sprintf('%s.t',basename);
[tran dt lat lon t0]=readsac(fname);

% vertical-component
fname=sprintf('%s.z',basename);
[up dt lat lon t0]=readsac(fname);

% time
len = length(rad)-1;
t = t0+dt*(0:len);

% scale input
rad = scale.*rad;
tran = scale.*tran;
up = -scale.*up; %positive down

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

