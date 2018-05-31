%
%  PLOTFILTER
%
%   Get and plot source time dependence function g(t).
%
%         [t,ef,nf,uf]=plotfilter( fname, fc, style, doplot, hld )
%
%        Input:
%        Output:
   function [t,ef,nf,uf]=plotfilter( fname, fc, style, doplot, hld )
if nargin < 5
  hld = 0;
end
if nargin < 4
  doplot = 1;
end
if nargin < 3
  style = 'k-';
end;

[t e n u]=readusgs(fname);
dt = t(2)-t(1);

if fc > 0
  [b a]=mybutter2(2*dt*fc);
  ef = myfiltfilt(b,a,e);
  nf = myfiltfilt(b,a,n);
  uf = myfiltfilt(b,a,u);
else
  ef = e;
  nf = n;
  uf = u;
end

if doplot == 1
   subplot(3,1,1);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,ef,style);

   rng = max(ef)-min(ef);
   axis([t(1) t(end) min(ef)-0.2*rng max(ef)+0.2*rng]);

   subplot(3,1,2);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,nf,style);

   rng = max(nf)-min(nf);
   axis([t(1) t(end) min(nf)-0.2*rng max(nf)+0.2*rng]);

   subplot(3,1,3);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(t,uf,style);

   rng = max(uf)-min(uf);
   axis([t(1) t(end) min(uf)-0.2*rng max(uf)+0.2*rng]);
end;

