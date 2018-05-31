%
% NORMUSGS
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios.Then compute the norm of the time series for t>t1
%
%              [rel_nrm, maxt_norm] = normusgs(filename, t1, doplot, linlog)
%
%       Input: filename - Name of receiever data file
%              t1: [optional], default t1=0. Otherwise, start time for evaluating the norm
%              doplot: [optional] default=0. If doplot=1, plot the norm as fcn of time
%              linlog: [optional] default=0 for linear scale. If linlog=1, use a 
%                      logarithmic scale
%               
function [rel_nrm, maxt_norm] = normusgs( filename, t1, doplot, linlog )

if nargin < 4
   linlog = 0;
end;

if nargin < 3
   doplot = 0;
end;

if nargin < 2
   t1 = 0;
end;

[t ux uy uz]=readusgs(filename);

[dum it]=min(abs(t-t1));

i1 = length(t);

nrm = sqrt(ux.^2 + uy.^2 + uz.^2);
maxt_nrm = max(nrm);

if (doplot == 1)
  if (linlog == 0)
    h=plot(t(it:i1), nrm(it:i1)/maxt_nrm, "b");
  else
    h=semilogy(t(it:i1), nrm(it:i1)/maxt_nrm, "b");
  endif
%set(h,'LineWidth',2.0)
  set(gca,'FontSize',20)
  axis tight;
endif

mnrm = max(nrm(it:i1));
rel_nrm = mnrm/maxt_nrm;
%printf("max norm of entire time series: %e\n", maxt_nrm);
%printf("For t>%e, relative max norm of time series: %e\n", t1, rel_nrm);
