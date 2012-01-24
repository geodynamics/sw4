%
%  FKWPPCOMP
%     Compute the difference between two 3-component time series. The first time series is stored in 3 sac files
%     (as output by the fk program). The second time series is given in one text file using the usgs format (as 
%     is output by the wpp program). The wpp output files is assumed to hold (x,y,z) components with z directed 
%     downwards. The fk time series are assumed to contain the radial, transverse and vertical (upwards) components.
%
%  APPROACH
%     The wpp time series is first rotated into radial, transverse and upwards components. The fk time series is then
%     interpolated onto the same time levels as the wpp time series. Finally, the vector L2 and max norms of the
%     difference is output together with the vector L2 and max norms of the fk time series.
%
%  USAGE:
%     [e2norm, emaxnorm, u2norm, umaxnorm]=fkwppcomp( fkbase, wppfile, plotit, loh3, sigma, strike )
%
%  ARGUMENTS:
%     Input:
%          fkbase:  base name for sac files. The actual files must be named fkbase.[rtz]. NOT used when loh3=1
%          wppfile: file name of wpp output file
%          plotit:  Plot the three components of the fk solution as well as the error (fk-wpp)
%          loh3:    Optional argument:
%                     0: (default), read output from fk. 1: read output from loh3exact(0.1)
%          sigma:   Optional argument: spread in Gaussian time function sigma=1/freq, freq is WPP frequency parameter
%                     default value: 0.1
%          strike:  Optional argument:
%                   strike angle [degrees] for the reciever location. Default: 53.1301
%     Output:
%          e2norm:    Vector L2-norm of difference
%          emaxnorm:  Vector max-norm of difference
%          u2norm:    Vector L2-norm of fk time series
%          umaxnorm:  Vector max-norm of fk time series
%
function [e2norm, emaxnorm, u2norm, umaxnorm]=fkwppcomp( fkbase, wppfile, plotit, loh3, sigma, strike )

if nargin < 6 % standard location of the reciever for the LOH1-3 test cases
  strike = 53.1301;
end

if nargin < 5
  sigma=0.1;
end

if nargin < 4
  loh3 = 0;
end

if loh3==0
% read fk files
  [tfk, radfk, tranfk, vertfk] = rtzfilter( fkbase, 0.01, 0 );
else
  [tfk, radfk, tranfk, vertfk] = loh3exact( sigma ); % sigma=0.1: change accordingly
end

% test
%plot(tfk,radfk,tfk,tranfk,tfk,vertfk)

% read wpp file
[tw, uxw, uyw, uzw]=readusgs( wppfile );

% strike angle defines radial and tangential components
ca = cos(strike*pi/180);
sa = sin(strike*pi/180);

% rotate wpp data
urw = ca*uxw + sa*uyw;
utw = -sa*uxw + ca*uyw;
uvw = uzw; % positive downwards

% test
%plot(tw, urw, tw, utw, tw, uvw)

% interpolate fk solution onto wpp time levels
[radi, trani, upi] = fkinterp(tfk, radfk, tranfk, vertfk, tw);

% test
%plot(tw, radi, tw, trani, tw, upi)

% compute errors

nt = length(tw);
%e2norm = sqrt((sumsq(radi-urw) + sumsq(trani-utw) + sumsq(upi-uvw))/nt);
%u2norm = sqrt((sumsq(radi) + sumsq(trani) + sumsq(upi))/nt);
e2norm = sqrt((sum((radi-urw).^2) + sum((trani-utw).^2) + sum((upi-uvw).^2))/nt);
u2norm = sqrt((sum(radi.^2) + sum(trani.^2) + sum(upi.^2))/nt);

ermax = max(abs(radi-urw));
etmax = max(abs(trani-utw));
evmax = max(abs(upi-uvw));

emaxnorm = max([ermax, etmax, evmax]);

urmax = max(abs(radi));
utmax = max(abs(trani));
uvmax = max(abs(upi));

umaxnorm = max([urmax, utmax, uvmax]);

if plotit == 1
  subplot(3,1,1);
%  plot(tw,radi,'k',tw,urw,'r',tw,radi-urw,'b');
%  legend('Radial fk','WPP','Difference');
  plot(tw,radi,'k',tw,urw,'r');
  legend('Radial fk','WPP');

  subplot(3,1,2);
%  plot(tw,trani,'k',tw,utw,'r',tw,trani-utw,'b');
%  legend('Transverse fk','WPP','Difference');
  plot(tw,trani,'k',tw,utw,'r');
  legend('Transverse fk','WPP');

  subplot(3,1,3);
%  plot(tw,upi,'k',tw,uvw,'r',tw,upi-uvw,'b');
%  legend('Vertical fk','WPP','Difference');
  plot(tw,upi,'k',tw,uvw,'r');
  legend('Vertical fk','WPP');

end
