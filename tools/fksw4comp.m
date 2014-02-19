%
%  FKSW4COMP
%     Compute the difference between two 3-component time series. The first time series is stored in 3 sac files
%     (as output by the fk program). The second time series is given in one text file using the usgs format (as 
%     is output by the sw4 program). The sw4 output files is assumed to hold (x,y,z) components with z directed 
%     downwards. The fk time series are assumed to contain the radial, transverse and vertical (upwards) components.
%
%  APPROACH
%     The sw4 time series is first rotated into radial, transverse and upwards components. The fk time series is then
%     interpolated onto the same time levels as the sw4 time series. Finally, the vector L2 and max norms of the
%     difference is output together with the vector L2 and max norms of the fk time series.
%
%  USAGE:
%     [e2norm, emaxnorm, u2norm, umaxnorm]=fksw4comp( fkbase, sw4file, plotit, tshift, loh, sigma, sac, strike )
%
%  ARGUMENTS:
%     Input:
%          fkbase:  base name for sac files. The actual files must be named fkbase.[rtz]. NOT used when loh=1
%          sw4file: file name of sw4 output file
%          plotit:  Plot the three components of the fk solution as well as the error (fk-sw4)
%          tshift:  Optional argument:
%                    shift fk time series by this amount (default value tshift=0)
%          loh:    Optional argument:
%                     0: (default), read output from fk. 
%                     1: read output from loh1exact
%                     3: read output from loh3exact
%          sigma:   Optional argument: only used for LOH3: spread in Gaussian time function sigma=1/freq, 
%                     freq is SW4 frequency parameter; default value: sigma=0.06
%          sac:     Optional argument: (default value sac=0)
%                     0: read sw4 results from 1 usgs formatted file with name 'sw4file'
%                     1: read sw4 results from 3 sac files with base name 'sw4file'
%          strike:  Optional argument:
%                   strike angle [degrees] for the reciever location. Default: 53.1301
%     Output:
%          e2norm:    Vector L2-norm of difference
%          emaxnorm:  Vector max-norm of difference
%          u2norm:    Vector L2-norm of fk time series
%          umaxnorm:  Vector max-norm of fk time series
%
function [e2norm, emaxnorm, u2norm, umaxnorm]=fksw4comp( fkbase, sw4file, plotit, tshift, loh, sigma, sac, strike )

if nargin < 8 % standard location of the reciever for the LOH1-3 test cases
  strike = 53.1301;
end

if nargin < 7
  sac=0;
end

if nargin < 6
  sigma=0.06;
end

if nargin < 5
  loh = 0;
end

if nargin < 4
  tshift = 0;
end

if (loh==0)
% read fk files
  [tfk, radfk, tranfk, vertfk] = rtzfilter( fkbase, 0.01, 0 );
  tfk = tfk + tshift;
elseif (loh == 3)
  [tfk, radfk, tranfk, vertfk] = loh3exact( sigma ); % sigma=0.1: change accordingly
elseif (loh == 1)
  [tfk, radfk, tranfk, vertfk] = loh1exact( sigma ); % sigma=0.1: change accordingly
end

% test
%plot(tfk,radfk,tfk,tranfk,tfk,vertfk)

% read sw4 file
if (sac==1)
% x file name
  fname=sprintf('%s.x', sw4file);
  [uxw dtw dum dum t0w] = readsac(fname);
  uxw = uxw';
% x file name
  fname=sprintf('%s.y', sw4file);
  [uyw dtw dum dum t0w] = readsac(fname);
  uyw = uyw';
% x file name
  fname=sprintf('%s.z', sw4file);
  [uzw dtw dum dum t0w] = readsac(fname);
  uzw = uzw';
% construct time levels
  ntw = length(uxw);
  tw = t0w + dtw*[0:ntw-1];
else
  [tw, uxw, uyw, uzw]=readusgs( sw4file );
end

% strike angle defines radial and tangential components
ca = cos(strike*pi/180);
sa = sin(strike*pi/180);

% rotate sw4 data
urw = ca*uxw + sa*uyw;
utw = -sa*uxw + ca*uyw;
uvw = uzw; % positive downwards

% test
%plot(tw, urw, tw, utw, tw, uvw)

% interpolate fk solution onto sw4 time levels
[radi, trani, upi] = fkinterp(tfk, radfk, tranfk, vertfk, tw);

tfkmax = max(tfk);
twmax  = max(tw);

% ignore any part of the sw4 solution that extends beyonf the fk solution
if twmax>tfkmax
  twmin = tw(1);
  dtw = (twmax-twmin)/(length(tw)-1);
  ntrunc = floor(1+(tfkmax-twmin)/dtw);
  
  nt = ntrunc;
else
  nt = length(tw);
end

% test
%plot(tw, radi, tw, trani, tw, upi)

% compute errors

%e2norm = sqrt((sumsq(radi-urw) + sumsq(trani-utw) + sumsq(upi-uvw))/nt);
%u2norm = sqrt((sumsq(radi) + sumsq(trani) + sumsq(upi))/nt);
e2norm = sqrt((sum((radi(1:nt)-urw(1:nt)).^2) + sum((trani(1:nt)-utw(1:nt)).^2) + sum((upi(1:nt)-uvw(1:nt)).^2))/nt);
u2norm = sqrt((sum(radi(1:nt).^2) + sum(trani(1:nt).^2) + sum(upi(1:nt).^2))/nt);

ermax = max(abs(radi(1:nt)-urw(1:nt)));
etmax = max(abs(trani(1:nt)-utw(1:nt)));
evmax = max(abs(upi(1:nt)-uvw(1:nt)));

emaxnorm = max([ermax, etmax, evmax]);

urmax = max(abs(radi(1:nt)));
utmax = max(abs(trani(1:nt)));
uvmax = max(abs(upi(1:nt)));

umaxnorm = max([urmax, utmax, uvmax]);

if plotit == 1
  subplot(3,1,1);
%  plot(tw,radi,'k',tw,urw,'r',tw,radi-urw,'b');
h=plot(tw(1:nt),radi(1:nt),'k',tw(1:nt),urw(1:nt),'r',tw(1:nt),radi(1:nt)-urw(1:nt),'b');
  legend('Radial fk','FD','Difference');
  set(h,"Linewidth",2)
  set(gca,"Fontsize",16);

  subplot(3,1,2);
%  plot(tw,trani,'k',tw,utw,'r',tw,trani-utw,'b');
h=plot(tw(1:nt),trani(1:nt),'k',tw(1:nt),utw(1:nt),'r',tw(1:nt),trani(1:nt)-utw(1:nt),'b');
  legend('Transverse fk','FD','Difference');
  set(h,"Linewidth",2)
  set(gca,"Fontsize",16);

  subplot(3,1,3);
%  plot(tw,upi,'k',tw,uvw,'r',tw,upi-uvw,'b');
h=plot(tw(1:nt),upi(1:nt),'k',tw(1:nt),uvw(1:nt),'r', tw(1:nt), upi(1:nt)-uvw(1:nt),'b');
  legend('Vertical fk','FD','Difference');
  set(h,"Linewidth",2)
  set(gca,"Fontsize",16);

end
