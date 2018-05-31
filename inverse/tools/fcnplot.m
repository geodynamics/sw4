%
%  FCNPLOT
%
%   Get and plot source time dependence function g(t).
%
%          [t,g]=fcnplot( fcn, freq, tmin, tmax, t0, col, doplot, dt )
%
%        Input:  fcn    - Name of source function (string), same name as
%                         used in SW4. Also some additional functions
%                         not in SW4, see code below.
%                freq   - Frequency parameter, scaling as in SW4.
%                tmin, tmax - Time interval 
%                t0     - Time shift, i.e., get and plot g(t-t0)
%                col    - Plotting color.
%                doplot -  0 -> No plot is produced.
%
%        Output: t - Times where g is sampled.
%                g - Source time function at the points given in t.
%
   function [t,g]=fcnplot( fcn, freq, tmin, tmax, t0, col, doplot, dt )
if nargin < 8
  dt = (tmax-tmin)/4999;
end
if nargin < 7
  doplot = 1;
end;
if nargin < 6
   col = 'k';
end;
if nargin < 5
   t0 = 0;
end;
omega = freq;
n = round((tmax-tmin)/dt);
t = tmin + (0:n)*dt;

if strcmp(fcn,'GaussianInt') == 1
%   g = omega/sqrt(2*pi)*erf(t)+0.5;
   g = 0.5*(1+erf(omega.*(t-t0)./sqrt(2)));
   yl= 0;
   yh= 1;
elseif strcmp(fcn,'RickerInt') == 1
   g = (t-t0).*exp(-(pi*omega.*(t-t0)).^2);
   yl= -0.2;
   yh=  0.2;
elseif strcmp( fcn,'Gaussian') == 1
   g  = omega/sqrt(2*pi)*exp(-((t-t0).*omega).^2/2);
   yl = 0.0;
   yh = omega/sqrt(2*pi);
elseif strcmp( fcn, 'Ricker' ) == 1
   g  = (2*(pi*omega.*(t-t0)).^2-1).*exp(-(pi*omega.*(t-t0)).^2);
   yl =-1.1;
   yh = 0.5;
elseif strcmp( fcn, 'Triangle' ) == 1
   for i=1:n+1
     ta = omega*(t(i)-t0);
     if ta>0 & ta < 1
        g(i) = 16/(pi*pi)*(sin(pi*ta)-(1/9)*sin(3*pi*ta)+(1/25)*sin(5*pi*ta)-(1/49)*sin(7*pi*ta));
     else
        g(i) = 0;
     end;
   end;
   yl = 0.0;
   yh = 2.0;
elseif strcmp( fcn, 'Sawtooth' ) == 1
   for i=1:n+1
     ta = omega*(t(i)-t0);
     if ta>0 & ta < 1
        g(i) =  8/(pi*pi)*(sin(2*pi*ta)-(1/9)*sin(6*pi*ta)+(1/25)*sin(10*pi*ta)-(1/49)*sin(14*pi*ta));
     else
        g(i) = 0;
     end;
   end;
   yl=-1;
   yh=1;
elseif strcmp( fcn, 'Ramp' ) == 1
   for i=1:n+1
     ta = omega*(t(i)-t0);
     if ta<0 
        g(i) = 0;
     elseif ta <= 1
        g(i) = 0.5*(1-cos(pi*ta));
     else
        g(i) = 1;
     end;
   end;
   yl=-0.1;
   yh=1.1;
elseif strcmp( fcn, 'DRamp' ) == 1
   for i=1:n+1
     ta = omega*(t(i)-t0);
     if ta<0 
        g(i) = 0;
     elseif ta <= 1
        g(i) = omega*0.5*(pi*sin(pi*ta));
     else
        g(i) = 0;
     end;
   end;
elseif strcmp( fcn, 'Smoothwave' ) == 1
   for i=1:n+1
     ta = omega*(t(i)-t0);
     if ta>0 & ta < 1
        g(i) = (2187/8)*ta.^3-(10935/8)*ta.^4+(19683/8)*ta.^5-(15309/8)*ta.^6+(2187/4)*ta.^7;
     else
        g(i) = 0;
     end;
   end;
   yl=-1.2;
   yh=1.2;
elseif strcmp( fcn, 'Brune' ) == 1
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta < 0
         g(i) = 0;
      else
         g(i) = 1-exp(-ta).*(1+ta);
      end;
   end;
   yl=-0.1;
   yh= 1.1;
elseif strcmp( fcn, 'DBrune' ) == 1
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta < 0
         g(i) = 0;
      else
         g(i) = omega*ta.*exp(-ta);
      end;
   end;
   yl=-0.1;
   yh= omega/2;
elseif strcmp( fcn, 'BruneSmoothed' ) == 1
   x0 = 2.31;
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta < 0
         g(i) = 0;
      elseif ta < x0
         g(i) = 1-exp(-ta).*(1+ta+0.5*ta.^2-1.5*ta.^3./x0+1.5*ta.^4./(x0*x0)-(1/(2*x0^3))*ta.^5);
      else
         g(i) = 1-exp(-ta).*(1+ta);
      end;
   end;
   yl=-0.1;
   yh= 1.1;
elseif strcmp( fcn, 'DBruneSmoothed' ) == 1
   x0 = 2.31;
%   x0 = 5;
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta < 0
         g(i) = 0;
      elseif ta <= x0
         dp = 1 + ta - 3*1.5*ta.^2./x0+4*1.5*ta.^3./(x0*x0)-(5/(2*x0^3))*ta.^4;
         g(i) = omega*exp(-ta).*((1+ta+0.5*ta.^2-1.5*ta.^3./x0+1.5*ta.^4./(x0*x0)-(1/(2*x0^3))*ta.^5)-dp);
      else
         g(i) = omega*ta.*exp(-ta);
      end;
   end;
   yl=-0.1;
   yh= omega/2;
elseif strcmp( fcn, 'VerySmoothBump' ) == 1
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta<0
         g(i) = 0;
      elseif ta>1 
         g(i) = 0;
      else
        g(i) = - 1024*ta.^10 + 5120*ta.^9 - 10240*ta.^8 + 10240*ta.^7 - 5120*ta.^6 + 1024*ta.^5;
      end;
   end;
   yl=-0.1;
   yh =1.1;
elseif strcmp( fcn, 'C6SmoothBump' ) == 1
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta<0
         g(i) = 0;
      elseif ta>1 
         g(i) = 0;
      else
        g(i) = 51480*ta.^7*(1-ta).^7;
      end;
   end;
   yl=-0.1;
   yh =3.5;
elseif strcmp( fcn, 'GaussianWindow' ) == 1
   ncyc = 5;
   g = sin(omega*t).*exp(-0.5*(omega*(t-t0)./ncyc).^2);
   yl=-1;
   yh= 1;
elseif strcmp( fcn, 'Liu' ) == 1
%  tt(i) = tmin + (i-1)/(n-1)*(tmax-tmin);
   tau = 2*pi/omega;
   tau1 = 0.13*tau;
   tau2 = tau-tau1;
   yl=-0.1;
   yh= 1.1;
   for i=1:n+1
     ta = (t(i)-t0);
     if ta < 0 
       g(i)= 0;
     elseif ta >= tau
       g(i)= 1;
     else
       ipi = 1.0/pi;
       cn = 1.0/(1.4*tau1+1.2*tau1*ipi + 0.3*tau2);
       if ta <= tau1 
	 g(i)= cn*(0.7*ta-0.7*tau1*ipi*sin(pi*ta/tau1)-1.2*tau1*ipi*(cos(0.5*pi*ta/tau1)-1));
       elseif ta <= 2*tau1 
	 g(i)= cn*(1.0*ta-0.3*tau1+1.2*tau1*ipi - 0.7*tau1*ipi*sin(pi*ta/tau1)+0.3*tau2*ipi*sin(pi*(ta-tau1)/tau2));
       elseif ta <= tau
	 g(i)= cn*(0.3*ta+1.1*tau1+1.2*tau1*ipi+0.3*tau2*ipi*sin(pi*(ta-tau1)/tau2));
       end;
     end;
   end;
elseif strcmp( fcn,'UnscaledGaussian') == 1
   g  = exp(-((t-t0).*omega).^2/2);
   yl = 0.0;
   yh = 1.0;
elseif strcmp( fcn,'TwoSine') == 1
   for i=1:n+1
      ta = omega*(t(i)-t0);
      if ta<0
         g(i) = 0;
      elseif ta>1 
         g(i) = 0;
      else
         g(i) = sin(2*pi*ta) - 0.5*sin(4*pi*ta);
      end;
   end;
   yl=-1.5;
   yh =1.5;
end;
if doplot == 1
   [h]=plot(t,g,col);
   set(h,'LineWidth',1.0);
   %get(h)
   set(gca,'FontSize',25)
   xlabel('t')
   ylabel('g(t)')
   title( [fcn ' \omega=' num2str(omega) ' t_0=' num2str(t0)]);
   axis([tmin tmax yl yh]);
end;

