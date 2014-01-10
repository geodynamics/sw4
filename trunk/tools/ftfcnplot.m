%
%  FTFCNPLOT
%
%    Plot Fourier transform of source time function.
%
%           [om,amp] = ftfcnplot( fcn, freq, tmin, tmax, t0, col, lin )
%
%        Input:  fcn    - Name of source function (string), same name as
%                         used in SW4. Also some additional functions
%                         not in SW4, see function 'fcnplot'.
%                freq   - Frequency parameter, scaling as in SW4.
%                tmin, tmax - Time interval 
%                t0     - Time shift.
%                col    - Plotting color.
%                lin    - 1: linear y-axis, 0: logarithmic y-axis
%
%        Output: om  - Frequencies [Hz]
%                amp - Amplitude of Fourier coefficients at om.
%
function [om fr]=ftfcnplot( fcn, freq, tmin, tmax, t0, col, lin )

[t,g] = fcnplot( fcn, freq, tmin, tmax, t0, 'k', 0 );
u = fft(g);
N = length(u);
dt = t(2)-t(1);
L = N*dt;

r = N/2;
f = u(1:r)./r;
% constant component is a factor of 2 larger
%f(1) = f(1)/2.;
% Scale with half the interval length 
f = f*0.5*L;
% We want frequency, not angular frequency
freqaxis = ((0:r-1))./L;

if lin==1
  plot(freqaxis,abs(f),col);
else
  eps=1e-99;
  semilogy( freqaxis,abs(f)+eps, col );
end
set(gca,'FontSize',15);
xlabel('frequency [Hz]')
ylabel('Magnitude of Fourier transform')
title(sprintf('Fourier transform of %s with freq=%g', fcn, freq))
fr = abs(f);
om = freqaxis;
