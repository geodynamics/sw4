% 
% BUTTERWORTH 2ND ORDER LOWPASS FILTER COEFFICIENTS
%
% [b,a] = mybutter2(Wn)
%
% Input: Wn = 2*dt*freq, where dt is the sampling interval of the time series (to be filtered),
%                        and freq is the corner frequency [Hz]
% Output: b(1:3): numerator coefficients
%         a(1:3): denominator coefficients with a(1)=1
%
function [b,a] = mybutter2(Wn)
fr = 2/Wn;
omegac = tan(pi/fr);
c = 1 + sqrt(2)*omegac + omegac^2;
b=0*(1:3);
a=0*(1:3);
b(1) = omegac^2 / c;
b(2) = 2*b(1);
b(3) = b(1);
a(1) = 1;
a(2) = 2*(omegac^2 - 1)/c;
a(3) = (1 - sqrt(2)*omegac + omegac^2)/c;
  
