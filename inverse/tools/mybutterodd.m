% 
% BUTTERWORTH BANDPASS DIGITAL FILTER COEFFICIENTS, corresponding to the REAL POLE
% for an odd order filter
%
% Note:  This routine is used to compute the coefficients in a second order section (SOS), 
% which can be used to evaluate a digital filter using the directform2.m script.
%
% [b1,a1] = mybutterodd(f1,f2,dt)
%
% Input: 
%        f1: low corner frequency [Hz]
%        f2: high corner frequency [Hz]
%        dt: time step [s] of the time series (to be filtered),
% 
% Output: 
%         b1(1:3): numerator coefficients,
%         a1(1:3): denominator coefficients with a1(1)=1.
%
function [b1,a1] = mybutterodd(f1,f2,dt)

%pre-warp the corner frequencies
om1 = tan(pi*dt*f1);
om2 = tan(pi*dt*f2);

b = om2 - om1;
p = om1*om2;

% pole #1
q = -1;

% analog bp filter coeff transfer fcn are saved as N(s)/D(s), 
% N(s) = n(1) + n(2)*s + n(3)*s^2
% D(s) = d(1) + d(2)*s + d(3)*s^2

% initialize storage
n1=0*[1:3];
d1=0*[1:3];

% These are the coefficients in SOS #1 (Numerator and denominator)
n1(2) = b;
d1(1) = p;
d1(2) = b;
d1(3) = 1;

% allocate space for output arrays
b1=0*[1:3];
a1=0*[1:3];

% transform analog to digital by the transformation AD(s) = (1-s)/(1+s)
% normalization factor
c = d1(1) + d1(2) + d1(3); 
% denominator
a1(1) = 1;
a1(2) = 2*(d1(1)-d1(3))/c;
a1(3) = (d1(1) - d1(2) + d1(3))/c;
% nominator
b1(1) = (n1(1) + n1(2) + n1(3))/c;
b1(2) = 2*(n1(1)-n1(3))/c;
b1(3) = (n1(1) - n1(2) + n1(3))/c;

  
