% 
% BUTTERWORTH BANDPASS DIGITAL FILTER COEFFICIENTS, given as two second order sections (SOS)
% that can be used to filter a time series using the direct form II algorithm in directform2.m
%
% Note: this routine uses a proptotype filter consisting of two complex 
% conjugated poles: q=exp(i alpha), qc=exp(-i alpha). It does not handle the
% case of a single real pole. The single real pole is handled by mybutterodd.m
%
% [b1,a1,b2,a2] = mybuttereven(f1,f2,dt,alpha)
%
% Input: 
%        f1: low corner frequency [Hz]
%        f2: high corner frequency [Hz]
%        dt: time step [s] of the time series (to be filtered),
%        alpha: angle of the pole [rad] (pi/2 < alpha < pi).
% 
% Output: 
%         b1(1:3): numerator coefficients, SOS 1
%         a1(1:3): denominator coefficients with a1(1)=1, SOS 1
%         b2(1:3): numerator coefficients, SOS 2
%         a2(1:3): denominator coefficients with a2(1)=1, SOS 2
%
function [b1,a1,b2,a2] = mybuttereven(f1,f2,dt,alpha)
if (alpha >= pi || alpha <= pi/2)
  printf("pole angle alpha=%e out of range\n");
  return
end
% imaginary unit
iu = 1i;

%pre-warp the corner frequencies
om1 = tan(pi*dt*f1);
om2 = tan(pi*dt*f2);

b = om2 - om1;
p = om1*om2;

% pole #1
q = cos(alpha) + iu*sin(alpha);

% analog bp filter coeff transfer fcn are saved as N(s)/D(s), 
% N(s) = n(1) + n(2)*s + n(3)*s^2
% D(s) = d(1) + d(2)*s + d(3)*s^2

% initialize storage
n1=0*[1:3];
d1=0*[1:3];
n2=0*[1:3];
d2=0*[1:3];

% roots of the two quadratics
s1 = (q*b + sqrt(q^2*b^2 - 4*p))/2;
s2 = (q*b - sqrt(q^2*b^2 - 4*p))/2;
s3 = (conj(q)*b + sqrt(conj(q)^2*b^2 - 4*p))/2;
s4 = (conj(q)*b - sqrt(conj(q)^2*b^2 - 4*p))/2;

% if Re(s1) >= 0 or Re(s2)>= the filter is unstable
if (real(s1) >= 0 || real(s2) >= 0)
  printf("WARNING: the analog filter has poles in the positive half-plane. s1=%e%+ei, s2=%e%+ei\n", real(s1), imag(s1), real(s2), imag(s2));
end
 
% check the algebra
%printf("P1: q*b=%e%+ei, s1+s2=%e%+ei, s1*s2=%e%+ei, p=%e\n", real(q*b), imag(q*b), real(s1+s2), imag(s1+s2), real(s1*s2), imag(s1*s2), p);
%printf("P2: conj(q)*b=%e%+ei, s3+s4=%e%+ei, s3*s4=%e%+ei, p=%e\n", real(conj(q)*b), imag(conj(q)*b), real(s3+s4), imag(s3+s4), real(s3*s4), imag(s3*s4), p);
%printf("Q1: 2*Re(s1) = %e, |s1|^2 = %e\n", 2*real(s1), abs(s1)^2)
%printf("Q2: 2*Re(s2) = %e, |s2|^2 = %e\n", 2*real(s2), abs(s2)^2)
%printf("Q1Q2, s^3: %e, s^2: %e, s^1: %e, s^0: %e\n", -2*(real(s1)+real(s2)), abs(s1)^2+abs(s2)^2+4*real(s1)*real(s2), -2*(real(s1)*abs(s2)^2 + real(s2)*abs(s1)^2), abs(s1)^2*abs(s2)^2);
%printf("P1P2, s^3: %e, s^2: %e, s^1: %e, s^0: %e\n", -b*(q+conj(q)), 2*p + b^2*q*conj(q), -b*p*(q+conj(q)), p^2);

% These are the coefficients in SOS #1 (Numerator and denominator)
n1(2) = b;
d1(1) = abs(s1)^2;
d1(2) = -2*real(s1);
d1(3) = 1;

% These are the coefficients in SOS #2 (Numerator and denominator)
n2(2) = b;
d2(1) = abs(s2)^2;
d2(2) = -2*real(s2);
d2(3) = 1;

% allocate space for output arrays
b1=0*[1:3];
a1=0*[1:3];
b2=0*[1:3];
a2=0*[1:3];

% transform analog to digital by the transformation AD(s) = (1-s)/(1+s)
% first SOS
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

% second SOS
% normalization factor
c = d2(1) + d2(2) + d2(3); 
% denominator
a2(1) = 1;
a2(2) = 2*(d2(1)-d2(3))/c;
a2(3) = (d2(1) - d2(2) + d2(3))/c;
% nominator
b2(1) = (n2(1) + n2(2) + n2(3))/c;
b2(2) = 2*(n2(1)-n2(3))/c;
b2(3) = (n2(1) - n2(2) + n2(3))/c;
  
