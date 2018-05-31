%
% ZEROPHASE2
%
% apply two second order sections, first forwards and then backwards
%
%  mf = zerophase2(b1,a1,b2,a2,u)
%
% Input: 
%        b1(3): vector with nominator coefficients (SOS #1)
%        a1(3): vector with denominator coefficients, a(1)=1 (SOS #1)
%        b2(3): vector with nominator coefficients (SOS #2)
%        a2(3): vector with denominator coefficients, a(1)=1 (SOS #2)
%        u(n): sequence to be filtered
% Output: mf(n): filtered sequence

function [mf] = zerophase2(b1,a1,b2,a2,u)
N=length(u);
mf = 0*(1:N);

# first SOS forwards
wn1 = 0;
wn2 = 0;
for i=1:N
  wn = u(i) - a1(2)*wn1 - a1(3)*wn2;
  mf(i) = b1(1)*wn + b1(2)*wn1 + b1(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# Second SOS forwards (over-writes mf)
wn1 = 0;
wn2 = 0;
for i=1:N
  wn = mf(i) - a2(2)*wn1 - a2(3)*wn2;
  mf(i) = b2(1)*wn + b2(2)*wn1 + b2(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# First SOS backward (over-writes mf)
wn1 = 0;
wn2 = 0;
for i=N:-1:1
  wn = mf(i) - a1(2)*wn1 - a1(3)*wn2;
  mf(i) = b1(1)*wn + b1(2)*wn1 + b1(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# Second SOS backward (over-writes mf)
wn1 = 0;
wn2 = 0;
for i=N:-1:1
  wn = mf(i) - a2(2)*wn1 - a2(3)*wn2;
  mf(i) = b2(1)*wn + b2(2)*wn1 + b2(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end
