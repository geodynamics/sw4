%
% ZEROPHASE3
%
% apply a filter described by three second order sections, first forwards and then backwards
%
   %  mf = zerophase3(b0,a0,b1,a1,b2,a2,u)
%
% Input: 
%        b0(3): vector with nominator coefficients (SOS #0)
%        a0(3): vector with denominator coefficients, a0(1)=1 (SOS #0)
%        b1(3): vector with nominator coefficients (SOS #1)
%        a1(3): vector with denominator coefficients, a1(1)=1 (SOS #1)
%        b2(3): vector with nominator coefficients (SOS #2)
%        a2(3): vector with denominator coefficients, a2(1)=1 (SOS #2)
%        u(n): sequence to be filtered
% Output: mf(n): filtered sequence

    function [mf] = zerophase3(b0,a0,b1,a1,b2,a2,u)
N=length(u);
mf = 0*(1:N);

# first SOS forwards
wn1 = 0;
wn2 = 0;
for i=1:N
  wn = u(i) - a0(2)*wn1 - a0(3)*wn2;
  mf(i) = b0(1)*wn + b0(2)*wn1 + b0(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# second SOS forwards
wn1 = 0;
wn2 = 0;
for i=1:N
  wn = mf(i) - a1(2)*wn1 - a1(3)*wn2;
  mf(i) = b1(1)*wn + b1(2)*wn1 + b1(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# third SOS forwards (over-writes mf)
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
  wn = mf(i) - a0(2)*wn1 - a0(3)*wn2;
  mf(i) = b0(1)*wn + b0(2)*wn1 + b0(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# second SOS backward (over-writes mf)
wn1 = 0;
wn2 = 0;
for i=N:-1:1
  wn = mf(i) - a1(2)*wn1 - a1(3)*wn2;
  mf(i) = b1(1)*wn + b1(2)*wn1 + b1(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end

# third SOS backward (over-writes mf)
wn1 = 0;
wn2 = 0;
for i=N:-1:1
  wn = mf(i) - a2(2)*wn1 - a2(3)*wn2;
  mf(i) = b2(1)*wn + b2(2)*wn1 + b2(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end
