%
% DIRECTFORM2
%
% evaluate one second order section of a 3 stage IIR filter
%
%  mf = directform2(b,a,u)
%
% Input: b(3): vector with nominator coefficients
%        a(3): vector with denominator coefficients, a(1)=1
%        u(n): sequence to be filtered
% Output: mf(n): filtered sequence

function [mf] = directform2(b,a,u)
N=length(u);
x1=0;
x2=0;
y1=0;
y2=0;
mf = 0*(1:N);
wn1 = 0;
wn2 = 0;
for i=1:N
  wn = u(i) - a(2)*wn1 - a(3)*wn2;
  mf(i) = b(1)*wn + b(2)*wn1 + b(3)*wn2;
  wn2 = wn1;
  wn1 = wn;
end
