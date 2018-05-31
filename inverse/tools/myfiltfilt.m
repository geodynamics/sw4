% 
% ZERO PHASE FORWARD AND BACKWARD FILTERING
%
% mf = myfiltfilt(b,a,u)
%
% Input: b(1:3): numerator filter coefficients
%        a(1:3): denominator filter coefficients 
%        u(1:N): time series to be filtered (N=length(u))
%
% Output: mf(1:N): filtered time series
%
function mf = myfiltfilt(b,a,u)
N=length(u);
x1=u(1);
x2=u(1);
y1=u(1);
y2=u(1);
%uf=0*(1:N);
mf=0*(1:N);
for i=1:N
  op = b(1)*u(i) + b(2)*x1 + b(3)*x2 - (a(2)*y1 + a(3)*y2);
  y2=y1;
  y1=op;
  x2=x1;
  x1=u(i);
%  uf(i)=op;
  mf(i)=op;
end

x1=mf(N);
x2=mf(N);
y1=mf(N);
y2=mf(N);

%x1=uf(N);
%x2=uf(N);
%y1=uf(N);
%y2=uf(N);

for i=N:-1:1
%  op = b(1)*uf(i) + b(2)*x1 + b(3)*x2 - (a(2)*y1 + a(3)*y2);
  op = b(1)*mf(i) + b(2)*x1 + b(3)*x2 - (a(2)*y1 + a(3)*y2);
  y2=y1;
  y1=op;
  x2=x1;
%  x1=uf(i);
  x1=mf(i);
  mf(i)=op;
end
