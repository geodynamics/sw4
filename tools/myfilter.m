function [mf] = myfilter(b,a,u)
N=length(u);
x1=0;
x2=0;
y1=0;
y2=0;
mf=0*(1:N);
for i=1:N
  op = b(1)*u(i) + b(2)*x1 + b(3)*x2 - (a(2)*y1 + a(3)*y2);
  y2=y1;
  y1=op;
  x2=x1;
  x1=u(i);
  mf(i)=op;
end
