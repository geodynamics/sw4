% Rectangular grid, x and y are vectors.
function gridplotrect(x,y,col,setwidth)
if nargin < 4
  setwidth = 0;
end;
if nargin < 3
  col = 'k';
end;
m=length(x);
n=length(y);
ox=ones(1,m);
oy=ones(1,n);
h=plot(x(1)*oy,y,col);
if setwidth==1
  set(h,'LineWidth',2.0);
end;
hold on
for i=2:m
  h=plot(x(i)*oy,y,col);
  if setwidth==1
    set(h,'LineWidth',2.0);
  end;
end;
for j=1:n
  h=plot(x,y(j)*ox,col);
  if setwidth==1
     set(h,'LineWidth',2.0);
  end;
end;
hold off

