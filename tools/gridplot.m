function gridplot(x,y,col,setwidth)
if nargin < 4
  setwidth = 0;
end;
if nargin < 3
  col = 'k';
end;
[m n]=size(x);
h=plot(x(1,:),y(1,:),col);
if setwidth==1
   set(h,'LineWidth',2.0);
end;
hold on
for i=2:m
   h=plot(x(i,:),y(i,:),col);
   if setwidth==1
     set(h,'LineWidth',2.0);
   end;
end;
for j=1:n
   h=plot(x(:,j),y(:,j),col);
   if setwidth==1
     set(h,'LineWidth',2.0);
   end;
end;
hold off
