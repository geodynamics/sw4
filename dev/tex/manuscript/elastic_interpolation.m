x = [1 2 3 4 5 6 7];
y1 = [1 1 1 1 1 1 1];
y2 = 2*y1;
y3 = 3*y1;
y4 = 4*y1;
y5 = 5*y1;
y6 = 6*y1;
y7 = 7*y1;
y = x;
x1 = y1;
x2 = y2;
x3 = y3;
x4 = y4;
x5 = y5;
x6 = y6;
x7 = y7;

figure(4)
plot(x,y1,'b',x,y2,'r--',x,y3,'b',x,y4,'r--',x,y5,'b',x,y6,'r--',x,y7,'b','linewidth',2)
hold on
plot(x1,y,'b',x2,y,'r--',x3,y,'b',x4,y,'r--',x5,y,'b',x6,y,'r--',x7,y,'b','linewidth',2)

f1 = plot([1],[1],'bo','MarkerFaceColor','b','Markersize',8);
plot([1],[3],'bo','MarkerFaceColor','b','Markersize',8)
plot([1],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([1],[7],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[1],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[3],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[7],'bo','MarkerFaceColor','b','Markersize',8)
plot([5],[1],'bo','MarkerFaceColor','b','Markersize',8)
plot([5],[3],'bo','MarkerFaceColor','b','Markersize',8)
plot([5],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([5],[7],'bo','MarkerFaceColor','b','Markersize',8)
plot([7],[1],'bo','MarkerFaceColor','b','Markersize',8)
plot([7],[3],'bo','MarkerFaceColor','b','Markersize',8)
plot([7],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([7],[7],'bo','MarkerFaceColor','b','Markersize',8)
axis([0 8 0 8])
f2 = plot([4],[4],'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2);
legend([f1 f2],'coarse grids', 'fine grid')
set(gca,'fontsize',16)

axis off

figure(3)
plot(x,y1,'b',x,y2,'r--',x,y3,'b',x,y4,'r--',x,y5,'b',x,y6,'r--',x,y7,'b','linewidth',2)
hold on
plot(x1,y,'b',x2,y,'r--',x3,y,'b',x4,y,'r--',x5,y,'b',x6,y,'r--',x7,y,'b','linewidth',2)

f1 = plot([3],[1],'bo','MarkerFaceColor','b','Markersize',8);
plot([3],[3],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([3],[7],'bo','MarkerFaceColor','b','Markersize',8)
axis([0 8 0 8])
f2 = plot([3],[4],'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2);
legend([f1 f2],'coarse grids', 'fine grid')
set(gca,'fontsize',16)

axis off

figure(2)
plot(x,y1,'b',x,y2,'r--',x,y3,'b',x,y4,'r--',x,y5,'b',x,y6,'r--',x,y7,'b','linewidth',2)
hold on
plot(x1,y,'b',x2,y,'r--',x3,y,'b',x4,y,'r--',x5,y,'b',x6,y,'r--',x7,y,'b','linewidth',2)

f1 = plot([1],[5],'bo','MarkerFaceColor','b','Markersize',8);
plot([3],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([5],[5],'bo','MarkerFaceColor','b','Markersize',8)
plot([7],[5],'bo','MarkerFaceColor','b','Markersize',8)
axis([0 8 0 8])
f2 = plot([4],[5],'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2);
legend([f1 f2],'coarse grids', 'fine grid')
set(gca,'fontsize',16)

axis off

figure(1)
plot(x,y1,'b',x,y2,'r--',x,y3,'b',x,y4,'r--',x,y5,'b',x,y6,'r--',x,y7,'b','linewidth',2)
hold on
plot(x1,y,'b',x2,y,'r--',x3,y,'b',x4,y,'r--',x5,y,'b',x6,y,'r--',x7,y,'b','linewidth',2)

axis([0 8 0 8])
f1 = plot([3],[5],'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2);
f2 = plot([3],[5],'bo','MarkerFaceColor','b','Markersize',8);
legend([f1 f2],'coarse grids', 'fine grid')
set(gca,'fontsize',16)

axis off

figure(5)
plot(x,y1,'r--',x,y2,'b',x,y3,'r--',x,y4,'b',x,y5,'r--',x,y6,'b',x,y7,'r--','linewidth',2)
hold on
plot(x1,y,'r--',x2,y,'b',x3,y,'r--',x4,y,'b',x5,y,'r--',x6,y,'b',x7,y,'r--','linewidth',2)

axis([0 8 0 8])

plot(x,y1,'rs',...
    x,y2,'rs',...
    x,y3,'rs',...
    x,y4,'rs',...
    x,y5,'rs',...
    x,y6,'rs',...
    x,y7,'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2)
plot(x1,y,'rs',...
    x2,y,'rs',...
    x3,y,'rs',...
    x4,y,'rs',...
    x5,y,'rs',...
    x6,y,'rs',...
    x7,y,'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2)

f1 = plot(x1(1),y(1),'rs','MarkerEdgeColor','r','Markersize',12,'linewidth',2);
f2 = plot([4],[4],'bo','MarkerFaceColor','b','Markersize',8,'linewidth',2);
legend([f1 f2], 'fine grids','coarse grid')
set(gca,'fontsize',16)

axis off

% Mass matrix
figure(6)
y = 30*ones(1,4);
x = [1 2 3 4];
plot(x(1),y(1),'mo','MarkerFaceColor','m','Markersize',10)
hold on
plot(x(2:4),y(2:4),'b^','MarkerFaceColor','b','Markersize',10)
y = 29*ones(1,5);
x = [1 2 3 4 5];
plot(x(1),y(1),'b^',x(3:5),y(3:5),'b^','MarkerFaceColor','b','Markersize',10)
plot(x(2),y(2),'mo','MarkerFaceColor','m','Markersize',10)
axis([1 30 1 30])
y = 28*ones(1,6);
x = [1 2 3 4 5 6];
plot(x(1:2),y(1:2),'b^',x(4:6),y(4:6),'b^','MarkerFaceColor','b','Markersize',10)
plot(x(3),y(3),'mo','MarkerFaceColor','m','Markersize',10)
for i = 27:-1:4
    y = i*ones(1,7);
    x = [28-i 29-i 30-i 31-i 32-i 33-i 34-i];
    plot(x(1:3),y(1:3),'b^',x(5:7),y(5:7),'b^','MarkerFaceColor','b','Markersize',10)
    plot(x(4),y(4),'mo','MarkerFaceColor','m','Markersize',10)
end
y = 3*ones(1,6);
x = [25 26 27 28 29 30];
plot(x(1:3),y(1:3),'b^',x(5:6),y(5:6),'b^','MarkerFaceColor','b','Markersize',10)
plot(x(4),y(4),'mo','MarkerFaceColor','m','Markersize',10)
y = 2*ones(1,5);
x = [26 27 28 29 30];
plot(x(1:3),y(1:3),'b^',x(5),y(5),'b^','MarkerFaceColor','b','Markersize',10)
plot(x(4),y(4),'mo','MarkerFaceColor','m','Markersize',10)
y = ones(1,4);
x = [27 28 29 30];
plot(x(1:3),y(1:3),'b^','MarkerFaceColor','b','Markersize',10)
plot(x(4),y(4),'mo','MarkerFaceColor','m','Markersize',10)
set(gca,'fontsize',16)

% each 'm*'
figure(7)
x = 1:21;
y1 = ones(1,21);
y2 = 2*ones(1,21);
y3 = 3*ones(1,21);
plot(x(1:9),y1(1:9),'ks',x(1:9),y2(1:9),'ks',x(1:9),y3(1:9),'ks','MarkerFaceColor','k','Markersize',10)
hold on
plot(x(10:12),y1(10:12),'ro',x(10:12),y2(10:12),'ro',x(10:12),y3(10:12),'ro','MarkerFaceColor','r','Markersize',10)
plot(x(13:21),y1(13:21),'ks',x(13:21),y2(13:21),'ks',x(13:21),y3(13:21),'ks','MarkerFaceColor','k','Markersize',10)
axis([0 22 -15 15])
axis off