clf
[d1 gy1 gz1]=readimagepatch('image.cycle=266.x=2500.div',1);
c1=readimagepatch('image.cycle=266.x=2500.curl',1);
subplot(1,2,1);
contourf(gy1,-gz1,d1);
axis equal
xlabel('X');ylabel('Y')
divmin=min(min(d1));
divmax=max(max(d1));
title(sprintf('%10.2e < div < %10.2e', divmin, divmax));

subplot(1,2,2);
contourf(gy1,-gz1,c1);
axis equal
xlabel('X');ylabel('Y')
curlmax=max(max(c1));
title(sprintf('|curl| < %10.2e', curlmax));
