
% Plot comparison between Prose and FD at
% station 10 of the LOH.1 problem with
% Qp=Qs=inf.

%  28Jun00

sig=0.06;
% Filename for exact solution
filename='LOH.1_prose3';
ReadUHS

indx=(1:2000);

figure(1)
clf
plot(t(indx),ra(indx),'b','LineWidth',1)
set(gca,'FontSize',14);
title('Semi-analytical solution of LOH.1, station 10. Radial component');
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);

figure(2)
clf
plot(t(indx),tr(indx),'b','LineWidth',1)
set(gca,'FontSize',14);
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);
title('Semi-analytical solution of LOH.1, station 10. Transverse component');

figure(3)
clf
plot(t(indx),ve(indx),'b','LineWidth',1)
set(gca,'FontSize',14);
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
title('Semi-analytical solution of LOH.1, station 10. Vertical component');
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);


