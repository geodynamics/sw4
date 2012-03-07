
% Plot comparison between Prose and FD at
% station 10 of the LOH.1 problem with 
% Qp=Qs=inf.

%  28Jun00

% Set path
%addpath Macintosh HD:MatlabScripts:PGE

sig=0.06;
Station=10;
% Files holding numerical solutions
% LOH.2
filename='cmugreen.loh.2';
Deconsave; % note that this routine uses the 4 first characters to create variables which are copied below
% makes cmug_time, cmug_rad, cmug_trans, cmug_vert

figure(2)

figure(1)
clf
plot(cmug_time,cmug_rad,'b-')
set(gca,'FontSize',14);
yMin=min(cmug_rad);
yMax=max(cmug_rad);
axis([0, 9, yMin - 0.1*(yMax-yMin), yMax + 0.1*(yMax-yMin)])
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);
title('Semi-analytical solution of LOH.2, station 10. Radial velocity component')

figure(2)
clf
plot(cmug_time,cmug_trans,'b-')
set(gca,'FontSize',14);
yMin=min(cmug_trans);
yMax=max(cmug_trans);
axis([0, 9, yMin - 0.1*(yMax-yMin), yMax + 0.1*(yMax-yMin)])
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);
title('Semi-analytical solution of LOH.2, station 10. Transverse velocity component')

figure(3)
clf
plot(cmug_time,cmug_vert,'b-')
set(gca,'FontSize',14);
yMin=min(cmug_vert);
yMax=max(cmug_vert);
axis([0, 9, yMin - 0.1*(yMax-yMin), yMax + 0.1*(yMax-yMin)])
h1=xlabel('Time (sec)');
set(h1,'FontSize',18);
h2=ylabel('Velocity (m/s)');
set(h2,'FontSize',18);
title('Semi-analytical solution of LOH.2, station 10. Vertical velocity component')
