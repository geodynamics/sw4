%
% HUDSONPLOT
%
% Plot as described by Hudson, Pearce, and Rogers in J.Geophys.Res. vol 94 (1989), pp.765-774.
%
%      function hudsonplot( mxx, mxy, mxz, myy, myz, mzz )
%
%      Input: mxx, mxy, mxz, myy, myz, mzz - Components of the moment tensor
%
function hudsonplot( mxx, mxy, mxz, myy, myz, mzz )
clrplot=0;
M = [mxx mxy mxz; mxy myy myz;mxz myz mzz];
[u,d]=eig(M);
mvals = sort([d(1,1) d(2,2) d(3,3)]);
mx=mvals(3);
mz=mvals(2);
my=mvals(1);
mavg=(mx+my+mz)/3;
mxp=mx-mavg;
myp=my-mavg;
mzp=mz-mavg;
if mzp > 0
  k = mavg/(abs(mavg)-myp);
  t = -2*mzp/myp;
elseif mzp == 0
  t = 0;
  k = mavg/(abs(mavg)+mxp);
else
  k = mavg/(abs(mavg)+mxp);
  t= 2*mzp/mxp;
end;
tau=t*(1-abs(k));

if tau*k < 0 
% 2nd or 4th quadrant:
  u=tau;
  v=k;
elseif tau >=0 && k >=0 
% 1st quadrant:
   if tau < 4*k
      u = tau/(1-tau*0.5);
      v = k/(1-tau*0.5);
   else
      u = tau/(1-2*k);	 
      v = k/(1-2*k);
   end;
else
% 3rd quadrant
   if tau > 4*k
      u = tau/(1+tau*0.5);
      v = k/(1+tau*0.5);
   else
      u = tau/(1+2*k);	 
      v = k/(1+2*k);
   end;
end;
if clrplot == 1
   clf;
end;
% outline boundaries in the u-v plane
bndcol='k';
diagcol='b';
set(gca,'FontSize',25);
h=plot([-4/3 0 4/3],[-1/3 1 1/3],bndcol);
set(h,'LineWidth',2.0);
hold on;
h=plot([-4/3 0 4/3],[-1/3 -1 1/3],bndcol);
set(h,'LineWidth',2.0);
% diagonal separating regions A and B
h=plot([-4/3 4/3],[-1/3 1/3],diagcol);
set(h,'LineWidth',2.0);
% Coordinate axis
h=plot([-1.3 1.3],[0 0],'k');
set(h,'LineWidth',1.0);
h=plot([0 0],[-1.1 1.1],'k');
set(h,'LineWidth',1.0);
% plot source
h=plot([u],[v],'*r','markersize',10);
h=plot([u],[v],'or','markersize',10);
xlabel('u');
ylabel('v');
axis('equal');
axis([-1.4 1.4 -1.2 1.2]);
hold off;

