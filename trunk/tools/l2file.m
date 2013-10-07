%
% L2FILE(fn1,stride1,fn2,stride2)
%
% Input:
%       fn1, fn2:    file names for USGS formatted time series
%       st1, st2:    strides to use when comparing the the two time series.
%
function diff=l2file(fn1,stride1,fn2,stride2)
[t1 ux1 uy1 uz1]=readusgs(fn1);
[t2 ux2 uy2 uz2]=readusgs(fn2);

n1=length(ux1);
n2=length(ux2);
diff=-999;
if (stride1 <= 0 || stride2 <= 0)
  disp(['ERROR stride1, stride2 must be positive, not: ', num2str(stride1), num2str(stride2)]);
  return
end
if (length(ux1(1:stride1:n1)) ~= length(ux2(1:stride2:n2)))
  disp(['ERROR: ux1 with stride=', num2str(st1) ' length=' num2str(length(ux1(1:stride1:n1))) ]);
  disp(['ux2 with stride=', num2str(st2), ' length=', num2str(length(ux2(1:stride2:n2)))]);
  return
end
ns = length(ux1(1:stride1:n1));
diff = sum((ux1(1:stride1:n1) - ux2(1:stride2:n2)).^2 + (uy1(1:stride1:n1) - uy2(1:stride2:n2)).^2 + (uz1(1:stride1:n1) - uz2(1:stride2:n2)).^2);
diff = sqrt(diff/ns);

% plot
tr=t1(1:stride1:n1);
subplot(3,1,1)
plot(tr,ux1(1:stride1:n1)-ux2(1:stride2:n2),"b")
legend("x-diff")
subplot(3,1,2)
plot(tr,uy1(1:stride1:n1)-uy2(1:stride2:n2),"g")
legend("y-diff")
subplot(3,1,3)
plot(tr,uz1(1:stride1:n1)-uz2(1:stride2:n2),"r")
legend("z-diff")
