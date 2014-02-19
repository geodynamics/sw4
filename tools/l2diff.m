%
% L2DIFF(ux1,uy1,uz1,st1,ux2,uy2,uz2,st2)
%
%
function diff=l2diff(ux1,uy1,uz1,st1,ux2,uy2,uz2,st2)
n1=length(ux1);
n2=length(ux2);
diff=-999;
if (st1 <= 0 || st2 <= 0)
  disp(['ERROR st1, st2 must be positive, not: ', num2str(st1), num2str(st2)]);
  return
end
if (length(ux1(1:st1:n1)) ~= length(ux2(1:st2:n2)))
  disp(['ERROR: ux1 with stride st1=' num2str(st1) ' length=' num2str(length(ux1(1:st1:n1))) ]);
  disp(['ux2 with stride st2=', num2str(st2), ' length=', num2str(length(ux2(1:st2:n2)))]);
  return
end
ns = length(ux1(1:st1:n1));
diff = sum((ux1(1:st1:n1) - ux2(1:st2:n2)).^2 + (uy1(1:st1:n1) - uy2(1:st2:n2)).^2 + (uz1(1:st1:n1) - uz2(1:st2:n2)).^2);
diff = sqrt(diff/ns);