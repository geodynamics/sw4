%-*-octave-*--
% CURVIGRID: Plot 1-D curvilinear mapping
%               
% function plotgrid(  )
%
function plotgrid( x, z, pstr)

  if nargin < 3
    pstr = "b";
  end
  
  plot(x,z,pstr);axis ij
  hold on;
  i = 10;
  n1 = size(z,1);
  x1 = [x', x']';
  z1 = [z(1,:)', z(n1,:)']';
%  plot([x(i), x(i)], [z(1,i), z(n1,i)], pstr);
  plot(x1, z1, pstr);
  hold off;

end




