%-*-octave-*--
% PLOTCURVI: Plot 1-D curvilinear mapping
%               
% function plotgrid(  )
%
function plotcurvi( x, z, pstr)

  if nargin < 3
    pstr = "b";
  end
  
  plot(x,z,pstr,x',z',pstr);axis ij

end




