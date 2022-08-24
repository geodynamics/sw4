%-*-octave-*--
% PLOTCART: Plot Cartesian grid
%               
% function plotgrid(  )
%
function plotcart( x, z, pstr )
  if nargin < 3
    pstr = "b";
  end

  nx = size(x,2);
  nz = size(z,2);

  xv = [x; x];
  zv = [z(1)*ones(1,nx); z(nz)*ones(1,nx)];
				%
  xh = [x(1)*ones(1,nz); x(nx)*ones(1,nz)];
  zh = [z; z];

  plot(xv, zv, pstr, xh, zh, pstr); axis ij

end
         
