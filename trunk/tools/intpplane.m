%
% INTPPLANE
%
%   Interpolate data from 2D axisymmetric data to the six sides of
%   a cube in three dimensions. The axisymmetric domain is in
%   (r,z)-coordinates, the cube is given by [-L/2,-L/2]x[-L/2,L/2]x[0,D]. 
%   Note the cube z axis is downward, while the 2D axisymmetric is upward, i.e.,
%   z_cube = -z_2D. Spacing is assumed equal in all directions.
%   The nearest points in the (r,z) domain are used to impose data
%   on the cube, the result will not be accurate if the cube
%   extends outside the (r,z) domain.
%
%        [uxo,uyo,uzo]=intpplane(rc,zc,vr,vz,L,D,hbox,fid,sout)
%        
%        
%       Input: rc, zc - Vectors of coordinates (1D arrays)
%              vr, vz - Velocities on the 2D domain (2D arrays).
%              L      - Size of domain.
%              D      - Depth of domain.
%              hbox   - Grid spacing in the three dimensional cube
%              fid    - The cube data is output to the already
%                       opened file belonging to this file descriptor.
%              sout   - Return data on side = sout, if zero no data
%                       is returned.
%       Output: uxo, uyo, uzo - Data on side sout. 
function [uxo,uyo,uzo]=intpplane( rc, zc, vr, vz, L, D, hbox, ...
					   fid, sout )
if nargin < 9
  sout = 0;
end;

% Number of points on cube
N    = L/hbox+1;
Nz   = D/hbox+1;

% Approximate grid sizes, the grid might be slightly non-uniform.
[M1 M2]=size(vr);
Lz = zc(M2)-zc(1);
Lr = rc(M1)-rc(1);
hzapp = Lz/(M2-1);
hrapp = Lr/(M1-1);

% Check that cube is inside domain
rtol = 0.05*Lr;
ztol = 0.05*Lz;
if L/2 > Lr+rtol
  disp(['Warning: Cube r-dimension outside range of given data']);
end;
if -D < zc(1)-ztol | zc(M2)+ztol< 0
  disp(['Warning: Cube z-dimension outside range of given data']);
end;

uxo = 0;
uyo = 0;
uzo = 0;

for side=1:4
   if side == 1 
     i = 1;
     x = (i-1)*hbox-L/2;
   elseif side == 2
     i = N;
     x = (i-1)*hbox-L/2;
   elseif side == 3
     j = 1;
     y = (j-1)*hbox-L/2;
   elseif side == 4 
     j = N;
     y = (j-1)*hbox-L/2;
   end;
   
   if side == 1 
      for k=1:Nz
	z = -(k-1)*hbox;
	kind = round((z+D)/hzapp+1);
	if z > zc(kind) 
	  while z > zc(kind) & kind < M2
	    kind = kind +1;
	  end;
	  kind = kind-1;
	else
	  while z <= zc(kind) & kind > 1
	    kind = kind-1;
	  end;
	end;
	wghz(k) = (zc(kind+1)-z)/(zc(kind+1)-zc(kind));
	kk(k)   = kind;
      end;
   end;

   for j=1:N
     if side == 1 | side == 2 
       y = (j-1)*hbox-L/2;
     else
       x = (j-1)*hbox-L/2;
     end;
     r = sqrt(x*x+y*y);
     jind = round(r/hrapp+1);
     if r > rc(jind) 
       while r > rc(jind) & jind < M1
	 jind = jind +1;
       end;
       jind = jind-1;
     else
       while r <= rc(jind) & jind > 1
	 jind = jind-1;
       end;
     end;
     wghr(j) = (rc(jind+1)-r)/(rc(jind+1)-rc(jind));
     rr(j) = jind;
     ca(j) = x/r;
     sa(j) = y/r;
   end;

   for k=1:Nz
     for j=1:N
	 vrint = wghz(k)*(wghr(j)*vr(rr(j),kk(k)) + (1-wghr(j))*vr(rr(j)+1,kk(k))) ...
	       +(1-wghz(k))*(wghr(j)*vr(rr(j),kk(k)+1) + (1-wghr(j))*vr(rr(j)+1,kk(k)+1));
	 vzint = wghz(k)*(wghr(j)*vz(rr(j),kk(k)) + (1-wghr(j))*vz(rr(j)+1,kk(k))) ...
	       +(1-wghz(k))*(wghr(j)*vz(rr(j),kk(k)+1) + (1-wghr(j))*vz(rr(j)+1,kk(k)+1));
	 ux(j,k) = vrint*ca(j);
	 uy(j,k) = vrint*sa(j);
	 uz(j,k) = -vzint;
	 fprintf(fid,'%12.5g %12.5g %12.5g\n',ux(j,k),uy(j,k),uz(j,k));
     end;
   end;
   if side == sout 
     uxo = ux;
     uyo = uy;
     uzo = uz;
   end;
end;

for side=5:6
  if side == 5 
    k = M2;
  elseif side == 6 
    k = 1;
  end;
  for j=1:N
    y = (j-1)*hbox-L/2;
    for i=1:N
      x = (i-1)*hbox-L/2;
      r = sqrt(x*x+y*y);
      jind = round(r/hrapp+1);
      if r > rc(jind) 
	while r > rc(jind) & jind < M1
	  jind = jind +1;
	end;
	jind = jind-1;
      else
	while r <= rc(jind) & jind > 1
	  jind = jind-1;
	end;
      end;
      wghr = (rc(jind+1)-r)/(rc(jind+1)-rc(jind));
      rr = jind;
      if abs(r)<1e-6
        ca = 1;
	sa = 0;
      else
	ca = x/r;
	sa = y/r;
      end;
      vrint   = (wghr*vr(rr,k) + (1-wghr)*vr(rr+1,k));
      vzint   = (wghr*vz(rr,k) + (1-wghr)*vz(rr+1,k));

      ux(i,j) = vrint*ca;
      uy(i,j) = vrint*sa;
      uz(i,j) = -vzint;
      fprintf(fid,'%12.5g %12.5g %12.5g\n',ux(i,j),uy(i,j),uz(i,j));
    end;
   end;
   if side == sout 
     uxo = ux;
     uyo = uy;
     uzo = uz;
   end;
end;

