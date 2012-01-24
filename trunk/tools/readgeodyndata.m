%
% Read 2D Geodyn data files into Matlab, scales the data into SI
% units, and also removes 'kinks' due to restarts in the Geodyn
% computation. Geodyn data consists of a number of files, where
% each files represents the time history at one point in space.
% The files are named name_i_j where i and j are spatial 
% coordinates, with i the r-direction and j-the z-direction.
% The files are assumed to form an almost uniform Cartesian grid in (i,j).
%
%     [t,ur,uz,r,z,rho,cl,ct,pl]=readgeodyndata( n, m, dir, file )
%
%          Input: n, m - Size of geodyn file grid.
%                 dir  - Directory of the geodyn files (with '/'
%                 appended).
%                 file - name of files (without the '_i_j' part).
%
%          Output: t  - The points in time, 1D array (assumed to be
%                       the same at all points in the file grid).
%                  ur, uz - r and z velocity components, 3D arrays
%                             ur(i,j,n) where i is r-direction, j
%                             is z-direction, and n is the time
%                             variable.
%                  r, z   - r and z coordinates of the Geodyn file
%                  grid, 3D arrays.
%                  rho, cl, ct - Density, shear wave speed and
%                  pressure wave speed, 3D arrays.
%                  pl - Plastic strain, 3D array.
%
function [t,ur,uz,r,z,rho,cl,ct,pl]=readgeodyndata( n, m, dir, file )

if nargin < 4
   file = 'WPP122';
end;
if nargin < 3
   dir  = './';
end;

% Assume time is the same everywhere
w = load([dir file '_1_1']);
tr = w(:,1);
[n1 n2]=size(w);

% Check if time is monotone...
timemonotone = 1;
i = 1;
while i<n1 & timemonotone == 1
   if tr(i+1)<tr(i)
      timemonotone = 0;
   end;
   i = i+1;
end;

% ..if not, filter out restarts:
if timemonotone == 0
   [il, iu]=olegremoverestarts(tr);
   rilen=length(il);
   il(rilen+1) = n1+1;
else;
  rilen = 0;
  il(1) = n1+1;
end;

for j=1:m
  for i=1:n
    w=load([dir file '_' num2str(i) '_' num2str(j)]);

    t(1:il(1)-1)      = tr(1:il(1)-1);
    ur(i,j,1:il(1)-1) = w(1:il(1)-1,5);
    uz(i,j,1:il(1)-1) = w(1:il(1)-1,6);
    r(i,j,1:il(1)-1)  = w(1:il(1)-1,8);
    z(i,j,1:il(1)-1)  = w(1:il(1)-1,9);
    rho(i,j,1:il(1)-1)= w(1:il(1)-1,4);
    cl(i,j,1:il(1)-1) = w(1:il(1)-1,2);
    ct(i,j,1:il(1)-1) = w(1:il(1)-1,3);
    pl(i,j,1:il(1)-1) = w(1:il(1)-1,7);

% If need to remove restarts, continue here:
    pos = il(1);
    for l=1:rilen
       for k = iu(l)+1:il(l+1)-1
	 t(pos) = tr(k);
	 ur(i,j,pos) = w(k,5);
	 uz(i,j,pos) = w(k,6);
         r(i,j,pos)  = w(k,8);
	 z(i,j,pos)  = w(k,9);
	 rho(i,j,pos)= w(k,4);
	 cl(i,j,pos) = w(k,2);
	 ct(i,j,pos) = w(k,3);	 
	 pl(i,j,pos) = w(k,7);
         pos = pos+1;
       end;
    end;
  end;
end;

% Transform to SI units
% Geodyn uses length scale = 1 mm
% Velocity scale = 1 km/s
% Time scale = 1 us
z  = z*1e-3;
r  = r*1e-3;
ur = ur*1e3;
uz = uz*1e3;
t  = t*1e-6;
