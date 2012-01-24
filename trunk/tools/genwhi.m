%
% Generate a whi-file from filtered geodyn 2D (r,z) data.
% 
%   genwhi( filename, tf, rf, zf, vrf, vzf, L, D, hbox, freq, toff )
%
%      Input : filename - Name of output file.
%              tf       - Time variable vector, uniform times of
%              filtered data
%              rf, zf - (r,z) coordinates of filtered data
%                            3D arrays with  rf(i,j,n) = r at
%                            position (i,j) and time tf(n).
%              vrf, vzf - r- and z- direction velocities of filtered
%                    data at position (i,j) at time tf(n).
%              L, D - Dimension of output box is [-L/2,L/2]^2x[0,D]
%              hbox - Spacing of data on box.
%              freq - Filter frequency used for tf,rf,zf,vrf,vzf
%              toff - Time offset due to filter padding. 
%
function genwhi( filename, tf, rf, zf, vrf, vzf, ...
				      L, D, hbox, freq, toff )

% Number of points in file
N = L/hbox+1;
Nz= D/hbox+1;

% Source center
xs = 0;
ys = 0;
zs = 122;

% Granite
rho = 2640;
vs  = 3120;
vp  = 5090;

dt = tf(2)-tf(1);
nsteps = length(tf);
% Write file header
fid = fopen(filename,'wt');
fprintf(fid,'external code=geodyn2d\n');
% Set origin to (L,L) assumes domain of size (3*L,3*L).
fprintf(fid,['grid faces=6 stepsize=%.10g nx=%i ny=%i nz=%i x0=%i y0=%i z0=0.0\n'],hbox,N,N,Nz,L,L);
fprintf(fid,'time timestep=%.10g nsteps=%i toff=%.10g\n',dt, nsteps,toff );

fprintf(fid,'material rho=%.10g vs=%.10g vp=%.10g \n',rho,vs,vp);
fprintf(fid,'source filter=%.10g x0=%.10g y0=%.10g z0=%.10g\n',freq,xs,ys,zs);
fprintf(fid,'begindata\n');

% Interpolate
for n=1:nsteps
%for n=1:nsteps+2*npad+1
   intpplane( rf(:,1,n), zf(1,:,n), vrf(:,:,n), vzf(:,:,n), L, D, ...
	      hbox, fid, 0 );
end;
