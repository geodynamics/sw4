%
% Filter and scale Geodyn 2D (r,z) data. 
% Data is obtained by first calling routine readgeodyndata.
%
% Data is ramped smoothly to zero at the first and last times, and
% padding depending on the filter frequency is added before the first
% time and afte the last time.
% 
%   [tf,vrf,vzf,ruf,zuf]=filtergeodyn(tim,rmat,zmat,vr,vz,freq,dt);
%
%      Input : tim    - Time variable vector
%              rmat, zmat - (r,z) coordinates of input data
%                            3D arrays with  rmat(i,j,n) - r at
%                            position (i,j) at time tim(n).
%              vr, vz - r- and z- direction velocities at
%                       position (i,j) at time tim(n).
%              freq   - Filter data to this frequency
%              dt     - Return data on uniformly spaced times t=n*dt
%
%      Output: tf, vrf, vzf, ruf, zuf - Input data variables, 
%               restricted to uniform times dt*n and 
%               filtered to frequency freq.
%              toff - Time offset due to padding of time array.
function [tf,vrf,vzf,ruf,zuf,toff]=filtergeodyn(tim,rmat,zmat,vr,vz,freq,dt);

% Number of steps on file
[M1 M2 nmax]=size(vr);

% New times for filtered data
tfinal = tim(nmax);
nsteps = floor(tfinal/dt);
tunif  = (0:nsteps)*dt;

% Pad to prepare for filter
padlen = 6/freq;
npad = round(padlen/dt);

% Ramp down data smoothly over q grid points:
q = 8;
for j=1:M2
  for i=1:M1
     vr1(1:nmax) = vr(i,j,:);
     vz1(1:nmax) = vz(i,j,:);
     r1(1:nmax) = rmat(i,j,:);
     z1(1:nmax) = zmat(i,j,:);

     vru = interp1(tim,vr1,tunif );
     vzu = interp1(tim,vz1,tunif );

     ru = interp1(tim,r1,tunif);
     zu = interp1(tim,z1,tunif);

     for n=nsteps-q+2:nsteps+1
        arg = (nsteps+1-n)/(q-1);
        wgh = arg^4*(35-84*arg+70*arg*arg-20*arg*arg*arg);
        vru(n) = vru(n)*wgh;
	vzu(n) = vzu(n)*wgh;
     end;
     vrup(1:npad) = zeros(size(1:npad));
     vrup(npad+1:nsteps+npad+1) = vru;
     vrup(nsteps+npad+2:nsteps+2*npad+1) = zeros(size(1:npad));

     vzup(1:npad) = zeros(size(1:npad));
     vzup(npad+1:nsteps+npad+1) = vzu;
     vzup(nsteps+npad+2:nsteps+2*npad+1) = zeros(size(1:npad));

     rup(1:npad) = ru(1)*ones(size(1:npad));
     rup(npad+1:nsteps+npad+1) = ru;
     rup(nsteps+npad+2:nsteps+2*npad+1) = ru(nsteps+1)*ones(size(1:npad));

     zup(1:npad) = zu(1)*ones(size(1:npad));
     zup(npad+1:nsteps+npad+1) = zu;
     zup(nsteps+npad+2:nsteps+2*npad+1) = zu(nsteps+1)*ones(size(1:npad));

     [b,a]=mybutter2(freq*dt*2);
     vruf=myfiltfilt(b,a,vrup);
     vzuf=myfiltfilt(b,a,vzup);

     nsize = nsteps+2*npad+1;
     vrf(i,j,1:nsize) = vruf(1:nsize);
     vzf(i,j,1:nsize) = vzuf(1:nsize);
     tf(1:nsize) = dt*((1:nsize)-(npad+1));
     ruf(i,j,1:nsize) = rup(1:nsize);
     zuf(i,j,1:nsize) = zup(1:nsize);
  end;
end;
toff = -npad*dt;
disp(['Number of time steps = ' num2str(nsteps+2*npad+1)]);
disp(['time offset = ' num2str(tf(1))]);
