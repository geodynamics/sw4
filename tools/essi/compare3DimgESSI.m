
% Script for comparing image data, essi data, and SAC files

% Open the 3D img files
dir = '/Users/hans/repos/sw4/tools/essi/M5.5_ESSI_srf.sw4output/';
file = 'm5.5_ESSI_srf.cycle=100.ux.3Dimg';
[im,x,y,z,t,timestring]=readimage3d(strcat(dir,file),1,1);

% x=4000, y=6000 -> i,j = 101,151
img.ux_11 = im(101,151,1);
img.ux_00 = im(96,146,1);
img.ux_02 = im(96,156,1);
img.ux_20 = im(106,146,1);
img.ux_22 = im(106,156,1);

% Open the 3D img files
dir = '/Users/hans/repos/sw4/tools/essi/M5.5_ESSI_srf.sw4output/';
file = 'm5.5_ESSI_srf.cycle=100.uy.3Dimg';
[im,x,y,z,t,timestring]=readimage3d(strcat(dir,file),1,1);

% x=4000, y=6000 -> i,j = 101,151
img.uy_11 = im(101,151,1);
img.uy_00 = im(96,146,1);
img.uy_02 = im(96,156,1);
img.uy_20 = im(106,146,1);
img.uy_22 = im(106,156,1);

% Op    en the 3D img files
dir = '/Users/hans/repos/sw4/tools/essi/M5.5_ESSI_srf.sw4output/';
file = 'm5.5_ESSI_srf.cycle=100.uz.3Dimg';
[im,x,y,z,t,timestring]=readimage3d(strcat(dir,file),1,1);

% x=4000, y=6000 -> i,j = 101,151
img.uz_11 = im(101,151,1);
img.uz_00 = im(96,146,1);
img.uz_02 = im(96,156,1);
img.uz_20 = im(106,146,1);
img.uz_22 = im(106,156,1);

% Have 5 station data, only compare the same time step
% S_11, x=4000, y=6000
file = 'S_11.x';
[u, dt, lat, lon, b, e, npts, year, jday, ...
    hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] ... 
    = readsac(strcat(dir,file));
sac.ux_11 = u(npts);

file = 'S_11.y';
[u, dt, lat, lon, b, e, npts, year, jday, ...
    hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] ... 
    = readsac(strcat(dir,file));
% Only compare the same time step
sac.uy_11 = u(npts);

file = 'S_11.z';
[u, dt, lat, lon, b, e, npts, year, jday, ...
    hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] ... 
    = readsac(strcat(dir,file));
% Only compare the same time step
sac.uz_11 = u(npts);

% For some reason, SAC files are single precision so only check 7 digits
assert(abs(sac.ux_11 - img.ux_11) < 1e-7 * max(abs(sac.ux_11),abs(img.ux_11)) );
assert(abs(sac.uy_11 - img.uy_11) < 1e-7 * max(abs(sac.uy_11),abs(img.uz_11)) );
assert(abs(sac.uz_11 - img.uz_11) < 1e-7 * max(abs(sac.uy_11),abs(img.uz_11)) );

% TODO - read in ESSI files using hdf5, compare results
file = 'm5.5_ESSI_srf.cycle=100.essi';
vel_0 = h5read(strcat(dir,file),'/vel_0 ijk layout');
vel_1 = h5read(strcat(dir,file),'/vel_1 ijk layout');
vel_2 = h5read(strcat(dir,file),'/vel_2 ijk layout');
% NB: matlab indices are reversed cycle,k,j,i
essi.ux_00 = vel_0(1,1,1,1);
essi.uy_00 = vel_1(1,1,1,1);
essi.uz_00 = vel_2(1,1,1);
essi.ux_20 = vel_0(1,1,1,11);
essi.uy_20 = vel_1(1,1,1,11);
essi.uz_20 = vel_2(1,1,1,11);
essi.ux_11 = vel_0(1,1,6,6);
essi.uy_11 = vel_1(1,1,6,6);
essi.uz_11 = vel_2(1,1,6,6);

assert(abs(essi.ux_11 - img.ux_11) < 1e-7 * max(abs(essi.ux_11),abs(img.ux_11)) );
assert(abs(essi.uy_11 - img.uy_11) < 1e-7 * max(abs(essi.uy_11),abs(img.uz_11)) );
assert(abs(essi.uz_11 - img.uz_11) < 1e-7 * max(abs(essi.uy_11),abs(img.uz_11)) );

