% Test for ESSI hdf5 output in SW4
% Assumes you have run the tests in the sw4/tool/essi directory:
%   sw4 M5.5_ESSI_srf.in (produces M5.5_ESSI_srf/ files)
%   sw4 berkeley-att-h100.in (produces berkeley-att-h100/ files)
clear;

basedir = pwd;

% ---------------------------------------------------------------------
% M5.5_ESSI_srf test
% ---------------------------------------------------------------------
testname = 'M5.5_ESSI_srf';
fprintf('%%% Beginning %s test\n', testname);

cycle = 10; % cycle at which we'll compare output
cycle_str = sprintf('%02d',cycle); % string cycle number
file = sprintf('%s/M5.5_ESSI_srf/m5.5_ESSI_srf.cycle=%s.essi', basedir, cycle_str);
if (exist(file))
    cycle_start_end = h5read(file,'/cycle start, end');
    assert(cycle >= cycle_start_end(1) && cycle <= cycle_start_end(2));
    cycle_ix = cycle+1; % 1-based indexing
    comp_str=['vel_0'; 'vel_1'; 'vel_2'];
    for comp=1:3
        vel = h5read(file,sprintf('/%s ijk layout',comp_str(comp,:)));
        % NB: matlab indices are reversed cycle,k,j,i
        index = cycle_ix - cycle_start_end(1);
        if (comp == 1)
            essi.ux_00 = vel(index,1,1,1);
            essi.ux_11 = vel(index,1,6,6);
            essi.ux_20 = vel(index,1,1,11);
        elseif (comp == 2)
            essi.uy_00 = vel(index,1,1,1);
            essi.uy_11 = vel(index,1,6,6);
            essi.uy_20 = vel(index,1,1,11);
        elseif (comp == 3)
            essi.uz_00 = vel(index,1,1,1);
            essi.uz_11 = vel(index,1,6,6);
            essi.uz_20 = vel(index,1,1,11);
        end
    end

    % Read in the appropriate 3Dimg file
    comp_str=['ux'; 'uy'; 'uz'];
    for comp=1:3
        file = sprintf('%s/M5.5_ESSI_srf/m5.5_ESSI_srf.cycle=%s.%s.3Dimg', ...
            basedir,cycle_str,comp_str(comp,:));
        if (~exist(file))
            warning('Unable to open .3Dimg file: %s', file); 
            clear img;
            break;
        end
        [im,x,y,z,t,timestring]=readimage3d(file,1,1);

        % x=4000, y=6000 -> i,j = 101,151
        if (comp == 1)
            img.ux_11 = im(101,151,1);
            img.ux_00 = im(96,146,1);
            img.ux_02 = im(96,156,1);
            img.ux_20 = im(106,146,1);
            img.ux_22 = im(106,156,1);
        elseif (comp == 2)
            img.uy_11 = im(101,151,1);
            img.uy_00 = im(96,146,1);
            img.uy_02 = im(96,156,1);
            img.uy_20 = im(106,146,1);
            img.uy_22 = im(106,156,1);
        elseif (comp == 3)
            img.uz_11 = im(101,151,1);
            img.uz_00 = im(96,146,1);
            img.uz_02 = im(96,156,1);
            img.uz_20 = im(106,146,1);
            img.uz_22 = im(106,156,1);
        end
    end

    % Have 5 station data, only compare the same time step
    % S_11, x=4000, y=6000
    comp_str=['x'; 'y'; 'z'];
    for comp=1:3
        file = sprintf('%s/M5.5_ESSI_srf/S_11.%s', basedir,comp_str(comp,:));
        if (~exist(file))
            warning('Unable to open rec file: %s', file);
            clear sac;
            break;
        end

        [u, dt, lat, lon, b, e, npts, year, jday, ...
            hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] ... 
            = readsac(file);
        if (comp == 1)
            sac.ux_11 = u(cycle_ix);
        elseif (comp == 2)
            sac.uy_11 = u(cycle_ix);
        elseif (comp == 3)
            sac.uz_11 = u(cycle_ix);
        end
    end

    % SAC files are single precision so only check 7 digits
    pass = false;
    if (exist('sac'))
        assert(abs(essi.ux_11 - sac.ux_11) < 1e-7 * max(abs(essi.ux_11),abs(sac.ux_11)) );
        assert(abs(essi.uy_11 - sac.uy_11) < 1e-7 * max(abs(essi.uy_11),abs(sac.uz_11)) );
        assert(abs(essi.uz_11 - sac.uz_11) < 1e-7 * max(abs(essi.uy_11),abs(sac.uz_11)) );
        fprintf('--> %s spot check of SAC and ESSI passes.\n', testname);
        pass = true;
    else
        disp('Skipping comparison of SAC and ESSI output.');
    end

    if (exist('img'))
        assert(abs(essi.ux_11 - img.ux_11) < 1e-14 * max(abs(essi.ux_11),abs(img.ux_11)) );
        assert(abs(essi.uy_11 - img.uy_11) < 1e-14 * max(abs(essi.uy_11),abs(img.uz_11)) );
        assert(abs(essi.uz_11 - img.uz_11) < 1e-14 * max(abs(essi.uy_11),abs(img.uz_11)) );
        fprintf('--> %s spot check of 3DIMG and ESSI passes.\n', testname);
        pass = true;
    else
        disp('Skipping comparison of 3DIMG and ESSI output.');
    end
    if (~pass)
        error('--> Did not pass/complete one or more %s tests.', testname);
    end
else
    warning('Skipping M5.5_ESSI_srf test, cannot open hdf5 file: %s', file);
end


% ---------------------------------------------------------------------
% berkeley-att-h100 test
% ---------------------------------------------------------------------
clear;
testname = 'berkeley-att-h100';
fprintf('%%% Beginning %s test\n', testname);

basedir = pwd;
cycle = 20; % cycle at which we'll compare output
cycle_str = sprintf('%02d',cycle); % string cycle number
% ESSI output for this run is all in one filename with cycle=00
file = sprintf('%s/berkeley-att-h100/berkeley.cycle=%s.essi', basedir, '00');
if (exist(file))
    % This has topography, so read that first
    z_topo = h5read(file,'/z coordinates');
    essi.z_000 = z_topo(1,1,1); % Note indices are local (k,j,i)
    essi.z_110 = z_topo(1,3,3);
    essi.z_200 = z_topo(1,1,5);
    essi.z_112 = z_topo(3,3,3);

    cycle_start_end = h5read(file,'/cycle start, end');
    assert(cycle >= cycle_start_end(1) && cycle <= cycle_start_end(2));
    cycle_ix = cycle+1; % 1-based indexing
    comp_str=['vel_0'; 'vel_1'; 'vel_2'];
    for comp=1:3
        vel = h5read(file,sprintf('/%s ijk layout',comp_str(comp,:))); 
        % NB: matlab indices are reversed cycle,k,j,i
        index = cycle_ix - cycle_start_end(1);
        if (comp == 1)
            essi.ux_000 = vel(index,1,1,1); % Indices are local (cycle,k,j,i)
            essi.ux_110 = vel(index,1,3,3);
            essi.ux_200 = vel(index,1,1,5);
            essi.ux_112 = vel(index,3,3,3);
        elseif (comp == 2)
            essi.uy_000 = vel(index,1,1,1);
            essi.uy_110 = vel(index,1,3,3);
            essi.uy_200 = vel(index,1,1,5);
            essi.uy_112 = vel(index,3,3,3);
        elseif (comp == 3)
            essi.uz_000 = vel(index,1,1,1);
            essi.uz_110 = vel(index,1,3,3);
            essi.uz_200 = vel(index,1,1,5);
            essi.uz_112 = vel(index,3,3,3);
        end
    end
    
    % Read in the appropriate 3Dimg file
    comp_str=['ux'; 'uy'; 'uz'];
    for comp=1:3
        file = sprintf('%s/berkeley-att-h100/berkeley.cycle=%s.%s.3Dimg', ...
            basedir,cycle_str,comp_str(comp,:));
        if (~exist(file))
            warning('Unable to open .3Dimg file: %s', file); 
            clear img;
            break;
        end
        [im,x,y,z,t,timestring]=readimage3d(file,2,1);
        img.z_110 = z(41,61,1); % Indices are global (i,j,k)
        img.z_000 = z(39,59,1);
        img.z_020 = z(39,63,1);
        img.z_200 = z(43,59,1);
        img.z_220 = z(43,63,1);
        img.z_112 = z(41,61,3);

        % x=.15,.3,.45 -> i = 9,18,27
        % y=.45,.5,.55 -> j = 27,30,33
        if (comp == 1)
            img.ux_110 = im(41,61,1);
            img.ux_000 = im(39,59,1);
            img.ux_020 = im(39,63,1);
            img.ux_200 = im(43,59,1);
            img.ux_220 = im(43,63,1);
            img.ux_112 = im(41,61,3);
        elseif (comp == 2)
            img.uy_110 = im(41,61,1);
            img.uy_000 = im(39,59,1);
            img.uy_020 = im(39,63,1);
            img.uy_200 = im(43,59,1);
            img.uy_220 = im(43,63,1);
            img.uy_112 = im(41,61,3);
        elseif (comp == 3)
            img.uz_110 = im(41,61,1);
            img.uz_000 = im(39,59,1);
            img.uz_020 = im(39,63,1);
            img.uz_200 = im(43,59,1);
            img.uz_220 = im(43,63,1);
            img.uz_112 = im(41,61,3);
        end
    end
    
    % Have 5 station data, only compare the same time step
    % S_11, x=.3, y=.5
    comp_str=['x'; 'y'; 'z'];
    for comp=1:3
        file = sprintf('%s/berkeley-att-h100/S_11.%s', basedir,comp_str(comp,:));
        if (~exist(file))
            warning('Unable to open rec file: %s', file);
            clear sac;
            break;
        end

        [u, dt, lat, lon, b, e, npts, year, jday, ...
            hour, min, sec, msec, cmpaz, cmpinc, idep, stnam] ... 
            = readsac(file);
        if (comp == 1)
            sac.ux_11 = u(cycle_ix);
        elseif (comp == 2)
            sac.uy_11 = u(cycle_ix);
        elseif (comp == 3)
            sac.uz_11 = u(cycle_ix);
        end
    end

    pass = false;
    if (exist('img'))
        assert(abs(essi.z_110 - img.z_110) < 1e-14 * max(abs(essi.z_110),abs(img.z_110)) );
        assert(abs(essi.z_112 - img.z_112) < 1e-14 * max(abs(essi.z_112),abs(img.z_112)) );
        fprintf('--> %s spot check of 3DIMG and ESSI z coordinate passes.\n', testname);

        assert(abs(essi.ux_110 - img.ux_110) < 1e-14 * max(abs(essi.ux_110),abs(img.ux_110)) );
        assert(abs(essi.uz_112 - img.uz_112) < 1e-14 * max(abs(essi.uz_112),abs(img.uz_112)) );
        fprintf('--> %s spot check of 3DIMG and ESSI velocities passes.\n', testname);
        pass = true;
    else
        fprintf('Skipping comparison of 3DIMG and ESSI output for %s.', testname);
    end
    
    if (exist('sac'))
        assert(abs(essi.ux_110 - sac.ux_11) < 1e-7 * max(abs(essi.ux_110),abs(sac.ux_11)) );
        assert(abs(essi.uy_110 - sac.uy_11) < 1e-7 * max(abs(essi.uy_110),abs(sac.uz_11)) );
        assert(abs(essi.uz_110 - sac.uz_11) < 1e-7 * max(abs(essi.uy_110),abs(sac.uz_11)) );
        fprintf('--> %s spot check of SAC and ESSI passes.\n', testname);
        pass = true;
    else
        disp('Skipping comparison of SAC and ESSI output.');
    end
    
    if (~pass)
        error('--> Did not pass/complete one or more %s tests.', testname);
    end


else
    warning('Skipping %s test, cannot open hdf5 file: %s', testname, file);
end

function foo = bar()
end

function foo2 = bar2()
end

