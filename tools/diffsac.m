%
% DIFFSAC
%
%  Read 2 receiever data files in SAC format and compare/diff
%
%              diffsac(basename1, basename2)
%
%       Input: basename1,2 - Name of receiever file, basename.x, basename.y, basename.z
%              dir1, dir2 - Path to the receiver files
%              
%               
function diffsac(pathbasename1, pathbasename2)

if nargin < 2
  error('Need 2 arguments');
end;
file(1:2,:) = [pathbasename1; pathbasename2];
sdir = ['x','y','z'];
names = {'u', 'dt', 'lat', 'lon', 'b', 'e', 'npts', 'year', 'jday', ...
    'hour', 'min', 'sec', 'msec', 'cmpaz', 'cmpinc', 'idep', 'stnam'}; 

% Stuff all the data into arrays of structures
for f = 1:2
for dir = 1:3
    d = struct([]);
    msg = sprintf('Reading basefile %s.%s:', file(f,:), sdir(1));
    disp(msg);
    [u, dt, lat, lon, b, e, npts, year, jday, hour, min, sec, msec, cmpaz, cmpinc, idep, stnam ] ...
        =readsac(sprintf('%s.%s', file(f,:), sdir(dir)));
    disp(' ');
    d(1).u=u; d(1).dt=dt; d(1).lat=lat; d(1).lon=lon;
    d(1).b=b; d(1).e=e; d(1).npts=npts; d(1).year=year; d(1).jday=jday; 
    d(1).hour=hour; d(1).min=min; d(1).sec=sec; d(1).msec=msec;
    d(1).cmpaz=cmpaz; d(1).cmpinc=cmpinc; d(1).idep=idep; d(1).stnam=stnam;
    data(f,dir) = d(1);
end
end

% Compare fields within same recording files
for f=1:2
for dir=2:3
    data1 = data(f,1);
    data2 = data(f,dir);
    msg = sprintf('Between basefiles %s.%s and .%s:', ...
        file(f,:), sdir(1), sdir(dir));
    disp(msg);
    for ival = 2:size(names,2)-1
        val1 = getfield(data1,names{ival});
        val2 = getfield(data2,names{ival});
        if (val1 ~= val2)
            msg = sprintf('  field %s differs:\n', names{ival});
            msg = [msg, '    ', num2str(val1), ' vs. ', num2str(val2)];
            disp(msg);
        end
    end
    disp(' ');
end
end

% Compare fields within same recording files
for dir=1:3
    data1 = data(1,dir);
    data2 = data(2,dir);
    msg = sprintf('Between basefiles %s.%s and %s.%s:', ...
        pathbasename1, sdir(dir), pathbasename2,sdir(dir));
    disp(msg);
    for ival = 2:size(names,2)-1
        val1 = getfield(data1,names{ival});
        val2 = getfield(data2,names{ival});
        if (val1 ~= val2)
            msg = sprintf('  field %s differs:\n', names{ival});
            msg = [msg, '    ', num2str(val1), ' vs. ', num2str(val2)];
            disp(msg);
        end
    end
    % Check the u field
    u1 = getfield(data1,names{1});
    u2 = getfield(data2,names{1});
    udiff = norm(u1-u2,inf);
    msg = sprintf('  u%s field max difference: %1.1e\n', sdir(dir), udiff);
    disp(msg);

end

end
