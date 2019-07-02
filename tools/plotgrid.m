%-*-octave-*--
%
% PLOTGRID
%
%     plotgrid( fil, machineformat, cvals )
%
%   Plots the image on file 'fil' with contour, using the contour levels cvals.
%   If cvals is ommitted, 21 countour levels will be obtained through the imageinfo fcn.
%
%   The boundary of each grid patch is outlined in black. A vector cvals can be 
%   obtained from function imageinfo.
%
%   Input:
%         fil:                       Name of image file
%         machineformat (optional):  Passed to fopen to read big endian, little endian, etc
%
function plotgrid( fil, machineformat )
  ## if nargin < 2
  ##   pstr = "b";
  ## end;
  if nargin < 2
    machineformat='native';
  end;

  col = ['k', 'r', 'g', 'b', 'm', 'c', 'y'];
  fd=fopen(fil,'r',machineformat);
  pr=fread(fd,1,'int');
  npatches=fread(fd,1,'int');
  t=fread(fd,1,'double');
  plane=fread(fd,1,'int');
  coord   =fread(fd,1,'double');
  mode    =fread(fd,1,'int');
  gridinfo=fread(fd,1,'int');
  fclose(fd);

  firstCurviPatch = npatches - gridinfo + 1;

  topCartesian = npatches - gridinfo;

  x1min=1e9;
  x2min=1e9;
  x1max=-1e9;
  x2max=-1e9;

  for b=1:topCartesian
    [im,x,y,z] = readimage(fil,b,0,machineformat);
    if plane==0
      plotcart(y,z,col(b));
      x1mi=min(y);
      x1ma=max(y);
      x2mi=min(z);
      x2ma=max(z);
    elseif plane==1
      plotcart(x,z,col(b));
      x1mi=min(x);
      x1ma=max(x);
      x2mi=min(z);
      x2ma=max(z);
    elseif plane==2
      plotcart(x,y,col(b));
      x1mi=min(x);
      x1ma=max(x);
      x2mi=min(y);
      x2ma=max(y);
    end
    if b==1
      hold on;
    end;

				% update global min/max
    x1min = min(x1min,x1mi);
    x1max = max(x1max,x1ma);
    x2min = min(x2min,x2mi);
    x2max = max(x2max,x2ma);
  end; % for

  for b=firstCurviPatch: npatches
    [im,x,y,z] = readimage(fil,b,0,machineformat);
    if plane==0
      plotcurvi(y,z,col(b));
      x1mi=min(min(y));
      x1ma=max(max(y));
      x2mi=min(min(z));
      x2ma=max(max(z));
    elseif plane==1
      plotcurvi(x,z,col(b));
      x1mi=min(min(x));
      x1ma=max(max(x));
      x2mi=min(min(z));
      x2ma=max(max(z));
    elseif plane==2
      plotcurvi(x,y,col(b));
      x1mi=min(x);
      x1ma=max(x);
      x2mi=min(y);
      x2ma=max(y);
    end
				% update global min/max
    x1min = min(x1min,x1mi);
    x1max = max(x1max,x1ma);
    x2min = min(x2min,x2mi);
    x2max = max(x2max,x2ma);
  end; % for curvilinear

  axis([x1min x1max x2min x2max]);
  hold off;
  axis ij; % flip z-axis to point downwards
  set(gca,"fontsize",16);
%  axis equal;
end %plotgrid

