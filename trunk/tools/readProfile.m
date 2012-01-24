%
% readProfile
%
%  Reads a vs profile binary file and outputs the grid spacing required at
%  the top and bottom of the domain (z=0.0 and z=maxz)
%
%        [vs, dh] = vs2blocks( fil )
%
%        Input:  fil   - Name of binary profile file.
%        Output: vs    - Average vs value in z direction
%                dh    - Original grid spacing
%                freq  - Original frequency from WPP source

function [vs, dh, freq]=readProfile(fil)

  fd=fopen(fil,'r');
   if fd < 0
      disp( ['Error: could not open file ' fil] );
   else
     %   Read the file
    dh = fread(fd, 1, 'double')
    freq = fread(fd, 1, 'double')
    nz= fread(fd, 1, 'int')
    %disp(['Grid Spacing: ' num2str(dh)])
    %disp(['Frequency: ' num2str(freq)])
    %disp(['Number points: ' num2str(nz)])
    vs = fread(fd, nz,'double');

    plot(vs);

    % PPW (points per wavelength) constant
    ppw = 10

    %  Now let's find the coarse grid spacing
    coarseDH = vs(nz)/(ppw*freq);
    fineDH = vs(1)/(ppw*freq);

    disp(['Coarse dh: ' num2str(coarseDH)]);
    disp(['Fine   dh: ' num2str(fineDH)]);
   end
   fclose(fd);
return