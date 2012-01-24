%
% WPPGETIMAGE
%
%  Reads one image slice file output from the WPP code
%
%      function [u] = wppgetimage( fil, plt )
%
%            Input: fil - Name of file where the output image is.
%                   plt - 0 --> no plotting, 1 --> draw simple color plot.
%            Output: u  - The image (data matrix)
%
  function [u] = wppgetimage( fil, plt )

  if nargin < 2 
     plt = 0;
  end;

  fd = fopen(fil,'r');
  if fd < 0 
     disp( ['Error: could not open file ' fil] );
  else
     dh= fread(fd, 1, 'float');
     ni= fread(fd,1,'int');
     nj= fread(fd,1,'int');
     [u,nr] = fread(fd,[ni nj],'float');
     if nr ~= ni*nj 
        disp( ['Warning, could not read all of the matrix']);
     end;
     fclose(fd);

     if plt == 1 
        clf;
        set(gca,'FontSize',15);
        surf(u);
        view(2);
     end;
  end;
