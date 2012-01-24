%
% IMAGESIZE
%
%     imagesize( fil )
%
%   Outputs basic grid size and number of grid points information about an image file 'fil'
%
function imagesize( fil)
fd=fopen(fil,'r');
pr=fread(fd,1,'int');
nb=fread(fd,1,'int');

for i=1:nb
      h(i)  = fread(fd,1,'double');
      ib(i) = fread(fd,1,'int');
      ie(i) = fread(fd,1,'int');
      jb(i) = fread(fd,1,'int');
      je(i) = fread(fd,1,'int');
      printf("i=%i, h=%e, ib=%i, ie=%i, jb=%i, je=%i\n", i, h(i), ib(i), ie(i), jb(i), je(i));
end;
fclose(fd);
