function [u]=readar(fname,c)
if nargin < 2
  c = 1;
end;
fd=fopen(fname,'r');
nc=fread(fd,1,'int');
ni=fread(fd,1,'int');
nj=fread(fd,1,'int');
nk=fread(fd,1,'int');
q1d=fread(fd,nc*ni*nj*nk,'double');
u=reshape(q1d(c:nc:end),ni,nj,nk);