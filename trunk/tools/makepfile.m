% MAKEPFILE
%
%
function makepfile()

wc1="absoluteseis_ALL_DATA_P*";
pnames=ls(wc1);
npfiles=size(pnames,1);

wc2="absoluteseis_ALL_DATA_S*";
snames=ls(wc2);
nsfiles=size(snames,1);

nlat=21;
nlon=21;
ndepth=20;

depth = [0 0.01 0.025 0.05 0.075 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 2.5 3 3.5];

p=zeros(nlon,nlat,ndepth);
s=zeros(nlon,nlat,ndepth);

%for k=1:1

% read all Vp values
for k=1:npfiles
% deblank removes trailing white space
  fn1 = deblank(pnames(k,:));
  rem = rindex(fn1,"km.xyz");
  start = length("absoluteseis_ALL_DATA_P_");
  nstr = fn1(start+1:rem-1);
  dp = str2num(nstr);
  idx = lookup(depth, dp);
  printf("Filename #%i, '%s', nstr='%s', dp=%f, depth(%i)=%f,  \n", k, fn1, nstr, dp, idx, depth(idx))

% open the file and read the Vp data
  fp = fopen(fn1,"r");

  for j=1:nlat
    for i=1:nlon
      [val cnt]=fscanf(fp,"%e %e %e", 3);
      lon(i) = val(1);
      lat(j) = val(2);
      p(i,j,idx) = val(3);
    endfor;
  endfor;

%  p(:,:,idx)

  fclose(fp);
endfor;

% read all Vs values
for k=1:nsfiles
% deblank removes trailing white space
  fn1 = deblank(snames(k,:));
  rem = rindex(fn1,"km.xyz");
  start = length("absoluteseis_ALL_DATA_S_");
  nstr = fn1(start+1:rem-1);
  dp = str2num(nstr);
  idx = lookup(depth, dp);
  printf("Filename #%i, '%s', nstr='%s', dp=%f, depth(%i)=%f,  \n", k, fn1, nstr, dp, idx, depth(idx))

% open the file and read the Vp data
  fp = fopen(fn1,"r");

  for j=1:nlat
    for i=1:nlon
      [val cnt]=fscanf(fp,"%e %e %e", 3);
      s(i,j,idx) = val(3);
    endfor;
  endfor;

%  s(:,:,idx)

  fclose(fp);
endfor;

% tmp: look at one plane
  p(:,:,1)

% write the pfile
dens=2.5;
delta = 0.004; % average step size

fp=fopen("eric.ppmod","w");

% header data
fprintf(fp,"Newberry cross-correlation\n");
fprintf(fp,"%e\n", delta);
fprintf(fp,"%i %e %e\n", nlat, lat(1), lat(nlat));
fprintf(fp,"%i %e %e\n", nlon, lon(1), lon(nlon));
fprintf(fp,"%i %e %e\n", ndepth, depth(1), depth(ndepth));
fprintf(fp,"-99 -99 -99 -99\n");
fprintf(fp,"F\n");

% data
for j=1:nlat
  for i=1:nlon
% first line of depth profile
    fprintf(fp, "%e %e %i\n", lat(j), lon(i), ndepth);
% material properties
    for k=1:ndepth
      fprintf(fp,"%i %e %e %e %e\n", k, depth(k), p(i,j,k), s(i,j,k), dens);
    endfor;
  endfor;
endfor;

fclose(fp);

end
