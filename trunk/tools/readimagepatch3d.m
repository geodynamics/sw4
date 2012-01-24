function [im3d,h,time,x,y,z]=readimagepatch3d( file, pnr )
fd=fopen(file,'r');
prec=fread(fd,1,'int');
npatches=fread(fd,1,'int');
time = fread(fd,1,'double');
if pnr < 1 | pnr > npatches 
   disp(['Error in input patch nr']);
   disp(['File contains ' num2str(npatches) ' patches']);
else
   offset = 0;
   for p=1:pnr-1
     fread(fd,4,'double');
     dim=fread(fd,6,'int');
     nvals = (dim(6)-dim(5)+1)*(dim(4)-dim(3)+1)*(dim(2)-dim(1)+1);
     offset  = offset + nvals;
   end;
   h=fread(fd,1,'double');
   xmin=fread(fd,1,'double');
   ymin=fread(fd,1,'double');
   zmin=fread(fd,1,'double');
   dim=fread(fd,6,'int');
   fseek(fd,(npatches-pnr)*(6*4+4*8),'cof');
   if( prec == 4 )
     fseek(fd,offset*4,'cof');
   else
     fseek(fd,offset*8,'cof');
   end;
   npts = (dim(6)-dim(5)+1)*(dim(4)-dim(3)+1)*(dim(2)-dim(1)+1);
   if prec == 4
      im3 = fread(fd,npts,'float');
   else
      im3 = fread(fd,npts,'double');
   end;
   im3d = reshape(im3,dim(2)-dim(1)+1,dim(4)-dim(3)+1,dim(6)-dim(5)+1);
   z = ((dim(5):dim(6))-1)*h+zmin;
   y = ((dim(3):dim(4))-1)*h+ymin;
   x = ((dim(1):dim(2))-1)*h+xmin;
end;
