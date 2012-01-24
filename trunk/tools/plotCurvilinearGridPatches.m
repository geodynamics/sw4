%
% READIMAGEPATCH
%
%  Reads one grid patch on an image file
%
%      [im,x,y]=readimagepatch( fil, inr )
%
%        Input: fil  - Name of image file.
%               filx - grid fileX
%               fily - grid fileY
%
function err=plotCurvilinearGridPatches( fil,filx,fily)

fd=fopen(fil,'r');
pr=fread(fd,1,'int');
ni=fread(fd,1,'int');

clear figure;
clear im;
clear zmin;

[im1,x,y]  = readimagepatch(fil,1);
[x1,x,y]   = readimagepatch(filx,1);
[y1,x,y]   = readimagepatch(fily,1);
maxT       = max(max(im1));
minT       = min(min(im1));  
maxX       = max(max(x1));       
minX       = min(min(x1));       
maxY       = max(max(y1));       
minY       = min(min(y1));         

if (ni >= 2)
   [im2,x,y]  = readimagepatch(fil,2);
   [x2,x,y]   = readimagepatch(filx,2);
   [y2,x,y]   = readimagepatch(fily,2);
   maxT       = max(max(max(im2)),maxT);
   minT       = min(min(min(im2)),minT);     
   maxX       = max(max(max(x2)),maxX);        
   minX       = min(min(min(x2)),minX);       
   maxY       = max(max(max(y2)),maxY);       
   minY       = min(min(min(y2)),minY);           
end
if (ni >= 3)
   [im3,x,y]  = readimagepatch(fil,3);
   [x3,x,y]   = readimagepatch(filx,3);
   [y3,x,y]   = readimagepatch(fily,3);
   maxT       = max(max(max(im3)),maxT);
   minT       = min(min(min(im3)),minT);     
   maxX       = max(max(max(x3)),maxX);        
   minX       = min(min(min(x3)),minX);       
   maxY       = max(max(max(y3)),maxY);       
   minY       = min(min(min(y3)),minY);             
end
if (ni >= 4)
   [im4,x,y]  = readimagepatch(fil,4);
   [x4,x,y]   = readimagepatch(filx,4);
   [y4,x,y]   = readimagepatch(fily,4);
   maxT       = max(max(max(im4)),maxT);
   minT       = min(min(min(im4)),minT); 
   maxX       = max(max(max(x4)),maxX);        
   minX       = min(min(min(x4)),minX);       
   maxY       = max(max(max(y4)),maxY);       
   minY       = min(min(min(y4)),minY);                 
end
if (ni >= 5)
   [im5,x,y]  = readimagepatch(fil,5);
   [x5,x,y]   = readimagepatch(filx,5);
   [y5,x,y]   = readimagepatch(fily,5);
   maxT       = max(max(max(im5)),maxT);
   minT       = min(min(min(im5)),minT);  
   maxX       = max(max(max(x5)),maxX);        
   minX       = min(min(min(x5)),minX);       
   maxY       = max(max(max(y5)),maxY);       
   minY       = min(min(min(y5)),minY);         
end
if (ni >= 6)
   [im6,x,y]  = readimagepatch(fil,6);
   [x6,x,y]   = readimagepatch(filx,6);
   [y6,x,y]   = readimagepatch(fily,6);
   maxT       = max(max(max(im6)),maxT);
   minT       = min(min(min(im6)),minT);
   maxX       = max(max(max(x6)),maxX);        
   minX       = min(min(min(x6)),minX);       
   maxY       = max(max(max(y6)),maxY);       
   minY       = min(min(min(y6)),minY);           
end
if (ni >= 7)
   [im7,x,y]  = readimagepatch(fil,7);
   [x7,x,y]   = readimagepatch(filx,7);
   [y7,x,y]   = readimagepatch(fily,7);
   maxT       = max(max(max(im7)),maxT);
   minT       = min(min(min(im7)),minT); 
   maxX       = max(max(max(x7)),maxX);        
   minX       = min(min(min(x7)),minX);       
   maxY       = max(max(max(y7)),maxY);       
   minY       = min(min(min(y7)),minY);          
end
if (ni >= 8)
   [im8,x,y]  = readimagepatch(fil,8);
   [x8,x,y]   = readimagepatch(filx,8);
   [y8,x,y]   = readimagepatch(fily,8);
   maxT       = max(max(max(im8)),maxT);
   minT       = min(min(min(im8)),minT); 
   maxX       = max(max(max(x8)),maxX);        
   minX       = min(min(min(x8)),minX);       
   maxY       = max(max(max(y8)),maxY);       
   minY       = min(min(min(y8)),minY);          
end
if (ni >= 9)
   [im9,x,y]  = readimagepatch(fil,9);
   [x9,x,y]   = readimagepatch(filx,9);
   [y9,x,y]   = readimagepatch(fily,9);
   maxT       = max(max(max(im9)),maxT);
   minT       = min(min(min(im9)),minT); 
   maxX       = max(max(max(x9)),maxX);        
   minX       = min(min(min(x9)),minX);       
   maxY       = max(max(max(y9)),maxY);       
   minY       = min(min(min(y9)),minY);         
end
if (ni >= 10)
   [im10,x,y] = readimagepatch(fil,10);
   [x10,x,y]   = readimagepatch(filx,10);
   [y10,x,y]   = readimagepatch(fily,10);
   maxT       = max(max(max(im10)),maxT);
   minT       = min(min(min(im10)),minT); 
   maxX       = max(max(max(x10)),maxX);        
   minX       = min(min(min(x10)),minX);       
   maxY       = max(max(max(y10)),maxY);       
   minY       = min(min(min(y10)),minY);         
end

contourf(x1,-y1,im1,linspace(minT,maxT));
axis([minX,maxX,-maxY,-minY]);
hold on 
if (ni >= 2)
   contourf(x2,-y2,im2,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 3)
   contourf(x3,-y3,im3,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 4)
   contourf(x4,-y4,im4,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 5)
   contourf(x5,-y5,im5,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 6)
   contourf(x6,-y6,im6,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 7)
   contourf(x7,-y7,im7,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 8)
   contourf(x8,-y8,im8,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 9)
   contourf(x9,-y9,im9,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
if (ni >= 10)
   contourf(x10,-y10,im10,linspace(minT,maxT));
   axis([minX,maxX,-maxY,-minY]);
   hold on
end
colorbar

err = 0;
