%
% READIMAGEPATCH
%
%  Plot all grid patches on an image file
%
%      err=plotCurvilinearGridPatches( fil,zplot)
%
%        Input: fil  - Name of image file.
%          int zplot - 1 if xy plane, 0 otherwise.
function err=plotCurvilinearGridPatches( fil,zplot)

fd=fopen(fil,'r');
pr=fread(fd,1,'int');
ni=fread(fd,1,'int');

clf
clear im;
clear zmin;

[im1,x1,y1]  = readimagepatch(fil,1);
maxT         = max(max(im1));
minT         = min(min(im1));  
maxX         = max(max(x1));       
minX         = min(min(x1));       
maxY         = max(max(y1));
minY         = min(min(y1));       

if (ni >= 2)
   [im2,x2,y2] = readimagepatch(fil,2);
   if (zplot==0)
      y1 = y1+max(max(y2));
   end
   maxY         = max(max(y1));       
   minY         = min(min(y1));         

   maxT       = max(max(max(im2)),maxT);
   minT       = min(min(min(im2)),minT);     
   maxX       = max(max(max(x2)),maxX);        
   minX       = min(min(min(x2)),minX);       
   maxY       = max(max(max(y2)),maxY);      
   minY       = min(min(min(y2)),minY);          
end
if (ni >= 3)
   [im3,x3,y3] = readimagepatch(fil,3);
   if (zplot==0)
      y1 = y1+max(max(y3));
      y2 = y2+max(max(y3));
   end
   maxT       = max(max(max(im3)),maxT);
   minT       = min(min(min(im3)),minT);     
   maxX       = max(max(max(x3)),maxX);        
   minX       = min(min(min(x3)),minX);       
   maxY       = max(max(max(y3)),maxY);       
   minY       = min(min(min(y3)),minY);             
end
if (ni >= 4)
   [im4,x4,y4]  = readimagepatch(fil,4);
   if (zplot==0)
      y1 = y1+max(max(y4));
      y2 = y2+max(max(y4));
      y3 = y3+max(max(y4));
   end
   maxT       = max(max(max(im4)),maxT);
   minT       = min(min(min(im4)),minT); 
   maxX       = max(max(max(x4)),maxX);        
   minX       = min(min(min(x4)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);                 
end
if (ni >= 5)
   [im5,x5,y5]  = readimagepatch(fil,5);
   if (zplot==0)
      y1 = y1+max(max(y5));
      y2 = y2+max(max(y5));
      y3 = y3+max(max(y5));
      y4 = y4+max(max(y5));
   end
   maxT       = max(max(max(im5)),maxT);
   minT       = min(min(min(im5)),minT);  
   maxX       = max(max(max(x5)),maxX);        
   minX       = min(min(min(x5)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);         
end
if (ni >= 6)
   [im6,x6,y6]  = readimagepatch(fil,6);
   if (zplot==0)
      y1 = y1+max(max(y6));
      y2 = y2+max(max(y6));
      y3 = y3+max(max(y6));
      y4 = y4+max(max(y6));
      y5 = y5+max(max(y6));
   end
   maxT       = max(max(max(im6)),maxT);
   minT       = min(min(min(im6)),minT);
   maxX       = max(max(max(x6)),maxX);        
   minX       = min(min(min(x6)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);           
end
if (ni >= 7)
   [im7,x7,y7]  = readimagepatch(fil,7);
   if (zplot==0)
      y1 = y1+max(max(y7));
      y2 = y2+max(max(y7));
      y3 = y3+max(max(y7));
      y4 = y4+max(max(y7));
      y5 = y5+max(max(y7));
      y6 = y6+max(max(y7));
   end
   maxT       = max(max(max(im7)),maxT);
   minT       = min(min(min(im7)),minT); 
   maxX       = max(max(max(x7)),maxX);        
   minX       = min(min(min(x7)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);          
end
if (ni >= 8)
   [im8,x8,y8]  = readimagepatch(fil,8);
   if (zplot==0)
      y1 = y1+max(max(y8));
      y2 = y2+max(max(y8));
      y3 = y3+max(max(y8));
      y4 = y4+max(max(y8));
      y5 = y5+max(max(y8));
      y6 = y6+max(max(y8));
      y7 = y7+max(max(y8));
   end
   maxT       = max(max(max(im8)),maxT);
   minT       = min(min(min(im8)),minT); 
   maxX       = max(max(max(x8)),maxX);        
   minX       = min(min(min(x8)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);          
end
if (ni >= 9)
   [im9,x9,y9]  = readimagepatch(fil,9);
   if (zplot==0)
      y1 = y1+max(max(y9));
      y2 = y2+max(max(y9));
      y3 = y3+max(max(y9));
      y4 = y4+max(max(y9));
      y5 = y5+max(max(y9));
      y6 = y6+max(max(y9));
      y7 = y7+max(max(y9));
      y8 = y8+max(max(y9));
   end
   maxT       = max(max(max(im9)),maxT);
   minT       = min(min(min(im9)),minT); 
   maxX       = max(max(max(x9)),maxX);        
   minX       = min(min(min(x9)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);         
end
if (ni >= 10)
   [im10,x10,y10] = readimagepatch(fil,10);
   if (zplot==0)
      y1 = y1+max(max(y10));
      y2 = y2+max(max(y10));
      y3 = y3+max(max(y10));
      y4 = y4+max(max(y10));
      y5 = y5+max(max(y10));
      y6 = y6+max(max(y10));
      y7 = y7+max(max(y10));
      y8 = y8+max(max(y10));
      y9 = y9+max(max(y10));
   end
   maxT       = max(max(max(im10)),maxT);
   minT       = min(min(min(im10)),minT); 
   maxX       = max(max(max(x10)),maxX);        
   minX       = min(min(min(x10)),minX);       
   maxY       = max(max(max(y1)),maxY);       
   minY       = min(min(min(y1)),minY);         
end

linsp = linspace(minT,maxT,25);
contourf(x1,-y1,im1,linsp);
hold on 
if (ni >= 2)
   contourf(x2,-y2,im2,linsp);
end
if (ni >= 3)
   contourf(x3,-y3,im3,linsp);
end
if (ni >= 4)
   contourf(x4,-y4,im4,linsp);
end
if (ni >= 5)
   contourf(x5,-y5,im5,linsp);
end
if (ni >= 6)
   contourf(x6,-y6,im6,linsp);
end
if (ni >= 7)
   contourf(x7,-y7,im7,linsp);
end
if (ni >= 8)
   contourf(x8,-y8,im8,linsp);
end
if (ni >= 9)
   contourf(x9,-y9,im9,linsp);
end
if (ni >= 10)
   contourf(x10,-y10,im10,linsp);
end
axis([minX,maxX,-maxY,-minY]);
colorbar
hold off;

err = 0;
