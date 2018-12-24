%
% PLOT3FILES
%
% Usage:  plot3files( basedir, basefile )
%
%       Input: basedir and basefile
%               
function plot3files( basedir, basefile )

  dir0 = sprintf("%s_output0", basedir)
  dir4 = sprintf("%s_output4", basedir)
  dird = sprintf("%s_data", basedir)
  fileout = sprintf("%s_out.txt", basefile)
  filedata = sprintf("%s.txt", basefile)

  file1=sprintf("%s/%s", dir0, fileout)
  file2=sprintf("%s/%s", dir4, fileout)
  file3=sprintf("%s/%s", dird, filedata)

  plotusgs(file1,"b",1);
  plotusgs(file2,"r",0);
  plotusgs(file3,"k",0);
end
