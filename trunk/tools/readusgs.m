%
% READUSGS
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios.
%
%              [t,ux,uy,uz,uxy,uxz,uyz]=readusgs(filename)
%
%       Input: filename - Name of receiever data file
%          
%       Output: t  - Time vector of the same length as the data vectors (1st column of data).
%               ux - X- or East-West direction data component (2nd column of data )
%               uy - Y- or North-South direction data component (3rd column of data )
%               uz - Z- or Up direction data component (4th column of data)
%               uxy, uxz, uyz - When strains are output (ux,uy,uz) are the diagonal
%                    components and these are the off diagonals.
%
%       Note: When the divergence is output, ux is the divergence and there is only one component.
%       Note: [ux,uy,uz] is [East-West,North-South,Up] or [X,Y,Z] components
%       depending on how the data file was written by WPP. 
%
function [t,u1,u2,u3,u4,u5,u6] = readusgs( fname )

fd=fopen(fname,'r');
for i=1:9
   lin = fgetl(fd);
end;
nc = str2num(lin(12:end));
for i=1:nc
   lin = fgetl(fd);
end;
q=fscanf(fd,'%f');
t=q(1:nc:end,1);
for c=2:nc
   eval(['u' int2str(c-1) '=q(' int2str(c) ':nc:end,1);']);
end;
fclose(fd);
