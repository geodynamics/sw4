%
% READUSGSOLD
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios. Old format, before allowing any number of columns.
%
%              [t,ux,uy,uz]=readusgsold(filename)
%
%       Input: filename - Name of receiever data file
%          
%       Output: t  - Time vector of the same length as the data vectors (1st column of data).
%               ux - X- or East-West direction data component (2nd column of data )
%               uy - Y- or North-South direction data component (3rd column of data )
%               uz - Z- or Up direction data component (4th column of data)
%               
%       Note: [ux,uy,uz] is [East-West,North-South,Up] or [X,Y,Z] components
%       depending on how the data file was written by WPP. 
%
function [t,ux,uy,uz] = readusgsold( fname )

fd=fopen(fname,'r');
for i=1:12,fgetl(fd);end;
q=fscanf(fd,'%f');
t=q(1:4:end,1);
ux=q(2:4:end,1);
uy=q(3:4:end,1);
uz=q(4:4:end,1);
fclose(fd);
