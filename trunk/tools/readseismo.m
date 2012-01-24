%
% READSEISMO
%
%  Read seismogram data in format specified by USGS (no header data)
%
%              [t,ue,un,uu]=readseismo(filename)
%
%       Input: filename - Name of receiever data file
%          
%       Output: t  - Time vector of the same length as the data vectors (1st column of data).
%               ue - East-West direction data component (2nd column of data )
%               un - North-South direction data component (3rd column of data )
%               uu - Vertical (positive up) direction data component (4th column of data)
%               
%       Note: [ue,un,uu] is usually [East-West,North-South,Up] but could be other components depending on how the data was written
%
function [t,ue,un,uu] = readseismo( fname )

fd=fopen(fname,'r');
q=fscanf(fd,'%f');
t=q(1:4:end,1);
ue=q(2:4:end,1);
un=q(3:4:end,1);
uu=q(4:4:end,1);
