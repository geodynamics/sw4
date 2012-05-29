%
% READUSGS
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios.
%
%              [t,ux,uy,uz,uxy,uxz,uyz]=readusgs(filename, evtlon, evtlat, rotate)
%
%       Input: filename - Name of receiever data file
%              evtlon - event longitude [optional]
%              evtlon - event longitude [optional]
%              rotate - transform to radial, tangential, up (=-z) components [optional]
%          
%       Output: t  - Time vector of the same length as the data vectors (1st column of data).
%               ux - X- or East-West direction data component (2nd column of data )
%               uy - Y- or North-South direction data component (3rd column of data )
%               uz - Z- or Up direction data component (4th column of data)
%               uxy, uxz, uyz - When strains are output (ux,uy,uz) are the diagonal
%                    components and these are the off diagonals.
%       if (rotate)
%          call heading.m to calculate the azimuthal angle
%          rotate components to radial, tangential components (but leave vertical=z)
%       endif
%
%       Note: When the divergence is output, ux is the divergence and there is only one component.
%       Note: [ux,uy,uz] is [East-West,North-South,Up] or [X,Y,Z] components
%       depending on how the data file was written by WPP. 
%
function [t,u1,u2,u3,u4,u5,u6] = readusgs( fname, evtlon, evtlat, rotate )

if nargin < 4
  rotate = 0;
  evtlat = 0;
  evtlon = 0;
endif

fd=fopen(fname,'r');
for i=1:7
   lin = fgetl(fd);
end;
% read actual loction 
loc = sscanf(lin(53:end),"%e %e",2);
reclon= loc(1)
reclat= loc(2)

lin = fgetl(fd);
lin = fgetl(fd);
% read the number of components
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

% rotate?
if (rotate)
  [azdeg az] = heading(evtlon, evtlat, reclon, reclat);
  azdeg
  ur = u1.*cos(az) + u2.*sin(az);
  ut = u1.*sin(az) - u2.*cos(az);
  u1 = ur;
  u2 = ut;
endif