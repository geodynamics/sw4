%
% READNORM
% Syntax:
% [t,u] = readnorm( filename )
%
%  Read receiever data in format specified by USGS and output the norm of the 3 solution components
function [t,u] = readnorm( filename )

[t ux uy uz]=readusgs(filename);
b1 = 1;
stride1=1;
n1=length(ux);

u = sqrt((ux(b1:stride1:n1)).^2 + (uy(b1:stride1:n1)).^2 + (uz(b1:stride1:n1)).^2);
