%
%  CONVSCALE
%
%   Convolve a time function and a solution operator 
%
%          [tc,w]=convscale( g, s, tg0, ts0, dt)
%
%        Input:  g      - Function g sampled at times tg0, tg0+dt, tg0+2*dt, ...
%                s      - Solution operator sampled at times ts0, ts0+dt, ts0+2*dt, ...
%                tg0    - start time for time series g
%                ts0    - start time for solution operator s
%                dt     - time step for time series and solution operator
%
%        Output: tc - Times where w is sampled.
%                w  - Solution: convolution between function g and solution operator s
%
   function [tc,w]=convscale( g, s, tg0, ts0, dt)
   w = conv(s,g)*dt;
   n = size(w)-1;
   tc = dt*(0:n) - tg0 + ts0;

