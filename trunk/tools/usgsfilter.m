%
% USGSFILTER
%
%  Filter seismic receiver data as done by USGS for the Hayward
%  fault earthquake scenarios.
%
%      [uf]=usgsfilter( dt, u, freq, dtf )
%
%           Input: dt   - Time step of data.
%                  u    - Data to be filtered, uniformly
%                         distributed with spacing dt.
%                  freq - (Optional) Filter to this cutoff frequency (Hz)
%                           Default is 0.25 Hz.
%                  dtf  - (Optional) Interpolate result to this spacing.
%                           Default is 0.05 s.
%           Output: uf - Filtered u, given with spacing dtf.
%       
function [uf]=usgsfilter( dt, u, freq, dtf )

if nargin < 4
  dtf = 0.05;
end;
if nargin < 3
  freq = 0.25;
end;
n = length(u);
tfinal = (n-1)*dt;
[b,a]=butter(2,freq*dt*2);
uif=filtfilt(b,a,u);
o = round(tfinal/dtf)+1;
uf=interp1((0:n-1)*dt,uif,(0:o-1)*dtf);

