%
% READFILTPLOT
%
%  Read seismogram data in format specified by USGS (no header data)
%
%              [tf,ef,nf,uf]=readfiltplot( filename, freq, dtf, colorstring, erasefirst )
%
%       Input: filename - Name of receiever data file (no header info)
%              freq - Corner frequency in Butterworth low-pass filtering
%              dtf - time step for low-passed time sequence
%              colorstring - plot style, for example 'r' to make red solid lines
%              erasefirst - 1 erases the current figure, 0 does a 'hold on' before plotting
%              
%          
%       Output: tf - Time vector with step size 'dtf', of the same length as the data vectors
%               ef - East-West data component 
%               nf - North-South data component
%               uf - Vertical (positive up) data component
%               
%
function [tf,ef,nf,uf] = readfiltplot( filename, freq, dtf, colorstring, erasefirst )

% read the raw data
fd=fopen(filename,'r');
q=fscanf(fd,'%f');
t=q(1:4:end,1);
e=q(2:4:end,1);
n=q(3:4:end,1);
u=q(4:4:end,1);

%filter
ef = usgsfilter(t(2)-t(1),e,freq,dtf);
nf = usgsfilter(t(2)-t(1),n,freq,dtf);
uf = usgsfilter(t(2)-t(1),u,freq,dtf);

% time sequence
N=length(ef);
tf=t(1) + (0:N-1)*dtf;

%plot
plotcomp(tf,ef,nf,uf,colorstring,erasefirst);
