%
% relmax
%
%  Read receiever data in format specified by USGS for the Hayward
%  fault earthquake scenarios.Then compute the norm of the time series for t>t1
%
%              [largest_rel_nrm] = relmax(directory)
%
%       Input: filename - Name of receiever data file
%              t1: [optional], default t1=0. Otherwise, start time for evaluating the norm
%              doplot: [optional] default=0. If doplot=1, plot the norm as fcn of time
%              linlog: [optional] default=0 for linear scale. If linlog=1, use a 
%                      logarithmic scale
%               
function [largest_rel_nrm] = relmax( directory )
t1 = 4.5;
wildcard=sprintf("%s/*.txt", directory);
filenames = dir(wildcard);
len = length(filenames);

rel_nrm = zeros(len,1);
fnames = zeros(len,1);
for k=1:len
  fname = sprintf("%s/%s", directory, filenames(k).name);
  rel_nrm(k) = normusgs( fname, t1);
end

[largest_rel_nrm im]=max(rel_nrm);

printf("Largest relative norm: %e, occured in the file: %s\n", largest_rel_nrm, filenames(im).name);