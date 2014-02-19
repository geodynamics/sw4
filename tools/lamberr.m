%
% LAMBERR
%
% Evaluate error in numerical solution by comparing to exact solution
%
% [merr msol]=lamberr(fname, dist, lambpath, verb)
%
% required arguments:
% dist:  must be positive
% fname: the name of the solution file
% lambpath: path to lamb1 executable
%
% optional argument:
% verb = 1 gives verbose output, verb=0 by default
%
function [merr msol]=lamberr(fname, dist, lambpath, verb)
if nargin < 4
  verb = 0;
end
% 
% read the numerical solution (assume path has been set to find readusgs.m
[t, ux, uy, uz]=readusgs(fname);

% number of entries
nt = length(t);

% number of time steps is one less
nt = nt - 1;
% end time
tmax = max(t);
% lamb1 path
%path="../../../optimize_v1.0";
% calling string
cmd = sprintf("%s/lamb1 -nsteps %i -dist %e -tmax %e", lambpath, nt, dist, tmax);

if (verb == 1)
   printf("cmd:'%s'\n", cmd);
end

% call lamb1 to generate exact solution
[stat output]=system(cmd);

if (verb == 1)
   printf("return status = %i\ncommand output:\n%s\n", stat, output);
end

% exact solution is in file uzex.dat
uzex = load("uzex.dat");
merr = max(abs(uzex(:,2)-uz));
msol = max(abs(uzex(:,2)));
