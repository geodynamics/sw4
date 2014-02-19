%
% ELAMBERR
%
% Evaluate error in several numerical solutions
%
% [merr]=elamberr(bdir, bname, lambpath, verb)
%
% required arguments:
% bdir:  base directory name, e.g. "lamb-sw4-nx"
% bname: base file name, e.g. "sw4-" or "wpp-"
% lambpath: path to lamb1 executable
%
% optional argument:
% verb = 1 gives verbose output, verb=0 by default
%
function [emat]=elamberr(bdir, bname, lambpath, verb)
if nargin < 4
  verb = 0;
end

% default resolutions
NN = [10 15 20 30 40 60 80];
% number of grid resolutions (7 for WPP)
nres=7;

% default distances
DD = [1 2 3 4 5 6 7 8 9 10];
ndist=10;

emat=zeros(nres,ndist);

for k=1:nres
  dirname=sprintf("%s%i", bdir, NN(k))
  for q=1:ndist
    fname=sprintf("%s/%s%i.txt", dirname, bname, DD(q));
    [merr msol]=lamberr(fname, DD(q), lambpath);
    if verb == 1
      printf("res=%i, dist=%i, filename=%s, max-err=%e\n", NN(k), DD(q), fname, merr);
    end
    emat(k,q) = merr;
  end
end