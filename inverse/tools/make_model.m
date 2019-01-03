% make_model
% Usage: make_model(sf1, sf2, sf3, fname)
% Input:
%       sf1: 1st scale factor (rho)
%       sf2: 2nd scale factor (mu or vs)
%       sf3: 3rd scale factor (mbda or vp)
%       fname (optional): file name, default 'one-pr10.bin'
function make_model(sf1, sf2, sf3, fname)
  if nargin < 4
    fname = 'onep-pr10.bin';
  end
  ngrids = 3;

  fileID = fopen(fname,'w');
  fwrite(fileID,ngrids,'int');
  fwrite(fileID, sf1, 'double');
  fwrite(fileID, sf2, 'double');
  fwrite(fileID, sf3, 'double');
  fclose(fileID);
end
