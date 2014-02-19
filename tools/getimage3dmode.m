%
% GETIMAGE3DMODE
%
%   Translate an integer image3d mode number to the corresponding mode string.
%
%   NOTE: These mode numbers are different from the 2D image mode numbers
%
%           str=getimage3dmode(mode)
%                Input: mode - Integer (currently 0 to 15)
%                Output: str - Name of mode.
function str=getimage3dmode(mode)
if mode == 0
   str='none';
elseif mode == 1
   str='ux';
elseif mode == 2
   str = 'uy';
elseif mode == 3
   str = 'uz';
elseif mode == 4
   str = 'rho';
elseif mode == 5
   str = 'lambda';
elseif mode == 6
   str = 'mu';
elseif mode == 7
   str = 'p';
elseif mode == 8
   str = 's';
elseif mode == 9
   str = 'gradRho';
elseif mode == 10
   str = 'gradMu';
elseif mode == 11
   str = 'gradLambda';
elseif mode == 12
   str = 'gradP';
elseif mode == 13
   str = 'gradS';
elseif mode == 14
   str = 'qp';
elseif mode == 15
   str = 'qs';
else
  str = 'unknown';
end;
