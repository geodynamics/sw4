%
% GETIMAGEMODESTR
%
%   Translate an integer image mode number to the corresponding mode string.
%
%           str=getimagemodestr(mode)
%                Input: mode - Integer (currently 0 to 32)
%                Output: str - Name of mode.
function str=getimagemodestr(mode)
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
   str = 'uxexact';
elseif mode == 10
   str = 'uyexact';
elseif mode == 11
   str = 'uzexact';
elseif mode == 12
   str = 'div';
elseif mode == 13
   str = 'curlmag';
elseif mode == 14
   str = 'divdt';
elseif mode == 15
   str = 'curlmagdt';
elseif mode == 16
   str = 'lat';
elseif mode == 17
   str = 'lon';
elseif mode == 18
   str = 'topo';
elseif mode == 19
   str = 'gridx';
elseif mode == 20
   str = 'gridy';
elseif mode == 21
   str = 'gridz';
elseif mode == 22
   str = 'uxerr';
elseif mode == 23
   str = 'uyerr';
elseif mode == 24
   str = 'uzerr';
elseif mode == 25
   str = 'magdudt';
elseif mode == 26
   str = 'hmagdudt';
elseif mode == 27
   str = 'hmaxdudt';
elseif mode == 28
   str = 'vmaxdudt';
elseif mode == 29
   str = 'mag';
elseif mode == 30
   str = 'hmag';
elseif mode == 31
   str = 'hmax';
elseif mode == 32
   str = 'vmax';
elseif mode == 33
   str = 'gradrho';
elseif mode == 34
   str = 'gradmu';
elseif mode == 35
   str = 'gradlambda';
elseif mode == 36
   str = 'gradp';
elseif mode == 37
   str = 'grads';
elseif mode == 38
   str = 'qp';
elseif mode == 39
   str = 'qs';
else
  str = 'unknown';
end;
