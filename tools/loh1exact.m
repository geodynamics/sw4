% function loh1exact
%
% usage: [t,ra,tr,ve] = loh1exact(sigma)
%
% Input:
%       sigma: (Optional argument, default value sigma=0.06) 
%              Spread in the Gaussian source-time function 
function [t,ra,tr,ve] = loh1exact(sigma)

if nargin < 1
  sigma = 0.06;
end

sig=sigma
% Filename for exact solution
filename='LOH.1_prose3';
ReadUHS

indx=(1:2000);

t = t(indx);
ra = ra(indx);
tr = tr(indx);
ve = ve(indx);


