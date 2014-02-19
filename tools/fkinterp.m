%
%  FKINTERP
%
%  USAGE:
%        [radi,trani,upi]=fkinterp( t, rad, tran, up, tw, style, doplot, hld )
%
%  ARGUMENTS:
%        Input: NOTE: the vectors t, rad, tran, and up must all have the same number of elements
%            t:       Vector of time values for the following three components
%            rad:     Vector of 1st component values (radial)
%            tran:    Vector of 2nd component values (transverse)
%            up:      Vector of 3rd component values (vertical)
%            tw:      Vector of monotonically increasing new time levels
%            style:   (Optional) Plot style, e.g., 'b-' for solid blue lines. Defaults to 'k-'.
%            doplot:  (Optional) 1 for plotting.                              Defaults to 0.
%            hld:     (Optional) 1 for doing 'hold on' before plotting.       Defaults to 0.
%        Output:
%            radi:    Vector of interpolated values for 1st component
%            trani:   Vector of interpolated values for 2nd component
%            upi:     Vector of interpolated values for 3rd component
%
function [radi,trani,upi]=fkinterp( t, rad, tran, up, tw, style, doplot, hld )
if nargin < 8
  hld = 0;
end
if nargin < 7
  doplot = 0;
end
if nargin < 6
  style = 'k-';
end;

% radial-component
radi =interp1(t,rad,tw,'cubic',0.0);

% transverse-component
trani =interp1(t,tran,tw,'cubic',0.0);

% vertical-component
upi =interp1(t,up,tw,'cubic',0.0);

if doplot == 1
   subplot(3,1,1);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(tw,radi,style);
   axis([tw(1) tw(end) min(radi) max(radi)]);

   subplot(3,1,2);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(tw,trani,style);
   axis([tw(1) tw(end) min(trani) max(trani)]);

   subplot(3,1,3);
   if hld==1
     hold on;
   else
     hold off;
   end;
   plot(tw,upi,style);
   axis([tw(1) tw(end) min(upi) max(upi)]);
end;

