%
% CONVERTJDAY
%
%    Convert a Julian day (1 to 366) into a day of a month.
%
%             [day,mon]=convertjday( jday, year )
%
%                Input: jday - Julian day
%                       year - Year of 'jday'
%                Output: day - Day in month of jday
%                        mon - Month of jday.
%
function [day,mon]=convertjday( jday, year )

if jday > 0 & jday < 367
   day = 1;
   jd  = 1;
   mon = 1;
   lastofmonth=[31 28 31 30 31 30 31 31 30 31 30 31];
% leapyear correction
   if ( mod(year,400) == 0 ) || ( (mod(year,4) == 0) && ~(mod(year,100) == 0) )
      lastofmonth(2) = 29;
   end;
   while jd < jday 
      jd = jd+1;
      day = day + 1;
      if day > lastofmonth(mon)
         day = 1;
         mon = mon + 1;
      end;
   end;
end;
