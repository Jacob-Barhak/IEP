function [ f ] = MyLog( x )
% Copyright (C) 2007,2009 The Regents of the University of Michigan
% Initially developed by Deanna Isaman, Jacob Barhak
%
% This file is part of the Indirect Estimation Prototype. The Indirect
% Estimation Prototype is free software: you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% The Indirect Estimation Prototype is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [ f ] = MyLog( x )
% This function can replace the log calculation for optimization purposes. 
% It calculates the natural logarithm for positive numbers larger than
% eps(0), which is the smallest number that can be represented. For numbers
% smaller than eps(0) it returns a linear function connected to that point 
% with c0 continuity. This allows calculating log with zero or 
% negative numbers, which allows log likelihood functions to be optimized.
% Mathematically, Mylog has this definition: 
% Mylog(x) = log(x)                                 where x>=eps(0)
%            log(realmin)                           where x>0 & x<eps(0)
%            log(eps(0))-eps(0)-x = log(eps(0))+x   where x<=0
% Note that near zero or negative numbers are no longer meaningful, nor may
% be numerically calculated correctly, and mylog(0) will not be -inf, yet 
% it avoids returning an error during optimization and will report declining
% numbers when x declines.
 
f=[];
for xx = x
    if xx>=eps(0),
        ff = log(xx);
    elseif xx<eps(0) && xx>0,
        ff = log(realmin);
        warning('MyLog:SmallNum','MyLog handles a positive number smaller than realmin');
    elseif xx<=0
        %disp('Info message: MyLog handles a negative number');
        ff = log(eps(0))-eps(0)-xx;
    elseif isnan(x),
        warning('MyLog:NaN','MyLog returns NaN for NaN input');
        %error('MyLog returns NaN for NaN input');
        ff = NaN;
    else    
        error('Assertion error in MyLog, x is in a range which is unhandled by the code');        
    end
    f(end+1) = ff;
end    
