function [ OutProb ] = MyCDF( Name, x, varargin )
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
% [ OutProb ] = MyCDF( Name, x, varargin )
% This function calculates the CDF for a distribution given in Name. It is
% based on the Matlab cdf function and encapsulates it while providing
% support specific to the prototype. Here are several differences:
% 1. Bernoulli distribution is recognized as 'ber' and defined as a 
%    sub case of Binomial distribution. The format for
%    this case is: MyCdf('Ber',0.3,0.7)
% 2. Geometric distribution is defined as pq^(x-1) rather than pq^x
% 3. Multinomial distribution is partially supported to represent a single
%    test. It is defined by a vector of probabilities for each possible
%    outcome. Multinomial distribution with several tests is not supported
%    by the prototype and therefore its CDF is not required. The format for
%    this case for outcome 2 for one test is: 
%    MyCdf('Mulnom',2, 1, [ 0.1 0.2 0.3 0.4])
% 4. If x is NaN the function returns zero to support the convention used
%    in table definition in the prototype.
 
% If x is NaN return zero
if isnan(x)
    OutProb = 0;
    return
end
 
% Default arguments to pass to cdf
xx = x;
NameArg = Name;
Arguments = varargin;
 
% In case of a geometric distribution, change x to x-1. Later the Matlab
% cdf will be used with this modified value
if strcmpi(Name,'Geo') ||  strcmpi(Name,'Geometric'),
    xx = x - 1;
end
 
% In case of a Bernoulli distribution, change it to a Binomial with n=1
if strcmpi(Name,'Ber') ||  strcmpi(Name,'Bernoulli'),
    NameArg = 'bino';
    Arguments = [1,varargin];
end
  
% If this point in the program was reached, just call Matlab cdf with the
% appropriate number of variables
switch length(Arguments),
    case 0,
        OutProb = cdf(NameArg, xx);
    case 1,
        OutProb = cdf(NameArg, xx, Arguments{1});
    case 2,
        OutProb = cdf(NameArg, xx, Arguments{1}, Arguments{2});
    otherwise,
        error('Invalid number of parameters. The prototype is not designed to support distributions with 3 or more input arguments. Please make sure the input data is correct.');
end
