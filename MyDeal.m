function [FuncOut, GradOut, HessOut] =  MyDeal (FuncIn, GradIn, HessIn)
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
% function [FuncOut, GradOut, HessOut] =  MyDeal (FuncIn, GradIn, HessIn)
% Given three input, the function returns 1,2, or 3 outputs depending on the
% number of variables requested at output. It can be used to encapsulate a
% function handle to regulate outputs.

if nargout >= 1,
    FuncOut = FuncIn;
end
if nargout >= 2,
    GradOut = GradIn;
end
if nargout >= 3,
    HessOut = HessIn;
end
if nargout > 3 || nargout < 1,
    error('cannot have more than 3 or less than 1 input arguments')
end
