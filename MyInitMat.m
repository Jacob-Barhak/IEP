function Mat = MyInitMat (Vec,Value)
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
%function Mat = MyInitMat (Vec,Value)
% This function replaces Matlab zeros and ones. It behaves in a similar 
% manner, except that it always assumes the input is a vector and it 
% creates a vector if the length of Vec is 1 rather than a matrix of size 
% Vec x Vec. All the output values of the matrix will be equal to Value
% No special consistency checks are made on the input variables
 
if nargin <1,
    % if no input arguments, just return zero
    Mat = 0;
    return
elseif nargin<2,
    % If value not specified, then assume to init to zero
    Value = 0;
end
 
if length(Value)~=1,
    error('Cannot init to a non scalar defined in the "Value"')
end
 
if length(Vec)==1,
    TempMat = ones([Vec 1]);
else
    TempMat = ones(Vec);
end    
 
Mat = TempMat*Value;
