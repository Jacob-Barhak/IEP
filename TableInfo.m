function [ Data, Dim, Sizes, Ranges] = TableInfo( InTable )
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
% [ Data, Dim, Sizes, Ranges] = TableInfo( InTable )
% The function parses the InTable information and returns the following:
% Data    : Multidimensional array that contains only the table data
% Dim     : A vector containing the list of dimension names
% DimSize : A vector corresponding to DIM with the length of each dimension
%           This should be the same as size(Data)
% Ranges  : A cell array that its cells contain the range bounds that split 
%           the table cell ranges in case of a dimension over a range. 
%           These numbers indicate indices to access the matrix in case of 
%           an indexed matrix, in which case the vector starts with NaN. 
%           Note that the Ranges matrix is of size length(Dim),max(Size)+1
 
 
% If InTable is not a structure then this is not a table
if ~isstruct(InTable),
    error('Input variable is not a table')
end
 
% extract the dimensions from the string
Dim = symvar(char(InTable.Dim));
 
Data = InTable.Data;
 
% Get the size information from the data
Sizes=size(Data);
 
% Take care of vectors that matlab shows as matrices
if Sizes(2)==1 && length(Sizes)==2,
    Sizes=Sizes(1);
end 
 
 
% Structure the dimension and the range data from the given vector by
% advancing an index and reading the data in that position in the vector.
 
% Check the validity of the sizes
NumOfExpecedDimParam = 2*length(Sizes)+sum(Sizes);
NumberOfUnityDimensions = (length(InTable.Dim) - NumOfExpecedDimParam) / 3;
 
 
if round(NumberOfUnityDimensions) ~= NumberOfUnityDimensions,
    % This means that the indices do not conform to the sizes
    error('Invalid dimension string - does not conform to data sizes');
end
 
 
% extend the sizes vector
if NumberOfUnityDimensions>0,
    Sizes(end+(1:NumberOfUnityDimensions))=1;
end
 
% Initialize Ranges
Ranges = {};
 
% check for consistency
if length(Dim)~=length(Sizes),
    error ('Number of Dimensions defined by data and by ranges is not matching');
end
 
i=1;
for d = 1:length(Sizes),
    if Dim{d} ~= char(InTable.Dim(i)),
        error('Dimension %s defined by data reads as %s in the range vector', Dim{d},char(InTable.Dim(i)));
   end
   Ranges{d} = double(InTable.Dim(i+1:i+Sizes(d)+1));
   i=i+Sizes(d)+2;
end
 


