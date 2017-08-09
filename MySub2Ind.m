function [ Index ] = MySub2Ind( TableSize, SubVec )
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
% function [ Index ] = MySub2Ind( TableSize, SubVec )
% The function is similar to Matlab's sub2ind function and uses it. 
% The difference is that the input is a single vector, rather than a
% separated list of variables.
% The function translates a subscript index vector SubVec into an Index of
% a vector that signifies the flattened multidimensional table while the table
% size is provided in TableSize.
 
% Just some consistancy checks
if any(SubVec>TableSize | SubVec<=0),
    error('The index is out of range for the given table in at least one dimension');
end
 
% If only one diemnsion is defined in the table size, just check 
% consistancy and return the input as it is already a flat index to a 
% vector
if length(TableSize)<2,
    if TableSize < SubVec,
        error('Index exceeds vector size')
    end
    Index = SubVec;
    return
end
 
% Access the table in the transition with the given indices
 
% First convert the index vector to a string with comma
% separated values.
CommaSeparatedIndex = regexprep(int2str(SubVec),'\ *',',');
CommandStr = sprintf('sub2ind(TableSize,%s);',CommaSeparatedIndex);
Index=eval(CommandStr);
 
