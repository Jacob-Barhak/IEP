function [ SubVec ] = MyInd2Sub( TableSize, Index )
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
% function [ SubVec ] = MyInd2Sub( TableSize, Index )
% The function is similar to Matlab's ind2sub function and uses it. 
% The difference is that the output is a single vector, rather than a
% separated list of variables.
% The function translates an Index of a vector that signifies the spread 
% multidimensional table into a subscript index vector SubVec corresponding 
% to the table size provided in Table size.
 
% Just some consistency checks
if Index>prod(TableSize),
    error('Subscript is out of index range for the given table');
end
 
% Create a variable vector for each dimension in TableSize
Str='';
for d = 1:length(TableSize),
    Str = [ Str 'i' num2str(d) ','] ;
end
% Delete the command at the end of the string
Str = Str(1:(end-1));
% Construct a vector
StrVec = ['[' Str ']'];
 
% Construct a call sub2ind and execute it
CommandStr = [ StrVec '=ind2sub(TableSize,Index);'];
eval (CommandStr);
 
% Convert the output to a vector
SubVec = eval (StrVec);
