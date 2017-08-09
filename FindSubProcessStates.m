function [PooledStates, IsSubProcess, NestingLevel] = FindSubProcessStates( AllStates, PoolingSubProcess, StatesFilter, CurrentNestingLevel)
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
% function [PooledStates, IsSubProcess, NestingLevel] = FindSubProcessStates( AllStates, PoolingSubProcess, StatesFilter, CurrentNestingLevel)
% The function drills down a sub-process provided in PoolingSubProcess and 
% returns all the states from AllStates associated with this sub-process in 
% PooledStates. 
% Note that the PoolingSubProcess itself is considered the pooling \
% sub-process and appears at the output. 
% NestingLevel is an output vector corresponding to pooled states that 
% holds the nesting level where the pooling sub process is level 0, its 
% children 1 grandchildren 2 etc.
% IsSubProcess is an output vector corresponding to PooledStates that has 1
% if the output state is a sub-process. Note that it is not logical although
% it is composed of 1's and 0's.
% StateFilter is a vector of states to be considered that is 
% used in the recursion to prevent infinite loops. If it is not defined, 
% then all states in AllStates are considered. CurrentNestingLevel is also
% used in the recursion to keep track of recursion level
 
if nargin < 3,
    % Meaning that States should be defined as all states 
    StatesFilter = 1:length(AllStates);
    CurrentNestingLevel=0;
end
 
if isempty(AllStates(PoolingSubProcess).Name),
    error('State %i is used in a subprocess or as an input variable - yet it is undefined. Provide a name to this state',PoolingSubProcess)
end
 
% First add the state to the output vector
PooledStates = PoolingSubProcess;
NestingLevel = CurrentNestingLevel;
 
if ~isempty(AllStates(PoolingSubProcess).PoolingData),
 
    % Possible children of a state are all the states in the pooling data
    % matrix where the prevalence is zero. If not zero, this means the 
    % state is a pooling state for a study rather than a model process.
    PossibleChildren = AllStates(PoolingSubProcess).PoolingData(1, ~AllStates(PoolingSubProcess).PoolingData(2,:));
 
    IsSubProcess = ~isempty (PossibleChildren);
    
    % Raise an error if the state is both a pooling state in a study and a
    % sub-process in a model
    if IsSubProcess && any(AllStates(PoolingSubProcess).PoolingData(2,:)),
        error ('State number %i is defined as a pooling state in a study and as a model subprocess',PoolingSubProcess);
    end
else
    IsSubProcess = 0;
    PossibleChildren=[];
end
 
% If the state is a pooling state, then loop through all the states in
% the sub-process and recursively call the function for each such state
for ChildState = PossibleChildren,
    
    % Remove output states from the remaining states list
    NewFilter = setdiff(StatesFilter, PooledStates);
    [ PooledStatesReturned, IsSubProcessReturned, NestingLevelReturned ] = FindSubProcessStates( AllStates, ChildState, NewFilter, CurrentNestingLevel+1);
    
    % Check if a state has already been used in another sub process
    % intersect the current resulting states with the previously
    % accumulated list. If a state was defined twice, the vector will not
    % be empty
    DoubleDefined = intersect(PooledStatesReturned,PooledStates);
    % Check for state definition consistency by checking that 
    if ~isempty(DoubleDefined)
        error('States [%s] are defined in more than one sub-process', num2str(DoubleDefined));
    end
       
    % Add the output to the results
    PooledStates = [PooledStates PooledStatesReturned];
    NestingLevel = [NestingLevel NestingLevelReturned];
    IsSubProcess = [IsSubProcess IsSubProcessReturned];
end
 
