function [ SubProbExpr , NameList] = ParseMarkovTerm( ProbExprIn, Index, Transitions)
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
% The function parses an expression passed to it in ProbExprIn. It 
% recognizes symbols that correspond to transition names provided in 
% Transitions. If a set of indices is provided in a vector Index, then 
% these indices are used to access a multidimensional table contained in the
% transition and replace the transition variable with the contents of the 
% appropriate table cell.
% The function returns the Transition names detected as a cell array in 
% NameList. It returns the substituted expression in SubProbExpr.
% Note that if Index is empty, no substitutions are made. This can help
% parse expressions other than members of Markov chain matrices.
% Note that ProbExprIn is either symbolic or a string. 
% The result returned is returned as a string
 
% Convert input variable to a string to treat symbolic objects
if isobject(ProbExprIn),
    ProbExpr = char(ProbExprIn);
else
    ProbExpr = ProbExprIn;
end
 
% First seperate all characters to an alphanumeric set and a delimiter set
AlphaNumeric = char([48:57,65:90,97:122,95]); % 0-9,A-Z,a-z,_
Delimiters = setdiff(1:255,AlphaNumeric);
 
 
% Init the reminder vector to the full string
Reminder = ProbExpr;
% Reset the output string to 
SubProbExpr='';
NameList = {};
 
% Loop until the string is empty
while (~isempty(Reminder)),
    
    % Collect a possible token consisting of alphanumeric characters
    [AlphaTok , AlphaReminder] = strtok(Reminder,Delimiters);
 
    % Collect a possible token consisting only of delimiters
    [OtherTok , OtherReminder] = strtok(Reminder,AlphaNumeric);
   
    % Since the implementation of strtok does not return position, check 
    % to see what token was first by the length of the reminder
    % Both reminders are empty at the end of the string. 
    if  ( length(AlphaReminder) < length(OtherReminder)) || ( isempty(AlphaReminder) && isempty(OtherReminder) && isempty(AlphaTok)),
        % This means a non alphanumeric token was first
        Tok = OtherTok;
        Reminder = OtherReminder;
    else
        % This means an alphanumeric token was first
        Tok = AlphaTok;
        Reminder = AlphaReminder;
 
        % Now examine the token
        
        % If the token is not a function from the supported function list, 
        % and if it is not a number then it is a token that should be in 
        % the output list. Note that the comparison is case sensitive 
        FunctionNames = {'exp', 'log' , 'log1p', 'log2', 'log10', 'logm', 'reallog' , 'sym'};
        if ~any(strcmpi(Tok,FunctionNames)) && isempty(str2num(Tok)),
            NameList{end+1} = Tok;
 
            % If the index is not empty, token substitution is requested 
            if ~isempty(Index),
                % At this point, it is assumed that the token contains a
                % transition name of the form T#_#_# where the first number
                % corresponds to the study id, the second to the start
                % state and the third to the end state. 
                % First, verify this format and extract these numbers
                ReadTrans = sscanf(Tok,'T%i_%i_%i');
                
                if length(ReadTrans) ~= 3,
                    error ('The term "%s" extracted from expression is not a transition, therefore the use of Index is invalid',Tok);
                end
                
                % Look for this transition from the input transition set.
                StudyCheck = [Transitions.StudyID] == ReadTrans(1);
                FromCheck = [Transitions.From] == ReadTrans(2);
                ToCheck = [Transitions.To] == ReadTrans(3);
                
                % The study index should satisfy all three indices
                TransIndex = find(StudyCheck & FromCheck & ToCheck);
                
                % Data consistancy checks
                if length(TransIndex)>1,
                    error('Input data has more than one transition defined as "%s"',Tok);                    
                end
 
                if isempty(TransIndex),
                    error('no transition is defined that corresponds to "%s"',Tok);                    
                end
                
                % Access the table in the transition with the given indices
                % Note that Index is a subscript index vector that is
                % converted to a single index into a flat vector
 
                FlatIndex = MySub2Ind( size(Transitions(TransIndex).Data), Index );
                CellData = Transitions(TransIndex).Data(FlatIndex);
                
                % Now replace the token with the cell data
                % If the data is a symbolic string - just replace the token
                % be sure to use parentheses
                if isobject(CellData)
                    Tok = ['(' char(CellData) ')'];
                elseif isa(CellData,'numeric')
                    Tok = ['(sym(' num2str(CellData) '))'];
                else
                    % Other types are unsupported, raise an error
                    error ('Unsupported data for token "%s". Data is "%s"',Tok,char(CellData));
                end
            end
        end
    end
 
    % Add the token to the output string
    SubProbExpr = [SubProbExpr char(Tok)];
    
end
 
