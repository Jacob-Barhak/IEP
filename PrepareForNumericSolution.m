function [ CombinedFuncHandle, VarVecOriginalName, ReplaceVarVec, InitialGuessModified, LowerBound, UpperBound ] = PrepareForNumericSolution (FunctionToPrepare, Coefficients, InitialGuess, SearchGlobal)
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
%function [ CombinedFuncHandle, VarVecOriginalName, ReplaceVarVec, InitialGuessModified, LowerBound, UpperBound ] = PrepareForNumericSolution (FunctionToPrepare, Coefficients, InitialGuess, SearchGlobal)
% The function prepares a function for numerical calculations by replacing
% variable names and returning a new function:
% Input: 
% FunctionToPrepare: The function to prepare using symbolic variables
% Coefficients: a vector of coefficients of interest
% InitialGuess: a vector of initial guesses corresponding to Coefficients
% SearchGlobal: a vector with initial guesses to be used for all parameters
% Output:
% CombinedFuncHandle: The modified function handle along with its Jacobian
% VarVecOriginalName: The original names of the variables in a vector
% ReplaceVarVec: The replacement names for the variables in a vector
% InitialGuessModified: A corresponding matrix of initial guesses 
% LowerBound: The lower bound for the parameters (trivially supported)
% UpperBound: The upper bound for the parameters (trivially supported)
 
% This is a constant that can be later changed
MinxAllowed = 1e-11;
 
% First extract the name of the variables
InputVarVec = symvar(char(FunctionToPrepare));
AllVarNum = length(InputVarVec);
 
% Convert the cell array to a single symbolic vector
SymVarVec = sym(zeros(AllVarNum,1));
for i=1:AllVarNum
   SymVarVec(i) = sym(char(InputVarVec(i)));
end
 
% Symbolically calculate the gradient vector / Jacobian matrix
PreparedGradient = jacobian(FunctionToPrepare,SymVarVec);
 
% Exchange the variables with the vector x
% At the same time calculate the bounds and construct the initial guess
% vector
ReplaceVarVec = {};
VarVec = {};
VarVecOriginalName = {};
LowerBound=[];
UpperBound=[];
InitialGuessModified=[];
% Define a set of initial guesses to start optimization
% Here we choose a small number close to zero, we increase these in orders
% of magnitude and scan the range up to 1. 
 
for i = 1:AllVarNum,
    % If this is a coefficient to be optimized provide it with a new name
    % otherwise raise an error
    OriginalName = char(InputVarVec(i));
    [IsInCoefficients LocInCoeficients] = ismember(OriginalName,Coefficients);
    if IsInCoefficients,
        % If the user defined the variable
        Index = length(VarVec)+1;
        VarName = ['x(' num2str(Index) ')'];
        VarVec{Index} = VarName;
        VarVecOriginalName{Index} = OriginalName;
        
        if ~isempty(InitialGuess),
            InitialGuessModified(1:size(InitialGuess,1),Index) = InitialGuess(:,LocInCoeficients);
        end
        % If the user specified values for global search, then
        if ~isempty(SearchGlobal),
            InitialGuessModified((size(InitialGuess,1)+1):(size(InitialGuess,1)+length(SearchGlobal)),i) = SearchGlobal;
        end
        % Probability coefficients will be bound between 0 an 1
        % Note that this is not evolved yet and may require more coding
        % For now, when the parameter name starts with 'P' then impose
        % probability boundaries of (0,1], otherwise use positive bounds 
        Low = MinxAllowed;
        if lower(OriginalName(1)) == 'p',
            High = 1;
        else
            High = inf;
        end
        LowerBound(Index)=Low;
        UpperBound(Index)=High;
        % Update the initial guess to be within bound limits
        InitialGuessModified(:,i)= min(max(InitialGuessModified(:,i),Low),High);
    else
        error('The coefficient %s is not registered as a coefficient', char(InputVarVec(i)))
    end
    % In the future, this may be used to incorporate bounds for the solver
    ReplaceVarVec{i} = VarName;
    LowerBound(i) = Low;
    UpperBound(i) = High;
end
 
Func = subs(FunctionToPrepare,InputVarVec,ReplaceVarVec);
JacobianMat = subs(PreparedGradient,InputVarVec,ReplaceVarVec);
 
% Replace the log function with the MyLog version of it to create a
% modified function for optimization. 
ModFunc = strrep(char(Func),'log','MyLog');
ModJacobianMat = strrep(char(sym(JacobianMat)),'log','MyLog');
 
 
% Take care of the matrix notation that will appear while converting the
% symbolic variable to a string
if ~isempty(ModFunc) && strcmp(ModFunc(1:min(length(ModFunc),7)),'matrix('),
    % Note that the double brackets are stripped as well as the matrix
    % leading brackets
    TruncMod = ModFunc(8:end-1);
else
    TruncMod = ModFunc;
end
 
 
if ~isempty(ModJacobianMat) && strcmp(ModJacobianMat(1:min(length(ModJacobianMat),7)),'matrix('),
    % Note that the double brackets are stripped as well as the matrix
    % leading brackets
    TruncModJacobianMat = ModJacobianMat(8:end-1);
else
    TruncModJacobianMat = ModJacobianMat;
end
 
CommandStr = ['@(x) PrintFuncValueDiagnostics(x,'  TruncMod , ',''X='','',F='',true)'];
FuncHandle = eval(CommandStr);
 
CommandStr = ['@(x) PrintFuncValueDiagnostics(x,'  TruncModJacobianMat , ',''X='','',Jacobian='',true)'];
JacobianHandle = eval(CommandStr);
 
CombinedFuncHandle = @(x) MyDeal(FuncHandle(x),JacobianHandle(x));
 
return
 

           

