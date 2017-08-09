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
% OptimizeLikelihood.m
% This script optimizes the Likelihood function L and outputs the
% optimization results to screen and to memory.
% This script should be called after the CalcLikelihood.m script that 
% constructs the likelihood expressions. Also additional input information
% is required beyond the Likelihood expression and should be defined in the
% input file such as intial guesses information and coefficeint names.
%
% REQUIRED INPUT:
% L - The Likelihood function to be optimized
% Coefficients - All the coefficient names as stings in a cell vector
% InitialGuess - The initial guess matrix, where each row represents a 
%                vector of guesses and each colum corresponds to a
%                different coefficeint defined in Coefficients
% SearchGlobal - A set of numbers that provides complimentary initial
%                guesses. Each number x in that vector is used to represent
%                an intial guess that all its coefficient values are that
%                number i.e. the initial guess can be defined as 
%                x*ones(1,length(Coefficeints))
%
% OUTPUT:
% VarVecOriginalName - a cell vector with strings representing the
%                      coefficient names optimized .
% BestX - The best optimization result vector. The members of this vector
%         correspond to the coefficients named in VarVecOriginalName. This
%         solution is the best solution for all initial guesses and defined
%         in the input. The criteria for selecting the best solution is
%         being a valid solution with the minimum minus log likelihood
% BestFuncVal - The minus minimum log likelihood obtained for the best
%               optimized solution.
% BestVariance - The Variance associated with the Best solution as it is
%                reported by the optimizer. It is calculated by inverting
%                the Hessian at the solution point
% AllResults - A vector containing the optimization results for all initial
%              Guesses. Each row in this vector contains the following
%              vector [BestSolution MinusLogLikelihood Variance(:)'].
%              Each row corresponds to an initial solution. The first rows
%              correspond to rows in the InitialGuess Vector. Then each row
%              corresponds to a single member from the SearchGlobal vector.
%              The best solution is indexed by BestGuessNum that also
%              appears in the output.


% Reset timer
tic
disp '===================================================='
disp '===================================================='
disp 'Optimizing the Log likelihood Expression'
 
% Set warning behavior
warning off all
lastwarn('')
 
LikelihoodToOptimize = -(L);
 
[ LikelihoodFuncHandle, VarVecOriginalName, ReplaceVarVec, InitialGuessModified, LowerBound, UpperBound ] = PrepareForNumericSolution (LikelihoodToOptimize, Coefficients, InitialGuess, SearchGlobal);
 
% Set optimization options
FunctionReplaceTol = 1e-11;
Options = optimset('fminunc'); 
Options = optimset(Options,'Display', 'iter');
Options = optimset(Options,'TolFun',1e-10);
Options = optimset(Options,'TolX',1e-10);
Options = optimset(Options,'FunValCheck', 'on');
Options = optimset(Options,'GradObj','on');
 
 
% Now loop through all initial guesses and try to optimize. If not successful
% due to an error try the next option. Each time record the results and
% finally select the best result from all attempts. 
BestX = nan;
BestFuncVal = inf;
BestGuessNum = nan;
BestVariance = nan;
 
AllResults=[];
 
for Guess = 1:size(InitialGuessModified,1),
 
    x0 = InitialGuessModified(Guess,:);
    
    % Try optimizing the function
    try 
        ReplacementFunctionReturns = LikelihoodFuncHandle (x0) ;
        OriginalFunctionReturns = subs (LikelihoodToOptimize,VarVecOriginalName,x0);
        if ~(abs(OriginalFunctionReturns - ReplacementFunctionReturns) < FunctionReplaceTol),
            Text = 'Initial guess does not produce similar results for the original and the replacement functions';
            disp (Text)
            error (Text)
        end
        [x,FuncVal,ExitFlag,Output,Lambda,GradRes,HessianRes] = fmincon(LikelihoodFuncHandle,x0,[],[],[],[],LowerBound,UpperBound,[],Options);
 
        ReplacementFunctionReturns = LikelihoodFuncHandle (x) ;
        OriginalFunctionReturns = subs (LikelihoodToOptimize,VarVecOriginalName,x);
        if FuncVal ~= ReplacementFunctionReturns,
            Text = 'Assertion Error - return value should be the same as calculated values';
            disp (Text)
            error (Text)
        end
        if ~(abs(OriginalFunctionReturns - ReplacementFunctionReturns) < FunctionReplaceTol),
            Text = 'Optimized value does not produce similar results for the original and the replacement functions';
            disp (Text)
            error (Text)
        end
       
        Variance = inv(full(HessianRes));
        % Store the results for 
        AllResults(Guess,:) = [x OriginalFunctionReturns Variance(:)'];
    catch
        % If unsuccesful, and an error was raised print an error message
        % and continue to the next guess in the loop
        disp ('The optimization failed for the following initial guess:');
        for i = 1:length(x0),
            disp ([VarVecOriginalName{i} ' = ' num2str(x0(i))]);
        end
        disp ('***Since an error was detected for this guess, no result is valid***');
        AllResults(Guess,:) = [x0*NaN NaN ones(1,length(x0)^2)*NaN];
        continue
    end
 
    
    % If received better results than the current best, note these as best. 
    if OriginalFunctionReturns<BestFuncVal,
        BestX = x;
        BestFuncVal = OriginalFunctionReturns;
        BestGuessNum = Guess;
        BestVariance = Variance;
    end
 
    disp (sprintf('The optimization initial guess was:'));
    for i = 1:length(x),
        disp ([VarVecOriginalName{i} ' = ' num2str(x0(i))]);
    end
 
    disp (sprintf('The Minimum Minus Log-Likelihood value is: %g', OriginalFunctionReturns));
 
    disp (sprintf('The Variance for the solution is:\n %s', mat2str(Variance)));
    
    disp (sprintf('The optimization results are:'));
    for i = 1:length(x),
        disp ([VarVecOriginalName{i} ' = ' num2str(x(i))]);
    end
 
    
end
    
disp (sprintf('-------------------------------------------------'));
 
if isnan(BestGuessNum),
    disp (sprintf('Optimization Failure for all initial guesses - Check the optimization function'));
else
 
    
    disp (sprintf('The Best optimization result was received for the following guess:'));
    for i = 1:length(BestX),
        disp ([VarVecOriginalName{i} ' = ' num2str(InitialGuessModified(BestGuessNum,i))]);
    end
 
    disp (sprintf('The Final Minimum Minus Log-Likelihood value is: %g', BestFuncVal));
 
    disp (sprintf('The Final Variance is:\n %s', mat2str(BestVariance)));
    
    disp (sprintf('The final optimization result is:'));
    for i = 1:length(BestX),
        disp ([VarVecOriginalName{i} ' = ' num2str(BestX(i))]);
    end
    
end
 
if DebugVec(2), 
    disp 'Other Optimization Related Information:'
    % Show the time required for optimization
    disp 'Optimization time was:'
    toc
end

disp '===================================================='
disp '===================================================='
 
