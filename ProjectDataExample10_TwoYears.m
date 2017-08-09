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
% This script contains an example to demonstrate the Indirect Estimation
% Prototype. For more information on this example, see the document
% ExamplesForIndirectEstimation.pdf for more information regarding this
% example.
% This script calls the following other scripts: InitParams.m,
% CalcLikelihood.m, and OptimizeLikelihood.m to generate the results. After
% running this script, it is possible to call LikelihoodGraphs.m to
% visualize the optimization results.

% Clear all variables
clear;

% Define parameters
Constants = {};
ConstantValues = [];
% Note that all index values must be defined for indexed params
Covariates = {'Gender','BP'};
Coefficients = { 'P01m','P01f'};
InitialGuess = [ 0.1 , 0.1 ];
SearchGlobal = [ 1e-11 2e-11 5e-11 1e-10 2e-10 5e-10 1e-9 2e-9 5e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 0.1 0.25 0.5 0.75 ];

% Initialize Parameter structure and create symbolic variables
InitParams

%%%% Define Population Set Information
PopColumns1(1) = struct ('Dim','Gender','Distribution','','OrderInSet',1);
PopColumns1(end+1) = struct ('Dim','BP','Distribution','','OrderInSet',2);
Pop(1) = struct ('Name', 'Data Population','Columns', PopColumns1 ,'Data',[1, 115; 0, 125; 0, 115; 1, 125] );

% Define the study length to solve the example for
StudyLengthConst=2;

% Define Studies. The model is the first study record with StudyLength=0
Studies(1) = struct ('Name', 'Model0', 'StudyLength', 0, 'MainProcess', 10, 'PopID', []);
Studies(2) = struct ('Name', 'Study1', 'StudyLength', StudyLengthConst, 'MainProcess', 0, 'PopID', 1);

% Define States 
States(1)  = struct ('Name', 'State0', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(2)  = struct ('Name', 'State1', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
%
States(10) = struct ('Name', 'MainProcess', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', [ 1 2 ; 0 0] );

% Define Tables
CoefficientTable1 = struct('Data',[P01f ; P01m],'Dim',[Gender,NaN,0,1]);

% Define Model and Study Transitions
Transitions(1) = struct ('StudyID',1,'From', 1, 'To', 2, 'Probability', CoefficientTable1 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',2,'From', 1, 'To', 2, 'Probability', [] ,'RegCoefVec', [0.1 0.2 0.05] ,'RegParamVec', StudyLengthConst*[1 Gender BP] , 'CovarianceMatrix', eye(3) , 'RegFunctionType', 1  );

% Build the Likelihood expression by calling a script
CalcLikelihood

% call the optimize Likelihood routine to extract the value
OptimizeLikelihood
