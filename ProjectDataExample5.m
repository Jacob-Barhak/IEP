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
Covariates = {};
Coefficients = { 'P01', 'P02', 'P12', 'P23', 'P34', 'P35', 'P45'};
InitialGuess = [ 0.04 , 0.04 , 0.1  , 0.015, 0.04 , 0.1  , 0.04 
                 0.05 , 0.03 , 0.1  , 0.01 , 0.04 , 0.11 , 0.04 
                 0.046, 0.026, 0.091, 0.017, 0.039, 0.119, 0.045];
SearchGlobal = [ 1e-11 2e-11 5e-11 1e-10 2e-10 5e-10 1e-9 2e-9 5e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 0.1 0.25 0.5 0.75 ];

% Initialize Parameter structure and create symbolic variables
InitParams

%%%% Define Population Set Information
PopColumns1(1) = struct ('Dim','Default','Distribution','','OrderInSet',0);
Pop(1) = struct ('Name', 'Default Population','Columns', PopColumns1 ,'Data',[0,1]);

% Define Studies. The model is the first study record with StudyLength=0
Studies(1) = struct ('Name', 'Model0', 'StudyLength', 0, 'MainProcess', 10, 'PopID', []);
Studies(2) = struct ('Name', 'Study1', 'StudyLength', 6, 'MainProcess', 0, 'PopID', 1);
Studies(3) = struct ('Name', 'Study2', 'StudyLength', 9, 'MainProcess', 0, 'PopID', 1);
Studies(4) = struct ('Name', 'Study3', 'StudyLength', 4, 'MainProcess', 0, 'PopID', 1);
Studies(5) = struct ('Name', 'Study4', 'StudyLength', 6, 'MainProcess', 0, 'PopID', 1);
Studies(6) = struct ('Name', 'Study5', 'StudyLength', 6, 'MainProcess', 0, 'PopID', 1);
Studies(7) = struct ('Name', 'Study6', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);
Studies(8) = struct ('Name', 'Study7', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);
Studies(9) = struct ('Name', 'Study8', 'StudyLength', 1, 'MainProcess', 0, 'PopID', 1);
Studies(10) = struct ('Name', 'Study9', 'StudyLength', 3, 'MainProcess', 0, 'PopID', 1);
Studies(11) = struct ('Name', 'Study10', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);
Studies(12) = struct ('Name', 'Study11', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);

% Define States 
States(1)  = struct ('Name', 'State0', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(2)  = struct ('Name', 'State1', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(3)  = struct ('Name', 'State2', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(4)  = struct ('Name', 'State3', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(5)  = struct ('Name', 'State4', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(6)  = struct ('Name', 'State5', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
%
States(9)  = struct ('Name', 'PooledState(0,1)', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', [1 2; 0.5 0.5]);
%
States(10) = struct ('Name', 'MainProcess', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', [ 1 2 3 4 5 6; 0 0 0 0 0 0] );

% Define Model and Study Transitions
Transitions(1) = struct ('StudyID',1,'From', 1, 'To', 2, 'Probability', P01 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 2, 'To', 3, 'Probability', P12 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 3, 'To', 4, 'Probability', P23 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 4, 'To', 5, 'Probability', P34 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 4, 'To', 6, 'Probability', P35 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 5, 'To', 6, 'Probability', P45 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 1, 'To', 3, 'Probability', P02 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', []  , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',2,'From', 1, 'To', 0, 'Probability', 79 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',2,'From', 1, 'To', 2, 'Probability', 12.871600960512 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',3,'From', 1, 'To', 0, 'Probability', 90 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',3,'From', 1, 'To', 2, 'Probability', 20.7539552818476 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',4,'From', 9, 'To', 0, 'Probability', 398 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',4,'From', 9, 'To', 3, 'Probability', 79.18810184 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',5,'From', 1, 'To', 0, 'Probability', 176 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',5,'From', 1, 'To', 3, 'Probability', 16.118751520256 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',6,'From', 2, 'To', 0, 'Probability', 49 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',6,'From', 2, 'To', 3, 'Probability', 22.959391 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',7,'From', 2, 'To', 0, 'Probability', 45 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',7,'From', 2, 'To', 3, 'Probability', 18.42795 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',8,'From', 3, 'To', 0, 'Probability', 202 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',8,'From', 3, 'To', 4, 'Probability', 19.4079990464 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',9,'From', 4, 'To', 0, 'Probability', 1000 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',9,'From', 4, 'To', 5, 'Probability', 40 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',10,'From', 4, 'To', 0, 'Probability', 231 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',10,'From', 4, 'To', 6, 'Probability', 109.24914 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',11,'From', 4, 'To', 0, 'Probability', 11929 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',11,'From', 4, 'To', 6, 'Probability', 7597.886603726 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%
Transitions(end+1) = struct ('StudyID',12,'From', 5, 'To', 0, 'Probability', 23 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',12,'From', 5, 'To', 6, 'Probability', 5.2030384375 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );

% Build the Likelihood expression by calling a script
CalcLikelihood

% call the optimize Likelihood routine to extract the value
OptimizeLikelihood
