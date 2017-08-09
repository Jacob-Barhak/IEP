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
Covariates = {'Gender','Smoke','Race','Age','SBP'};
Coefficients = { 'P01','P02','P12','P14','P23m','P23f','P32','P34'};
InitialGuess = [  0.05,0.032, 0.05, 0.05, 0.05, 0.05, 0.09, 0.07];
SearchGlobal = [ 1e-11 2e-11 5e-11 1e-10 2e-10 5e-10 1e-9 2e-9 5e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 0.1 0.25 0.5 0.75 ];

% Initialize Parameter structure and create symbolic variables
InitParams

%%%% Define Population Set Information
PopColumns1(1)     = struct ('Dim','Gender','Distribution','Ber(0.5)','OrderInSet',0);
Pop(1) = struct ('Name', 'Dummy Population','Columns', PopColumns1 ,'Data',[]);

PopColumns2(1)     = struct ('Dim','Gender','Distribution','','OrderInSet',1);
PopColumns2(end+1) = struct ('Dim','Smoke','Distribution','','OrderInSet',2);
PopColumns2(end+1) = struct ('Dim','Race','Distribution','','OrderInSet',3);
PopColumns2(end+1) = struct ('Dim','Age','Distribution','','OrderInSet',4);
PopColumns2(end+1) = struct ('Dim','SBP','Distribution','','OrderInSet',5);
Pop(2) = struct ('Name', 'Data Population','Columns', PopColumns2 ,'Data',[1 0 0 54 115;1 0 1 56 125;0 1 0 54 115;0 0 0 56 125] );

% Define Studies. The model is the first study record with StudyLength=0
Studies(1) = struct ('Name', 'Model0', 'StudyLength', 0, 'MainProcess', 10, 'PopID', []);
Studies(2) = struct ('Name', 'Study1', 'StudyLength', 10, 'MainProcess', 0, 'PopID', 1);
Studies(3) = struct ('Name', 'Study2', 'StudyLength', 7, 'MainProcess', 0, 'PopID', 1);
Studies(4) = struct ('Name', 'Study3', 'StudyLength', 10, 'MainProcess', 0, 'PopID', 2);
Studies(5) = struct ('Name', 'Study4', 'StudyLength', 2, 'MainProcess', 0, 'PopID', 1);
Studies(6) = struct ('Name', 'Study5', 'StudyLength', 2, 'MainProcess', 0, 'PopID', 1);
Studies(7) = struct ('Name', 'Study6', 'StudyLength', 1, 'MainProcess', 0, 'PopID', 1);
Studies(8) = struct ('Name', 'Study7', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);
Studies(9) = struct ('Name', 'Study8', 'StudyLength', 7, 'MainProcess', 0, 'PopID', 1);
Studies(10) = struct ('Name', 'Study9', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);
Studies(11) = struct ('Name', 'Study10', 'StudyLength', 5, 'MainProcess', 0, 'PopID', 1);

% Define States 
States(1)  = struct ('Name', 'State0', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(2)  = struct ('Name', 'State1', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(3)  = struct ('Name', 'State2', 'IsEvent', 1, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(4)  = struct ('Name', 'State3', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(5)  = struct ('Name', 'State4', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', []);
States(10) = struct ('Name', 'MainProcess', 'IsEvent', 0, 'IsSplit', 0, 'JoinerOfSplit', 0, 'PoolingData', [ 1 2 3 4 5; 0 0 0 0 0] );

% Define Tables
CoefficientTableP23 = struct('Data',[P23f ; P23m],'Dim',[Gender,NaN,0,1]);
CoefficientTableP24 = struct('Data',[1-P23f ; 1-P23m],'Dim',[Gender,NaN,0,1]);
TableStartPopulation2 = struct('Data',[1897; 2643],'Dim',[Gender,NaN,0,1]);
TableObservedCount2 = struct('Data',[78.61168 152.9931951205 223.211967802713 289.373249867647 351.609452192458; 125.17248 243.609316788 355.417866346905 460.765719948907 559.863713921029],'Dim',[Gender,NaN,0,1,Time,NaN,2,4,6,8,10]);
TableStartPopulation4 = struct('Data',[163 ; 312],'Dim',[Gender,NaN,0,1]);
TableObservedCount4 = struct('Data',[48.9 ; 62.4],'Dim',[Gender,NaN,0,1,Time,NaN,1]);
TableObservedCount5 = struct('Data',[1.46;2.8616;6.7388965376],'Dim',[Time,NaN,1,2,5]);
TableStartPopulation8 = struct('Data',[181 ; 287],'Dim',[Gender,NaN,0,1]);
TableObservedCount8 = struct('Data',[22.3378411472491 ; 32.8260810667131],'Dim',[Gender,NaN,0,1,Time,NaN,5]);

% Define Model and Study Transitions
% Model First
Transitions(1) = struct ('StudyID',1,'From', 1, 'To', 2, 'Probability', P01 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 1, 'To', 3, 'Probability', P02 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 2, 'To', 3, 'Probability', P12 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 2, 'To', 5, 'Probability', P14 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 3, 'To', 4, 'Probability', CoefficientTableP23 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 3, 'To', 5, 'Probability', CoefficientTableP24 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 4, 'To', 3, 'Probability', P32 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',1,'From', 4, 'To', 5, 'Probability', P34 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study1
Transitions(end+1) = struct ('StudyID',2,'From', 1, 'To', 0, 'Probability', 1138 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',2,'From', 1, 'To', 2, 'Probability', 48.7255388591079 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study2
Transitions(end+1) = struct ('StudyID',3,'From', 1, 'To', 0, 'Probability', 890 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',3,'From', 1, 'To', 3, 'Probability', 172.110040233589 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study3
Transitions(end+1) = struct ('StudyID',4,'From', 1, 'To', 0, 'Probability', TableStartPopulation2 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',4,'From', 1, 'To', 4, 'Probability', TableObservedCount2 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study4
Transitions(end+1) = struct ('StudyID',5,'From', 2, 'To', 0, 'Probability', 569 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',5,'From', 2, 'To', 3, 'Probability', 54.055 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study5
Transitions(end+1) = struct ('StudyID',6,'From', 2, 'To', 0, 'Probability', 569 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',6,'From', 2, 'To', 5, 'Probability', 68.100765 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study6
Transitions(end+1) = struct ('StudyID',7,'From', 3, 'To', 0, 'Probability', TableStartPopulation4 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',7,'From', 3, 'To', 5, 'Probability', TableObservedCount4 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study7
Transitions(end+1) = struct ('StudyID',8,'From', 4, 'To', 0, 'Probability', 73 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',8,'From', 4, 'To', 3, 'Probability', TableObservedCount5 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study8
Transitions(end+1) = struct ('StudyID',9,'From', 4, 'To', 0, 'Probability', 169 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',9,'From', 4, 'To', 3, 'Probability', 21.0026880998605 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study9
Transitions(end+1) = struct ('StudyID',10,'From', 4, 'To', 0, 'Probability', 78 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',10,'From', 4, 'To', 3, 'Probability', 7.2004647936 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
%Study10
Transitions(end+1) = struct ('StudyID',11,'From', 4, 'To', 0, 'Probability', TableStartPopulation8 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );
Transitions(end+1) = struct ('StudyID',11,'From', 4, 'To', 5, 'Probability', TableObservedCount8 ,'RegCoefVec', [] ,'RegParamVec', [] , 'CovarianceMatrix', [] , 'RegFunctionType', []  );

% Build the Likelihood expression by calling a script
CalcLikelihood

% call the optimize Likelihood routine to extract the value
OptimizeLikelihood
