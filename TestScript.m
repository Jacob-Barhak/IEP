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
% This script runs all test examples provided consecutively and records the
% output transcript in the file TestingOutput.txt
% Use this script as a test to verify that the system is working properly
% after installation or after any modification.
% For additional details, about the examples see the following document:
% ExamplesForIndirectEstimation.pdf


clear;
clc;

diary off
delete TestingOutput.txt
diary TestingOutput.txt

%Create results section
Results = {};


% Run the following examples
ExpectedResults = {};
ExampleFileNames = {};
ExampleFileNames{end+1}='ProjectDataExample1';
ExpectedResults{end+1} = [0.0716822332774442];
ExampleFileNames{end+1}='ProjectDataExample2';
ExpectedResults{end+1} = [0.1 0.2];
ExampleFileNames{end+1}='ProjectDataExample3';
ExpectedResults{end+1} = [0.1 0.2];
ExampleFileNames{end+1}='ProjectDataExample4';
ExpectedResults{end+1} = [0.1 0.2 0.3];
ExampleFileNames{end+1}='ProjectDataExample5';
ExpectedResults{end+1} = [0.03 0.01 0.1 0.02 0.04 0.2 0.05];
ExampleFileNames{end+1}='ProjectDataExample6';
ExpectedResults{end+1} = [0.0967842995713971];
ExampleFileNames{end+1}='ProjectDataExample7';
ExpectedResults{end+1} = [0.3 0.4];
%ExampleFileNames{end+1}='ProjectDataExample7_OnlyOneStudy';
%ExpectedResults{end+1} = [0.3 0.4]; % multiple possible solutions
ExampleFileNames{end+1}='ProjectDataExample8';
ExpectedResults{end+1} = [0.126419535263701 0.0716822332774442];
ExampleFileNames{end+1}='ProjectDataExample9';
ExpectedResults{end+1} = [0.170183283369612];
ExampleFileNames{end+1}='ProjectDataExample9_OneYear';
ExpectedResults{end+1} = [0.90483741803596];
ExampleFileNames{end+1}='ProjectDataExample10_OneYear';
ExpectedResults{end+1} = [0.0023133231471746951	0.0018939888023991047];
ExampleFileNames{end+1}='ProjectDataExample10_TwoYears';
ExpectedResults{end+1} = [2.8362402018089483e-06 1.9011877737673544e-06];
ExampleFileNames{end+1}='ProjectDataExample10_ThreeYears';
ExpectedResults{end+1} = [4.8691513043763734e-09 2.672246868229422e-09];
ExampleFileNames{end+1}='ProjectDataExample11_OneYear';
ExpectedResults{end+1} = [0.50325950521854423	0.59330328062324944];
ExampleFileNames{end+1}='ProjectDataExample11_TwoYears';
ExpectedResults{end+1} = [0.5031043626979157	0.59317626067050078];
ExampleFileNames{end+1}='ProjectDataExample11_ThreeYears';
ExpectedResults{end+1} = [0.50294936539903712	0.59304935961526495];
ExampleFileNames{end+1}='ProjectDataExample12_OneYear';
ExpectedResults{end+1} = [0.00597431015512395 0.00653004287191838];
ExampleFileNames{end+1}='ProjectDataExample12_TwoYears';
ExpectedResults{end+1} = [0.00620653272582317 0.00678083530428653];
ExampleFileNames{end+1}='ProjectDataExample12_ThreeYears';
ExpectedResults{end+1} = [0.00645076680543022 0.007044275112843];
ExampleFileNames{end+1}='ProjectDataExample13';
ExpectedResults{end+1} = [0.2 0.3];
ExampleFileNames{end+1}='ProjectDataExample14';
ExpectedResults{end+1} = [0.2 0.3];
ExampleFileNames{end+1}='ProjectDataExample15';
ExpectedResults{end+1} = [0.2 0.3];
ExampleFileNames{end+1}='ProjectDataExample16';
ExpectedResults{end+1} = [0.2 0.3];
ExampleFileNames{end+1}='ProjectDataExample17';
ExpectedResults{end+1} = [0.2 0.3 0.4];
ExampleFileNames{end+1}='ProjectDataExample17_OneYear';
ExpectedResults{end+1} = [0.2 0.3 0.4];
ExampleFileNames{end+1}='ProjectDataExample18';
ExpectedResults{end+1} = [0.2 0.3 0.4];
ExampleFileNames{end+1}='ProjectDataExample18_OneYear';
ExpectedResults{end+1} = [0.2 0.3 0.4];
ExampleFileNames{end+1}='ProjectDataExample19';
ExpectedResults{end+1} = [0.2 0.3 0.4 0.5];
ExampleFileNames{end+1}='ProjectDataExample19_OneYear';
ExpectedResults{end+1} = [0.2 0.3 0.4 0.5];
ExampleFileNames{end+1}='ProjectDataExample20';
ExpectedResults{end+1} = [0.2 0.3 0.4 0.5 0.6];
ExampleFileNames{end+1}='ProjectDataExample20_OneYear';
ExpectedResults{end+1} = [0.2 0.3 0.4 0.5 0.6];
ExampleFileNames{end+1}='ProjectDataExample21';
ExpectedResults{end+1} = [0.1];
ExampleFileNames{end+1}='ProjectDataExample22';
ExpectedResults{end+1} = [0.1 0.2 0.4 0.3];
ExampleFileNames{end+1}='ProjectDataExample23';
ExpectedResults{end+1} = [0.005 0.03 0.05 0.05 0.75 0.02 0.02];
ExampleFileNames{end+1}='ProjectDataExample24';
ExpectedResults{end+1} = [0.005 0.03 0.05 0.05 0.75 0.02 0.02];
ExampleFileNames{end+1}='ProjectDataExample25';
ExpectedResults{end+1} = [0.005 0.03 0.05 0.05 0.75 0.02 0.02];
ExampleFileNames{end+1}='ProjectDataExample26';
ExpectedResults{end+1} = [0.005 0.03 0.05 0.05 0.7 0.8 0.02 0.02];
ExampleFileNames{end+1}='ProjectDataExample27';
ExpectedResults{end+1} = [0.005 0.025 0.035 0.05 0.05 0.7 0.8 0.02 0.015 0.025];



IndexToFilesToRun = 0;
while IndexToFilesToRun<length(ExampleFileNames),
    IndexToFilesToRun = IndexToFilesToRun+1;
    save TempFile IndexToFilesToRun ExampleFileNames Results ExpectedResults
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************
    disp (['   SCRIPT FILE NAME - ' ExampleFileNames{IndexToFilesToRun}])
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************
    disp *******************************************************

    % Clear the maple work space
    % maple restart
    % Run the script
    eval (ExampleFileNames{IndexToFilesToRun});

    % Since data was cleared reload old data and continue
    load TempFile

    % Calculate the error from the expected sulotion 
    ErrorFromExpectedSolution = BestX - ExpectedResults{IndexToFilesToRun};
    disp 'Error From The Expected Solution is:'
    disp (mat2str(ErrorFromExpectedSolution))
    disp '#######################################################'
    disp '#######################################################'
    % stores the results back
    Results{IndexToFilesToRun,1} = VarVecOriginalName;
    Results{IndexToFilesToRun,2} = BestX;
    Results{IndexToFilesToRun,3} = BestFuncVal;
    Results{IndexToFilesToRun,4} = BestVariance;    
    Results{IndexToFilesToRun,5} = ExampleFileNames{IndexToFilesToRun};
    Results{IndexToFilesToRun,6} = ErrorFromExpectedSolution;
    
end

diary off

save ('TestResultsMatFile','Results')

% Collect the errors and show results in a plot
Errors = [];
for i = 1:length(Results),
    Errors = [Errors Results{i,6}];
end
ErrorsSorted = sort(abs(Errors));
CalcCDF = [];
i=0;
for XVal = ErrorsSorted,
    i=i+1;
    CalcCDF(i) = sum(XVal>=ErrorsSorted)/length(ErrorsSorted);
end
figure;
semilogx(ErrorsSorted,CalcCDF,'-','LineWidth',4)
title('Cumulative Frequency Distribution Plot of the Absolute Error')
xlabel('Absolute Error')
ylabel('Observed Cumulative Probability')

