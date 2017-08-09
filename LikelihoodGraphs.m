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
% LikelihoodGraphs.m
% Show results graphically after optimization is finished
% The function displays the Likelihood graphs near the optimum. Since this
% is a multivariate function with possibly more than 2 dimensions, then
% only one or two dimensions are shown at once, while setting all other
% dimensions 
% to the optimum values.
%
% INPUT:
% L - The Likelihood function to be optimized
% VarVecOriginalName - a cell vector with strings representing the
%                      coefficient names optimized .
% BestX - The best optimization result vector. The members of this vector
%         correspond to the coefficients named in VarVecOriginalName. This
%         solution is the best solution for all initial guesses and defined
%         in the input. The criteria for selecting the best solution is
%         being a valid solution with the minimum minus log likelihood
%
% OUTPUT:
% 1) Graphs displayed in different figure windows
% 2) If requested, PDF files containing these graphs
 
% Configuration parameters - Change these to control display
 
% ShowAllCombinations:
% A Booleean. If true, then graphs are generates for all possible
% combinations of pairs of variables. If false, then a limited number of
% graphs are created that considers only consecutive pairs of variables in
% a cyclic list.
ShowAllCombinations = false;
 
% UseSurfShow
% An integer defining the graph type to use from the following:
% 1 = Mesh, 2 = Surface, 3 = Surface with contour plots
UseSurfShow = 2;
 
% TransperancyValue
% the transperancy of surfaces from 0 to 1, where 0.0 is totaly tranperant
TransperancyValue = 0.5;
 
% FilePrefix
% The string representing the filename prefix for the output pdf filenames
% containing the graph
FilePrefix='MockupScaledDownCovMat';
 
% RangeWidth
% Determines the local range around the optimum value for which to show
% results on the graph. This number is subtracted and added to the optimum
% values to determine the axis min and max bounds for each dimension.
RangeWidth = 0.1;
 
% MinBound and MaxBound
% Minimus and maximum bound for each axis. Overrides RangeWidth in case of
% a conflict.
MinBound = 1e-20;
MaxBound = 1;
 
% StartGraphVersion , EndGraphVersion
% The system supports graphs calculated via one of the following:
% 1) symbolic calculations, 2) VPA precision, 3) Numeric calculations.
% Since in some cases, symbolic calculations may not be able to produce
% graphs properly, several techniques to calculate the graphs may be
% attempted until successful according to the order above.
% StartGraphVersion and EndGraphVersion indicate at which graph version to
% start and at which to end the attempts.
StartGraphVersion = 1;  
EndGraphVersion = 3;  
 
% PrecisionOfVPA
% Indicates the VPA precision if case of Graph version 2
PrecisionOfVPA = 32;
 
% NumberOfTicks
% A number representing the number of mesh elements in one dimension. Note
% that this number is used only if a numeric graph is generated using Graph
% Version 3.
NumberOfTicks = 50; 
 
 
NumberOfVariables = length(VarVecOriginalName);
NumberOfVariablesToSet = max(NumberOfVariables-2,0);
NumberOfVariablesToShow = NumberOfVariables - NumberOfVariablesToSet;
 
if NumberOfVariablesToSet==0,
    ChooseIndices = [];
else
    % check user preference as to how many graphs to show
    if ShowAllCombinations,
        % If the user requests all possible combinations, create these
        ChooseIndices = nchoosek(1:NumberOfVariables,NumberOfVariablesToSet);
    else
        % Otherwise, if the user requests a small number of graphs, create 
        % a graph where consecutive pairs of parameters are displayed
        ChooseIndices = [];
        for i = 1:NumberOfVariablesToSet+1
            ChooseIndices(end+1,:) = [1:(i-1), (i+NumberOfVariablesToShow): NumberOfVariables];        
        end % for i
        % complete the cyclic cycle for two parameters and more
        if NumberOfVariablesToShow == 2,
            ChooseIndices(end+1,:) = 2:(NumberOfVariables-1);        
        end        
    end
end
 
NumberOfGraphs = max(1,size(ChooseIndices,1));
for GraphNum = 1:NumberOfGraphs,
    GraphVersion = StartGraphVersion;
    while GraphVersion <= EndGraphVersion && GraphVersion >0,
        if length(ChooseIndices)>0,
            VaraiblesToSet = VarVecOriginalName(ChooseIndices(GraphNum,:));
            ValuesToSetTo = BestX(ChooseIndices(GraphNum,:));
            SubstitutedLikelihood = subs(L,VaraiblesToSet,ValuesToSetTo);
        else
            SubstitutedLikelihood = L;
        end
        % if using VPA, change this now
        if GraphVersion == 2,
            SubstitutedLikelihood = vpa(SubstitutedLikelihood,PrecisionOfVPA);
        end
        VariablesToBeShown = symvar(char(SubstitutedLikelihood));
        ValuesOfPeak = [];
        EasyShowRange = [];
        ShowRange = [];
        SymVariablesToBeShown = sym([]);
        IndicesOfVarsToShow = [];
        for VarToShowIndex = 1:length(VariablesToBeShown),
            VarToShow = VariablesToBeShown(VarToShowIndex);
            IndexOfVar = strmatch(VarToShow, VarVecOriginalName, 'exact');
            IndicesOfVarsToShow = [IndicesOfVarsToShow IndexOfVar];
            SymVariablesToBeShown(end+1) = sym(char(VarToShow));
            Peak = BestX(IndexOfVar);
            ValuesOfPeak(end+1) = Peak;
            MinRange = max(min(Peak-RangeWidth,MaxBound),MinBound);
            MaxRange = max(min(Peak+RangeWidth,MaxBound),MinBound);
            ShowRange = [ShowRange ; MinRange : ((MaxRange-MinRange)/NumberOfTicks) : MaxRange ];
            EasyShowRange = [EasyShowRange MinRange MaxRange];
        end % for VarToShowIndex
 
        TitleString = char(SymVariablesToBeShown);
        if strcmp(TitleString(1:min(length(TitleString),7)),'matrix('),
            % If a matrix strip the leading and the trailing matrix signs
            % and brackets
            TitleString = TitleString(9:end-2);
        end
        TitleString = ['Log(L(' TitleString '))' ];
        disp (TitleString)
 
 
        if GraphVersion<=2,
            PeakValue = subs(SubstitutedLikelihood,SymVariablesToBeShown,ValuesOfPeak);
            disp (SubstitutedLikelihood)
        else
            PeakValue = -LikelihoodFuncHandle(BestX);
        end
 
        figure
        % If only one variable exits
        if length(IndicesOfVarsToShow) == 1,
            if GraphVersion<=2,
                try
                    ezplot(SubstitutedLikelihood, EasyShowRange)
                catch
                    GraphVersion = GraphVersion +1;
                    disp 'Cannot Show due to numerical problems - trying another version'
                    close
                    continue
                end
            else
                LikelihoodVec = [];
                AxisTickX = [];
                for i = 1:length(ShowRange),
                    ValX = ShowRange(i);
                    SubsVec = BestX;
                    SubsVec(IndicesOfVarsToShow) = ValX;
                    AxisTickX(i) = ValX;
                    LikelihoodVec(i) = -LikelihoodFuncHandle(SubsVec);
                end % for i
                plot(AxisTickX, LikelihoodVec)
                xlabel (char(SymVariablesToBeShown(1)));
                ylabel ('Log L');
            end
            Title (TitleString);
            hold on
            plot(ValuesOfPeak(1), 0 , ValuesOfPeak(1), PeakValue, 'k:o', 'MarkerSize',6', 'LineWidth', 1.5)
            hold off
        else
            if GraphVersion<=2,
                try
                    if UseSurfShow == 1,
                        ezsurf(SubstitutedLikelihood, EasyShowRange);
                    elseif UseSurfShow == 2,
                        ezsurfc(SubstitutedLikelihood, EasyShowRange);
                    else
                        ezmesh (SubstitutedLikelihood, EasyShowRange);
                    end
                catch
                    GraphVersion = GraphVersion +1;
                    disp 'Cannot Show due to numerical problems - trying another version'
                    close
                    continue
                end                    
            else
                LikelihoodVec = [];
                AxisTickX = [];
                AxisTickY = [];
                for i = 1:length(ShowRange(1,:)),
                    for j = 1:length(ShowRange(2,:)),
                        ValX = ShowRange(1,i);
                        ValY = ShowRange(2,j);
                        SubsVec = BestX;
                        SubsVec(IndicesOfVarsToShow(1)) = ValX;
                        SubsVec(IndicesOfVarsToShow(2)) = ValY;
                        AxisTickX(i,j) = ValX;
                        AxisTickY(i,j) = ValY;
                        LikelihoodVec(i,j) = -LikelihoodFuncHandle(SubsVec);
                    end %for j
                end % for i
                if UseSurfShow == 1,
                    surf(AxisTickX, AxisTickY, LikelihoodVec);
                elseif UseSurfShow == 2,
                    surfc(AxisTickX, AxisTickY, LikelihoodVec);
                else
                    mesh (AxisTickX, AxisTickY, LikelihoodVec);
                end
                xlabel (char(SymVariablesToBeShown(1)));
                ylabel (char(SymVariablesToBeShown(2)));
                zlabel ('Log L');
            end
            h = findobj(gca,'Type','surface');
            set(h,{'FaceAlpha'},{TransperancyValue})
            Title (TitleString);
            hold on
            plot3(ValuesOfPeak(1),ValuesOfPeak(2), PeakValue , ValuesOfPeak(1),ValuesOfPeak(2), PeakValue, 'g:o','MarkerSize',6, 'LineWidth', 1.5);
            hold off
        end    
 
        if length(FilePrefix)>1,
            try
                FigurePrintFileName = [FilePrefix TitleString];
                eval (['print -dpdf ' FigurePrintFileName])
            catch
                disp 'could not print to file - check for file locks or for illegal characters in title/filename'
            end
        end
        
        GraphVersion = 0;
        
    end % of while GraphVersion
end  % GraphNum = 1:NumberOfGraphs
