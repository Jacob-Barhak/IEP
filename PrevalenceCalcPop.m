function [ Prevalence ] = PrevalenceCalcPop( Dim, Range, Pop )
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
% function [ Prevalence ] = PrevalenceCalcPop( Dim, Range, Pop )
% The function calculates the Prevalence of a table cell defined by the
% dimension cell array Dim and the associated Range within the 
% population Pop. Range is defined as a matrix of size #Dim x 2 were each 
% row has the low and high range values for dimension Dim
 
 
% For easy reading use H and L to represent the high and low range values
% in the Range Matrix.
L = Range(:,1);
H = Range(:,2);
 
 
% Check if population is defined by data or by distribution
if ~isempty(Pop.Data),
    % Population is defined by data
 
    % Init individual compatibility vector
    Compatible=ones(size(Pop.Data,1),1);
    ValidData=Compatible;
    
    % Loop through the dimensions. 
    for i=1:length(Dim),
        
        if strcmp(Dim(i),'Dummy')
            % if a dummy dimension appears, ignore it. This simulates a
            % Dummy dimension in each population that conforms to any data.
            continue;
        end
        
        % Find the column index corresponding to the given input dimensions
        Cells = {Pop.Columns.Dim};
        DimInd = find(strcmp(Cells,Dim(i)));
 
        % Consistency check
        if isempty (DimInd),
           error ('Dimension %s does not exist in population %s', char(Dim(i)), Pop.Name) 
        end
   
        % Handle the case where the low bound is NaN. This may happen in
        % case of a first cell in an indexed dimension. In this case 
        % replace the low bound with minus infinity
        if isnan(L(i)),
           L(i)=-inf; 
        end
        
        % Perform the comparison for this dimension and save the result in
        % a logical vector 
        DimCompatible = (Pop.Data(:,DimInd)>L(i)) & (Pop.Data(:,DimInd) <=H(i));
        
        % Construct a similar logical vector that signifies valid data,
        % i.e. is there was no NaN in a specific dimension. Note that
        % previous operations with NaN would have returned a zero
        DimValidData = ~isnan(Pop.Data(:,DimInd));
        
        % Update individual compatibility and data validness
        Compatible = Compatible & DimCompatible;
        ValidData = ValidData & DimValidData;
        
    end
 
    % The prevalence is the average value of the compatibility vector.
    % Note that the individual count is the sum of the compatibility vector
    % and when divided by the size of the population represented by the
    % vector size, the prevalence is calculated. This is possible since
    % the vector is of Booleans. Note that only valid data entries will be
    % considered.
    Prevalence = mean(Compatible(ValidData));
    
else
    % Population is defined by distribution
 
    % Initialize the prevalence multiplication temporary variable to 1
    MultTemp = 1;
    for i=1:length(Dim),
 
        if strcmp(Dim(i),'Dummy')
            % if a dummy dimension appears, ignore it. This simulates a
            % Dummy dimension in each population that conforms to any data.
            continue;
        end
        
        % Find the column index corresponding to the given input dimensions
        Cells = {Pop.Columns.Dim};
        DimInd = find(strcmp(Cells,Dim(i)));
 
        % Consistency check
        if isempty (DimInd)
           error ('Dimension %s does not exist in population %s', Dim, Pop.Name) 
        end
 
        % First extract the distribution information by extracting its name
        % in one string. Another string will keep the rest of the 
        % information in the parenthesis and keeping the rest of the 
        % variables together along with the separating comma characters 
        % as this parameters string will be reused later.
        DistrInfo = strread(Pop.Columns(DimInd).Distribution,'%s','delimiter','()');
        DistrType = DistrInfo(1);
        DistrParams = DistrInfo(2);
 
        % The distribution parameters should be passed to the function as
        % defined in the population vector. For this an expression is
        % created that calls MCDF with the distribution type and
        % parameters. Although this requires constructing a command and its
        % evaluation to support the variable number of parameters, this
        % expression is essentially the following line:
        % MyCDF(DistrType,x,DistrParams) 
        
        % Prepare command strings
        ExpressionL = sprintf('MyCDF(''%s'',%g,%s);', char(DistrType) , L(i) , char(DistrParams) );
        ExpressionH = sprintf('MyCDF(''%s'',%g,%s);', char(DistrType) , H(i) , char(DistrParams) );
 
        % Evaluate strings to expressions
        ComponentL=eval(ExpressionL);
        ComponentH=eval(ExpressionH);
        
        % Use expressions in the prevalence formula
        MultTemp = MultTemp * (ComponentH - ComponentL);
        % Note that low bounds that are NaN are handled internally by MyCDF
    end
    
    Prevalence = MultTemp;
    
end
    
