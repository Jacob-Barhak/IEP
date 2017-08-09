function [ MeanValue ] = MeanCalcPop( Dim, Range, Pop )
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
% [ MeanValue ] = MeanCalcPop( Dim, Range, Pop )
% The function calculates MeanValue which is the mean of the dimension Dim
% within a given value Range in the population Pop.
% If Range is [-inf inf] then the function returns the mean from the
% distribution and the mean of the relevant dimension from the population
% set. Regardless of the range, when population data is used, individuals
% with no data in the given dimension are ignored.
 
% For easy reading use H and L to represent the high and low range values
L = Range(1);
H = Range(2);
 
% Find the column index corresponding to the given input dimensions
 
Cells = {Pop.Columns.Dim};
DimInd = find(strcmp(Cells,Dim));
 
% Consistency check
if isempty (DimInd)
   error ('Dimension %s does not exist in population %s', Dim, Pop.Name) 
end
 
% Check if population is defined by data or by distribution
if ~isempty(Pop.Data),
    % Population is defined by data
 
    % Handle the case where the low bound is NaN. This may happen in
    % case of a first cell in an indexed dimension. In this case 
    % replace the low bound with minus infinity
    if isnan(L),
       L = -inf; 
    end
 
    
    % Find all individuals that correspond to the range in the input 
    % dimension and have data defined for the given dimension.
 
    % Construct a logical vector that will index all individuals with data
    % in the given dimension. 
    ValidIndividuals = ~isnan(Pop.Data(:,DimInd));
    % Then extract the valid data - the dimension column in the table
    ValidData=Pop.Data(ValidIndividuals,DimInd);
    
    % Now see if the valid data is in the range
    DataInd = find( ValidData > L & ValidData <= H) ;
 
    % In case no individual found in this range
    if isempty(DataInd)
        error ('No individuals in population %s correspond to the given range [%g %g] in dimension %s. Consider exchanging the population with another', Pop.Name,L,H,Dim);
    end
    
    % Calculate the mean
    MeanValue = mean (ValidData(DataInd));
    
else
    % Population is defined by distribution
        
    % First extract the distribution information
    DistrInfo = strread(Pop.Columns(DimInd).Distribution,'%s','delimiter',',()');
    DistrType = DistrInfo(1);
 
    if strcmpi(DistrType,'Ber') || strcmpi(DistrType,'Bernouli'),
        %The format is Bernouli(p) therefore the second token is p
        p=sscanf(cell2mat(DistrInfo(2)),'%g');
        if all(Range == [-inf inf]),
            MeanValue = p;
        else
            MeanValue = H;
        end
        % Just make a consistency check
        if MeanValue>1 || MeanValue<0 ,
            error('The mean value calculated for Bernouli is not in the range - Input range was [%g %g]', L,H);
        end
        
    elseif strcmpi(DistrType,'Bino') || strcmpi(DistrType,'Binomial'),
 
        % Handle the case where the low bound is NaN. This may happen in
        % case of a first cell in an indexed dimension. In this case 
        % replace the low bound with the low bound of the Binomial
        % distribution support, i.e. L=0
        if isnan(L),
           L = 0; 
        end
 
        % The format is Binomial(n,p) 
        % Therefore the second token is N and the third is p
        n=sscanf(cell2mat(DistrInfo(2)),'%g');  
        p=sscanf(cell2mat(DistrInfo(3)),'%g');
        if all(Range == [-inf inf]),
            MeanValue = n*p;
        else
 
            % This code uses the initial expression defined
            Nom=0;
            Denom=0;
            for i=L:H,
                Element = nchoosek(n,i)*p^i*(1-p)^(n-i);
                Nom=Nom+i*Element;
                Denom=Denom+Element;
            end
 
            MeanValue = Nom/Denom;
        end
        
    elseif strcmpi(DistrType,'Geo') || strcmpi(DistrType,'Geometric'),
        % The format is Geometric(p)
        % Therefore the second token is p
        p=sscanf(cell2mat(DistrInfo(2)),'%g');
        
        % Handle the case where the low bound is NaN. This may happen in
        % case of a first cell in an indexed dimension. In this case 
        % replace the low bound with the low bound of the Geometric 
        % distribution support, as it is defined here i.e. L=1
        if isnan(L),
           L = 1; 
        end
 
        % If the range is defined as infinite
        if all(Range == [-inf inf]),
            MeanValue = 1/p;
        % Otherwise Check if the range is valid
        elseif (L<1 || H<1),
            % Return NaN if L and H are out of the support range of the
            % geometric distribution.
            error('The range [%g,%g] that was provided to the geometric distribution is out of the support range',L,H);
        else % otherwise use the following equation
            MeanValue = ( -(H*p+1)*(1-p)^H + (L*p+1-p)*(1-p)^(L-1) )/( p*((1-p)^(L-1) - (1-p)^H ));
        end
 
    elseif strcmpi(DistrType,'unif') || strcmpi(DistrType,'Uniform'),
        % The format is Uniform(a,b)
        % Therefore the second token is a and the third is b
        a=sscanf(cell2mat(DistrInfo(2)),'%g');
        b=sscanf(cell2mat(DistrInfo(3)),'%g');
 
        % Handle the case where the low bound is NaN. This may happen in
        % case of a first cell in an indexed dimension. In this case 
        % replace the low bound with the low bound of the Uniform 
        % distribution support, as it is defined by the user
        if isnan(L),
           L = a; 
        end
 
        MeanValue = (max(L,a)+min(H,b))/2;
        
        % Consistency check
        if MeanValue>H || MeanValue<L,
           error('The Range [%g,%g] is defined outside the distribution boundaries [%g,%g]',L,H,a,b);
        end
 
    elseif strcmpi(DistrType,'norm') || strcmpi(DistrType,'Normal'),
        % The format is norm(Mean,STD)
        % Therefore the second token is Mean and the third is STD
        Mean=sscanf(cell2mat(DistrInfo(2)),'%g');
        STD=sscanf(cell2mat(DistrInfo(3)),'%g');
        % Convert to normalized variables
        Zh=(H-Mean)/STD;
        Zl=(L-Mean)/STD;
        NormalizedMean = (pdf('norm',Zl,0,1)-pdf('norm',Zh,0,1)) / (MyCDF('norm',Zh,0,1) - MyCDF('norm',Zl,0,1));
        MeanValue =  NormalizedMean * STD + Mean;
        
    else
        % Unsupported distribution name
        error('This distribution function is not supported by the prototype');
    end
end
