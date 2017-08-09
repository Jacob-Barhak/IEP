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
% CalcLikelihood.m
% This script implements part of the Lemonade Method.
% This script builds the log Likelihood expression.
% This code requires that all the proper definitions of the model, studies,
% and associated population be defined already before the script starts.
% In a sense, this script acts as a function that outputs the log
% likelihood L. Since it is not encapsulated as a function, use this code 
% cautiously. It should be called once just after defining the input data.
% Note that this script builds the Log likelihood and does not optimize it.
% The optimization script OptimizeLikelihood.m should be called afterwards.
%
% INPUT:
%
% Constants - All the Constant names as stings in a cell vector.
% ConstantValues - A vector defining the values of the constants  
%                  defined in Constants.
% Covariates - All the Covariate names as stings in a cell vector.
% Coefficients - All the coefficient names as stings in a cell vector
% Params - The union of Covariates Coefficients and Constants. Note that
%          this should be defined before running the script and there
%          should be a set of symbolic variables with the same names
%          defined in Params. The script InitParams.m generates these
%          properly.
% Pop - A cell vector of structures defining the population sets. The index
%       in the vector is the population ID and each structure contains the 
%       following fields:
%       'Name' - A string containing the name of the population set
%       'Columns' - A cell vector of structures representing the columns to
%                   be used in this population set with the following
%                   fields:
%                   'Dim' - A string containing the dimension name. This
%                           should be defined as a Covariate.
%                   'Distribution' - If the population has data, then this 
%                                    field is ignored. Otherwise, it should
%                                    contain the distribution that defines
%                                    this column in the population set. For
%                                    more information about these
%                                    distributions see MyCDF.m .
%                                    Distributions will be strings of one 
%                                    of the following forms:
%                                    Ber(p) - for Bernoulli
%                                    Normal(mean,STD) - for Normal Gaussian 
%                                    bino(n,p) - for Binomial distribution
%                                    Geo(p) - for geometric  distribution
%                                    Uniform(a,b) - for Uniform 
%                   'OrderInSet' - Required for population sets with data
%                                  and defines the column number in the 
%                                  data matrix.
%       'Data' - This field contains a matrix with the population data,
%                Each row in the matrix corresponds to data describing an 
%                individual from this population. Each column corresponds
%                to a column defined in columns and located according to 
%                the OrderInSet field.
% States - A cell vector of structures defining the states and
%          sub-processes used in the study. The index in the vector is the
%          Study ID and each structure contains the following fields:
%          'Name' - A string containing the name of the state
%          'IsEvent' - A Boolean that indicates if this is an event state.
%                      An event state is instantaneous and is handled
%                      differently during optimization.
%          'IsSplit' – Future implementation. Not currently supported. Use 
%                      Zero for a value for this field.
%          'JoinerOfSplit' - Future implementation. Not currently supported. 
%                            Use Zero for a value for this field.
%          'PoolingData' - A two row matrix containing information about
%                          states nested within this state. This structure
%                          takes care of both pooled states and
%                          sub-processes. In a regular state that does not
%                          pool any state an empty matrix, i.e. [], is used.
%                          Otherwise the first row of the matrix will hold
%                          the State IDs of all states nested in the
%                          Sub-Process/Pooling State. The second row will
%                          will contain the prevalence values associated
%                          with the pooled states in the first row for a
%                          pooling states. The second row will be all zeros
%                          for a sub-process.
% Studies - A cell vector of structures defining the model and the studies.  
%           The index in the vector is the Study ID and each structure 
%           contains the following fields:
%          'Name' - A string containing the name of the Study
%          'StudyLength' -  Defines the study length in time units. For a
%                           model, this number will be 0. Note that there
%                           should be only one model in the list and this
%                           zero identifies it.
%          'MainProcess' - The State ID corresponding to the sub-process
%                          that describes the models main process. For a
%                          study, this number can be zero. For a model, it
%                          should be a sub-process ID that contains all the
%                          states in the model. Currently, nesting of
%                          sub-process states in other sub-processes is not
%                          supported as this is a implementation issue. 
%                          Nevertheless, the main sub-process still should
%                          be defined for a model.
%          'PopID' - The population ID of the population to be used with
%                    this study during optimization. For a model, this can
%                    be empty, i.e. []. For a study, this population set
%                    should provide sufficient information to accommodate
%                    both study and model. Basically, this means that any
%                    covariate defined by the study and by the model should
%                    be defined properly in the population set.
% Transitions - A cell vector of structures defining the transitions 
%               between the states in all models and studies. Each 
%               structure contains the following fields:
%               'StudyID' - The Study ID indicating the study from the
%                           Studies vector that represents the Study/Model 
%                           to which this transition belongs.
%               'From' - The State ID indicating the state in the States
%                        vector from which this transition starts. Note
%                        that this state should belong to the model/study
%                        defined for this transition in StudyID.
%               'To' - The State ID indicating the state in the States
%                      vector to which this transition leads, or zero
%                      in case this transition holds initial population
%                      counts for a study. Note that this state should
%                      belong to the model/study defined for this
%                      transition in 'StudyID'.
%               'Probability' - Holds probability related information
%                               associated with this transition. The
%                               information stored in the field will change
%                               depending on the transition type: 
%                               1) For a study this field will hold the
%                                  start population in one of the following
%                                  forms:
%                                   1a) Initial population counts - If the 
%                                       study is based on categorical data
%                                       and the 'To' field is zero, it 
%                                       contains initial population count 
%                                       of individuals starting at the 
%                                       'From' state in that study. This 
%                                       data can be or provided as a number
%                                       in Table Structure (see below).
%                                   1b) If the study based on categorical 
%                                       data and the 'To' state is non zero
%                                       it contains incident counts of
%                                       Individuals starting from the 
%                                       'From' state and ending in 'To'
%                                       state during the study length.
%                                       This data can be provided as a
%                                       number or in Table Structure (see
%                                       Below). Note that the table can
%                                       contain the Time dimension, in
%                                       which case the incident counts will
%                                       be associated with the time of the
%                                       table cell. 
%                                   1c) If the study is based on Regression
%                                       data then this field should be
%                                       empty.
%                               2) For a Model this field will hold the
%                                  probability of progression from the
%                                  'From' State to the 'To' state. This
%                                  information can be provided as a
%                                  number or a coefficient / expression 
%                                  containing a coefficient / constants.
%                                  The expression should use the symbolic
%                                  variables rather than be written as a
%                                  string as the system uses symbolic math
%                                  to work with the expression. Note that
%                                  Expressions using covariates are a
%                                  future implementation issue. Also note 
%                                  that the field may hold a table 
%                                  structure containing simple expressions
%                                  within the table cells. See below more
%                                  information on Table Structures.
%               'RegCoefVec' - Contains the regression coefficient vector
%                              for regression studies. Must be empty for
%                              Model transitions and categorical study
%                              transitions. The vector will contain a set
%                              of numbers reported in the literature as
%                              regression coefficients and will be
%                              associated with the parameters provided in
%                              the field 'RegParamVec'.
%               'RegParamVec' - Contains the regression parameters vector
%                               for regression studies. Must be empty
%                               for Model transitions and categorical study
%                               transitions. The vector will contain a set
%                               of members reported in the literature 
%                               as regression parameters and will be
%                               associated with the regression coefficients
%                               provided in the field 'RegCoefVec'. Each
%                               member is either a Covariate, an
%                               expression involving a Covariate or a
%                               Number/Constant in case a Bias value is
%                               required in the regression formula.
%               'CovarianceMatrix' - Contains the regression covariance
%                                    matrix for regression studies. Must be
%                                    empty for Model transitions and
%                                    categorical study transitions. The
%                                    matrix members are numbers reported in
%                                    the literature and will be associated
%                                    with the parameters provided in the
%                                    field 'RegCoefVec'.
%               'RegFunctionType' - Contains a constant defining the
%                                   regression function type. Must be
%                                   empty for Model transitions and
%                                   categorical study transitions.
%                                   Currently, the following Regression
%                                   function Types are supported and are
%                                   represented by the following Integers:
%                                   1 - Exponential form, i.e. 
%                                       exp(-RegParamVec*RegCoefVec')
%                                   2 - Linear form, i.e. 
%                                       RegParamVec*RegCoefVec'
%                                   3 - UKPDS form, i.e.
%                                       1-exp(-(1-d^t)/(1-d)*prod(a.^b))
%                                       where:
%                                       a = RegCoefVec(1:(end-1))
%                                       b = RegParamVec(1:(end-1))
%                                       d = RegCoefVec(end)
%                                       t = RegParamVec(end) and should be
%                                           the same as the Study Length
%                                       Additional details may be found in:
%                                       Stevens, R. Kothari, V. Adler, A.
%                                       Stratton I. (2001). The UKPDS risk
%                                       engine: A model for the risk of
%                                       coronary heart disease in type II
%                                       diabetes UKPDS 56. Clinical
%                                       Science, 101, 671-679.
%                                   11 - One Minus Exponential form, i.e. 
%                                        1-exp(-RegParamVec*RegCoefVec')
% 
% Table Structures: Some transitions may require information in table form
%                   that clusters data according to values of one or more
%                   Covariates. The system supports representation of data
%                   in multiple-dimensions. The Information should be
%                   supplied to the system in a structure having the 
%                   following fields:
%                   'Data' - This field will hold a multidimensional array 
%                            that contains numbers or symbolic expressions.
%                            Matlab conventions are used to access the
%                            array. Note that tables with one dimension,
%                            i.e. vectors, should be column vectors. The
%                            data should be ordered such that the index 
%                            vector to access the table member should
%                            correspond to the order of dimensions defined
%                            in the field 'Dim'
%                   'Dim' - This field contains a vector defining the
%                           dimensions and the dimension ranges associated
%                           with the dimension. For each dimension d in
%                           'Data' of size n, this vector will contain n+2
%                           members. The first member will be the name of
%                           the dimension. It must be previously defined
%                           as a parameter and its symbolic variable should
%                           exist. Also, in the case of studies, it is
%                           possible to use the symbolic variable Time to
%                           define the dimension name. The next n+1
%                           variables will define the Range of this
%                           dimension and how it splits space. Each range
%                           can describe an indexed dimension or a
%                           continuous range. Here are details:
%                           Indexed range - The first member of the range
%                                           will be NaN and the next n
%                                           values in the range will define
%                                           the index values. The table can
%                                           be accessed by these discrete
%                                           index values by making the
%                                           check that x = Range(I+1),
%                                           where x will be the value
%                                           supplied to access the table
%                                           and I will be the dimension to
%                                           be accessed in the array
%                                           defined in 'Data'.
%                           Continuous Range - All n+1 range values will be
%                                              used to define range 
%                                              boundary separators. so that
%                                              the index to access the
%                                              table can be deduced by 
%                                              Range(I) < x <= Range(I+1)
%                                              where x will be the value
%                                              supplied to access the table
%                                              and I will be the dimension 
%                                              to be accessed in the array
%                                              defined in 'Data'.
% 
% OUTPUT:
%
% L - The log likelihood expression constructed by combining the model and
%     all the studies and the associated population sets. This symbolic
%     expression can then be optimized using OptimizeLikelihood.m to
%     extract the optimal values for the coefficients used in it.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Record the time when starting the program.
tic
 
% Define a constant 
% Perturbation size for regression gradient is set here as a constant
h = 1e-8;
 
% Define additional debug prints in the following line
% The debug Vector determines what debug data to print
DebugVec = zeros(1:10);
DebugVec(1) = 1; % Print output information for each study
DebugVec(2) = 1; % Print output information while constructing Likelihood
 
% Handle Future implementation. The user should add the following line 
% FUTURE_IMPLEMENTATION = false
% before the running this script if the use of future implementation 
% options are wanted.
if ~ismember('FUTURE_IMPLEMENTATION',who),
    FUTURE_IMPLEMENTATION = false;
end
 
%%%% Verify Study/Model/Population correlation %%%%
% This part of the script verifies that the input data is valid 
 
SufficientInfo = true;
 
% Find the model ID in the studies collection
Vec = cell2mat({Studies.StudyLength});
ModelID = find(Vec==0);
 
%Check that there is only one model
if length(ModelID)>1,
    error ('More than one model is defined, make sure only one study has StudyLength=0')
elseif length(ModelID)<1,
    error ('No Model is defined, make sure there is one study with StudyLength=0')
end
 
% Find list of indices related to the model in the Transition structure
Vec = cell2mat({Transitions.StudyID});
ModelTransInd = find(Vec==ModelID);
 
for StudyID = 1:length(Studies), % Study loop
    % reset dimensions
    D=[];
    
    % Find list of indices related to the study in the Transition structure
    Vec = cell2mat({Transitions.StudyID});
    StudyInd = find(Vec==StudyID);
    
    % Add model transition indices to the study indices
    Ind = [StudyInd ModelTransInd];
    
    % Traverse all the transitions related to the study
    for t = Ind,
        % If the transition is a table
        if isstruct(Transitions(t).Probability),
            % Extract the dimensions from the table info
            [Data, Dim, Sizes, Ranges] = TableInfo( Transitions(t).Probability );
        % else if the transition is an expression in a model
        elseif isobject(Transitions(t).Probability),
            % Find all parameters
            Params = symvar(char(Transitions(t).Probability));
            % Remove coeficients from this list
            WithoutCoefficients = setdiff(Params,Coefficients);
            % Remove Constants from this list
            Dim = setdiff(WithoutCoefficients,Constants);
        % else if the transition has a single probability in it
        elseif isnumeric(Transitions(t).Probability) && ~isempty(Transitions(t).Probability),
            % No dimensions in this case
            Dim = [];
        % else if the transition has regression data in it
        elseif isempty(Transitions(t).Probability) && ~isempty(Transitions(t).RegParamVec)  && ~isempty(Transitions(t).RegCoefVec) && ~isempty(Transitions(t).CovarianceMatrix) && ~isempty(Transitions(t).RegFunctionType),
            Dim = symvar(char(Transitions(t).RegParamVec)) ;
        elseif Transitions(t).From ~= 0, 
            % Raise an error as the data is unrecognized
            % Note that symbolic coefficients are such
            error('Invalid data defined in transition %s from state %s to state %s ',Studies(Transitions(t).StudyID).Name,States(Transitions(t).From).Name,States(Transitions(t).To).Name);
        end
        
        % Collect the new dimension in the dimension list D
        D=union(D,Dim); 
        
    end  % End Transition Loop
 
    % Remove time from the list as it will not be in the population
    D = setdiff(D,'Time');
    
    % Find list of indices related to the study population in the relevant
    % population column list
    PopID = Studies(StudyID).PopID;
    
    % If a population is defined in the study
    if ~isempty(PopID) && PopID ~= 0,
        % Extract all columns relevant to this population
        Columns = {Pop(PopID).Columns.Dim};
        MissingDim=setdiff(D,Columns);
        if ~isempty(MissingDim);
            MsgStr = sprintf ('Error : Population set "%s" does not contain the dimensions %s required by study %s and by the model', Pop(PopID).Name, char(MissingDim) ,Studies(StudyID).Name);
            disp(MsgStr);
            SufficientInfo = false;
        end
    end 
end % End Study Loop
 
if ~SufficientInfo,
   error ('There is no sufficient population information to estimate the project'); 
end
 
%%%%% End of Verify Study/Model/Population correlation %%%%
 
 
%%%% Construct Log Likelihood %%%%%
 
% Init general Log Likelihood term to 0
L=sym(0);
FunctionListStudyLogLikelihood = {};
FunctionListCellLogLikelihood = {};
FunctionListTransLogLikelihood = {};
FunctionListStudyCellLogLikelihood = {};

% Mark all studies as unused
StudiesUsed=zeros(1,length(Studies));
 
%%%%%% Start main loop %%%%%
 
% Find all sub-processes in the model
[PooledStates, NestingLevel, IsSubProcess] = FindSubProcessStates(States, Studies(ModelID).MainProcess);
SubProcessIndices = PooledStates(logical(~IsSubProcess));
 
% The main subprocess loop
for SubProcess = SubProcessIndices,  
    % Find all states in the sub-process. Ignore study pooling states and
    % treat only pooling states with prevalence zero
    StatesInSubProcessUnsorted = States(SubProcess).PoolingData(1,~States(SubProcess).PoolingData(2,:));
    
    % Sort this vector for convenience in reading and handling the code
    StatesInSubProcess = sort(StatesInSubProcessUnsorted);
    % Find the indices of all event states
    EventStateIDs = find([States.IsEvent] == 1);
    % Intersect the lists to find the event states in the model sub-process
    EventStateIDsInSubProcess = intersect (EventStateIDs, StatesInSubProcess);
    
    %%%%% Construct the R, E matrices %%%%%
    
    n=length(StatesInSubProcess);
    % Init the non event probability matrix to zeros and make it symbolic
    R=sym(zeros(n,n));
    % Init the event probability matrix to unity and make it symbolic
    E=sym(eye(n,n));
    
    % Create an index map to help find an index of a state in a Sub-Process
    ModelMatIndMap = zeros(1,length(States));
    ModelMatIndMap(StatesInSubProcess) = 1:length(StatesInSubProcess);
    
    % Find all transitions associated with the states in the Sub Process.
    % First, find transitions into Sub Process states
    TransIntoSubProcInd=find(ismember([Transitions.To],StatesInSubProcess));
    ModelTransIntoSubProcInd=intersect(TransIntoSubProcInd,ModelTransInd);
    % Then find transitions out of Sub Process states
    TransFromSubProcInd=find(ismember([Transitions.From],StatesInSubProcess));
    ModelTransFromSubProcInd=intersect(TransFromSubProcInd,ModelTransInd);
    % Combine results to find all the transitions in the Sub-Process
    SubProcTransInd = intersect(ModelTransFromSubProcInd, ModelTransIntoSubProcInd);
    
    % Just a consistency check: make sure that transitions that lead into
    % a sub-process and out of a sub-process lead from a splitter state or 
    % into a joiner state.
    PartSubProcTransInd = setdiff(union(ModelTransFromSubProcInd, ModelTransIntoSubProcInd),SubProcTransInd);
    % If empty, no transition is split inside and outside the sub-process
    if (~isempty(PartSubProcTransInd)),
        if any(~ismember( [Transitions(PartSubProcTransInd).From] , StatesInSubProcess ) && ~States([Transitions(PartSubProcTransInd).From]).IsSplit),
            error ('Transitions %s are transitions originating from outside this sub-process, yet their originating states are not a splitter state', num2str(PartSubProcTransInd));
        end
        if any(~ismember( [Transitions(PartSubProcTransInd).To] , StatesInSubProcess ) && States([Transitions(PartSubProcTransInd).To]).JoinerOfSplit ~= 0),
            error ('Transitions %s are transitions originating from outside this sub-process, yet their originating states are not a splitter state', num2str(PartSubProcTransInd));
        end
    end
    % The consistency check ends here
    
    % For each relevant transition create a symbolic variable and populate
    % the Matrix with it
    for t = SubProcTransInd,
        % Raise an error if the transition is a self loop
        if Transitions(t).From == Transitions(t).To
           error('Transition %i is a self loop',t);
        end
        
        % Create a symbolic variable of the form T#_#_# where the first 
        % number corresponds to the study id, the second to the start
        % state and the third to the end state. 
 
        VarStr = sprintf ('T%i_%i_%i',Transitions(t).StudyID, Transitions(t).From, Transitions(t).To);
        CommandStr = sprintf ('syms %s',VarStr);
        eval(CommandStr);
 
       
        % Determine the index position in the P matrix for this symbol by
        % extracting the index from the StatesInSubProcess vector
        i=ModelMatIndMap(Transitions(t).From);
        j=ModelMatIndMap(Transitions(t).To);
 
        % check if the transition is from an event state
        if ~ismember(Transitions(t).From, EventStateIDsInSubProcess)
            % not an event state - update the non event probability matrix
            
            % Update the matrix i,j element to the transition name
            CommandStr = sprintf ('R(i,j) = %s;',VarStr);
            eval(CommandStr);
        else
            % An event state - update the event probability matrix
            % First, reset the diagonal
            E(i,i)=0;
            % Update the matrix i,j element to the transition name
            CommandStr = sprintf ('E(i,j) = %s;',VarStr);
            eval(CommandStr);
        end
        
    end % Transition loop
    
    % Update the diagonal of the non event probability matrix to reflect 
    % the fact that this matrix represents a Markov model, meaning that the
    % diagonal completes the probability of each row to 1
    for i=1:n,
        R(i,i)=1-sum(R(i,setdiff(1:n,i)));
    end
 
    %%%%%%% end of R,E construction %%%%%%%%
    
    %%%%%% Modify P for each study %%%%%%%
    
    % Loop through all studies
    
    % Find the ID of all studies - without the model
    AllStudies = find([Studies.StudyLength]);
 
    % Loop through all studies
    for StudyID=AllStudies,
        
        % Just a progress notification print
        disp '--------------------------------------'
        disp (sprintf('Processing Study #%i - "%s" ',StudyID,Studies(StudyID).Name));
        disp ''
 
        
        % Extract all the states participating in the study by traversing
        % transitions and collecting states and pooled state through these
       
        % Init lists to nothing
        RegularStates=[];
        PooledStates=[];
        Prevalences=[];
        PoolingStateInfo=[];
        PoolingStatesUsed=[];
 
        
        % Extract transitions associated with the study
        StudyTransInd = find([Transitions.StudyID]== StudyID);
        
        % Loop through these transitions and collect states
        % While looping perform consistency checks to avoid invalid
        % definitions of pooling states
        % PooledStates will hold the states, while Prevalence will
        % hold the prevalence values associated with these states. 
        % PoolingState will hold the state that pooled the states
 
        for t = StudyTransInd,
            % Repeat twice for From and for To state
            for  s = [Transitions(t).From Transitions(t).To],
                % Ignore Null and states that were previously processed 
                if s ~= 0 && ~ismember(s,[PooledStates,PoolingStatesUsed,RegularStates]),
                    % If this is a pooling state
                    if ~isempty(States(s).PoolingData),
                        PoolingData = States(s).PoolingData(:,States(s).PoolingData(2,:)>0);
                        % A consistency check if the pooled states are 
                        % valid for this model
                        if any(ismember(PoolingData(1,:),[PooledStates,PoolingStatesUsed,RegularStates])),
                            error('The pooled States list for study %s is invalid',Studies(StudyID).Name);
                        end
                        % Add the pooled states to the list
                        PooledStates = [PooledStates, PoolingData(1,:)];
                        Prevalences = [Prevalences, PoolingData(2,:)];
                        PoolingStateInfo = [PoolingStateInfo, ones(1,size(PoolingData,2))*s];
                        PoolingStatesUsed = [PoolingStatesUsed s];
                    else % No pulled data
                        % Collect the From state in a set with no pooled states
                        RegularStates=[RegularStates s];
                    end
                elseif any(ismember(PooledStates,[RegularStates PoolingStatesUsed])) || any(ismember(RegularStates, [PooledStates PoolingStatesUsed])) || any(ismember(PoolingStatesUsed, [PooledStates RegularStates])),
                    error('The pooled States list for study %s is invalid',Studies(StudyID).Name);
                end
            end
        end
 
        % Sort the states for ease of debugging
        [PooledStates SortInd] = sort(PooledStates);
        Prevalences = Prevalences(SortInd);
        PoolingStateInfo = PoolingStateInfo(SortInd);
        PoolingStatesUsed = sort(PoolingStatesUsed);
        RegularStates = sort(RegularStates);
        
        % This is the state list without pooled states and with all the
        % pooling states in it. Note that it will be used to provide index
        % for the matrices associated with the study later on. 
        StudyStates = [RegularStates PoolingStatesUsed];
        
        % Check if the study is fully contained in the sub-process
        if ~all(ismember([RegularStates PooledStates],StatesInSubProcess)),
            continue; % to the next study in the loop
        else
            % If the study has not been used, mark as used
            if ~StudiesUsed(StudyID), 
                StudiesUsed(StudyID)=1; 
            else % Otherwise raise an assertion error
                error ('Assertion error - Study %s, has already been used during the loop',Studies(StudyID).Name);
            end
        end
 
        
            
        % Check if all transitions are event states
        % First extract all outcome states that are not zero
        OucomeStatesWithNULL = [Transitions(StudyTransInd).To];
        OucomeStates = OucomeStatesWithNULL(OucomeStatesWithNULL~=0);
        IsOutcomeEventStates = [States(OucomeStates).IsEvent];
 
        A = E;
 
        if all(~IsOutcomeEventStates),
            % If all outcomes are non event states 
            % use the Event Including Model Matrix
            
            % Currently do nothing as A is initialized before 
            
        elseif length(IsOutcomeEventStates)==1 && IsOutcomeEventStates == 1,
            % If there is only one outcome state, 
            % sse the Event Ignoring Model Matrix
 
            k=ModelMatIndMap(OucomeStates);
            
            % Create the modified event matrix
            A(k,:) = sym(zeros(1,n));
            A(k,k) = sym (1);
         
        else
            % If the study uses more than one event state as an outcome 
            % or combines non event and event states, raise an Error
            error ('Study #%i - "%s" uses more than one event state as an outcome or combines both event and non event states as study outcomes. This is not currently supported by the system. ', StudyID, Studies(StudyID).Name);
        end
        
 
        % Calculate the consecutive event probability matrix.
        % Extract the number of event states in the sub process
        EventStateNum = length(EventStateIDsInSubProcess);
        MultA = A^(EventStateNum);
 
        % Consistency check to establish there are no event loops
        SelfLoops = (diag(MultA) ~= 0) & (diag(A) == 0);
        if any(SelfLoops),
            EventStateID = find (ismember(SelfLoops,ModelMatIndMap));
            error('Consecutive Event Transitions lead to a self loop in event states #%s',mat2str(EventStateIDs));
        end
 
        % Calculate the Event ignoring Probability matrix P by adding the 
        % Regular and the event probability matrices. Note that the 
        % diagonal is incorrect and should be modified later on
        P = (R)*MultA;
        
        
        %%%% Construct the K and K1 matrices
        % Init variables
        m = length(StudyStates);
        K = zeros(m,n);
        K1 = zeros(n,m);
 
        % Create an index map for this study
        StudyMatIndMap = zeros(1,length(States));
        StudyMatIndMap(StudyStates) = 1:length(StudyStates);
        
        % Loop through all sub-process states,
        for i = 1:length(StatesInSubProcess),
            % For each sub-process state 
            
            % Compare the state to the list of study states. The result
            % vector is 1 in the position associated with the study index j
            j = StudyMatIndMap(StatesInSubProcess(i));
            % Check if the state is in the study regular states
            if j~=0,
                K(j,i)=sym(1);
                K1(i,j)=sym(1);
            else % Otherwise it is a pooled state
                % Find this state in the list that includes pooled states
                PooledID = find(StatesInSubProcess(i) == PooledStates);
                % If the process state is not in the study - continue 
                if ~isempty(PooledID), 
                    % Extract the Prevalence and the pooling state
                    PoolingState = PoolingStateInfo(PooledID);
                    StatePrevalence = Prevalences(PooledID);
                    % Extract the index of the pooling state
                    z = StudyMatIndMap(PoolingState);
                    % Update the K matrices
                    K(z,i)=sym(StatePrevalence);
                    K1(i,z)=sym(1);
                end
            end
 
            
            %%%%% Handle Sink states %%%%%
            % Traverse the study transitions,
            for t = StudyTransInd,
                % If the transition outcome is state i
                if (Transitions(t).To == StatesInSubProcess(i)),
                    % Reset row i to 1 in the diagonal and zero elsewhere
                    P(i,:)=sym(0);
                    P(i,i)=sym(1);
                end
            end % Transition loop
        end % States in sub process 
        
        % A consistency check: verify that the sum of prevalences in K is 1
        VerifyVec = sum(K,2);
        NotOnes = find(VerifyVec ~= 1);
        
        if ~isempty(NotOnes),
            StateNames = {States(StudyStates(NotOnes)).Name};
            StateNamesFormat = StateNames{1};
            for i=2:length(StateNames),
                StateNamesFormat = [StateNamesFormat, ', ', StateNames{i}];
            end
            error('The prevalence values in study %s do not sum to 1 in the following study states %s', Studies(StudyID).Name, StateNamesFormat);
        end
        
        %%%%%% End of constructing K and K1 Matrices
        
        %%%%%% Take care of retrospective studies
        
        % Traverse study transitions
        for t = StudyTransInd,
            % Check if this is a retrospective study transition indicated
            % by no probability and no regression data 
            if isempty(Transitions(t).Probability) && isempty(Transitions(t).RegParamVec),
                if ~FUTURE_IMPLEMENTATION,
                    error('The code below represents future implementation intended to support retrospective studies and has not been tested. To use this code, the variable FUTURE_IMPLEMENTATION should be set to true before running this code')
                end
                % Just a series of consistency checks
                if (Transitions(t).To == 0) || ~isempty(Transitions(t).RegCoefVec)|| ~isempty(Transitions(t).CovarianceMatrix),
                    error('Invalid definitions of transition %i',t);
                end
                    
                % Extract the index in the P matrix 
                j = ModelMatIndMap(Transitions(t).To);
                if Transitions(t).From == 0,
                    % If the from transition is NULL, the entire column 
                    % should be reset. Note that the diagonal should be
                    % modified. This will be performed later
                    P(:,j)=sym(0);
                else % this means a single member in the matrix reset
                    i = ModelMatIndMap(Transitions(t).From);
                    P(i,j)=sym(0);
                end
            end
        end  % Transitions loop
        
        
        % Update the diagonal of the matrix to reflect the fact that this
        % matrix represents a Markov model. This time it is necessary.
        for i=1:n,
            P(i,i)=sym(1)-sum(P(i,setdiff(1:n,i)));
        end
 
        %%%%%% End of Taking care of retrospective studies
 
        
        %%%%% Calculate timed probability matrices %%%%%
        
        % Init the first year matrices
        Pt={};
        Yt={};
        Pt{1} = P;
        Yt{1} = K*P*K1;
        
        % Loop through the rest of the study years and calculate matrices
        for i=2:Studies(StudyID).StudyLength,
            Pt{i}=Pt{i-1}*P;
            Yt{i}=K*Pt{i}*K1;
        end  % Time loop 
        
        %%%%% End of calculate timed probability matrices %%%%%
        
        %%% Find all relevant model transitions for this study %%%
        
        % Reset the transition names list
        RelevantTransNames = [];
        
        % Traverse study transitions
        for t = StudyTransInd,
            % Exclude transitions to and from Null states
            if Transitions(t).From ~= 0 && Transitions(t).To ~= 0,
                % Locate i,j indices from study transitions
                i=StudyMatIndMap(Transitions(t).From);
                j=StudyMatIndMap(Transitions(t).To);
                
                % Extract the i,j element for the time of the matrix
                MatElement = Yt{Studies(StudyID).StudyLength}(i,j);
                
                % Parse the expression to get the list of transition names
                [ SubProbExpr , NameList] = ParseMarkovTerm( MatElement, []);
                
                % Collect the transition names
                RelevantTransNames = union (RelevantTransNames,NameList);
                
            end
        end  % Transitions loop
 
        % Convert the collected names list to indices
        
        RelevantTransInd=zeros(1,length(RelevantTransNames));
        % For each model transition name
        for i = 1:length(RelevantTransNames),
            ReadTrans = sscanf(RelevantTransNames{i},'T%i_%i_%i');
            % The transition should be associated with the study while the
            % states should be the same as in the model
            StudyCheck = [Transitions.StudyID] == ReadTrans(1);
            FromCheck = [Transitions.From] == ReadTrans(2);
            ToCheck = [Transitions.To] == ReadTrans(3);
                
            % The study index should satisfy all three index values
            RelevantTransInd(i) = find(StudyCheck & FromCheck & ToCheck);
            
            if (RelevantTransInd(i)==0),
                error('Assertion error - The system extracted a non existing transition from the matrix - "%s", this is an internal system error', RelevantTransNames{i});
            end
        end % model transition name loop
        
        %%% End of Find all relevant model transitions for this study %%%
 
        %%%% Calculate the unified Study-Model table dimension set %%%%
 
        % Reset the study dimension list
        Ds={};
        StudyTableRanges={};
        
        % Determine if a Study is a Regression or a Categorical Study.       
        IsRegressionStudy = false;
        IsCategoricalStudy = false;
        % Traverse all Transitions in the study 
        for t = StudyTransInd,
              IsRegressionStudy = IsRegressionStudy || ~isempty(Transitions(t).RegCoefVec);
              IsCategoricalStudy = IsCategoricalStudy || ~isempty(Transitions(t).Probability);
        end % of Study Transitions Loop
 
        % Perform a Consistency Check to see if all study transitions are 
        % of the same kind.
        if ~xor (IsRegressionStudy , IsCategoricalStudy),
            error ('Study %s is defined ambigously by its transitions. A study should be either a categorical study or a regression study',Studies(StudyID).Name);
        end
 
        
        % Traverse all study transitions
        FirstTransition = true;
        for t = StudyTransInd,
            if IsCategoricalStudy,
                % If this is not an actual table - meaning it is a constant
                if ~isstruct(Transitions(t).Probability),
                    % Calculate the dimensions
                    Data = Transitions(t).Probability;
                    Dim={};
                    Dim{1} = 'Dummy';
                    Sizes = 1;
                    Ranges{1}=[-inf inf];
                else
                    % Calculate the dimensions
                    [Data, Dim, Sizes, Ranges] = TableInfo( Transitions(t).Probability );
                end
                % Exclude time
                [IsTimeDim TimeIndex] = ismember('Time',Dim);
                Dim = Dim(setdiff(1:length(Dim),TimeIndex));
                if isempty(Dim),
                    % In the case only a time dimension exists in the table
                    % then it was removed, then create the dummy dimension
                    Dim{1} = 'Dummy';
                    Sizes = 1;
                    Ranges{1}=[-inf inf];
                else
                    % If other dimensions existed
                    Sizes = Sizes(setdiff(1:length(Dim),TimeIndex));
                    Ranges = Ranges(setdiff(1:length(Dim),TimeIndex));
                end
                % Check if Ds is already defined
                if FirstTransition,
                    % If not defined, copy it from the set of dimensions
                    Ds=Dim;
                    StudyTableRanges=Ranges;
                    FirstTransition = false;
                else
                    % This means comparison of dimensions is needed.
                    % traverse all dimensions and see if these are defined 
                    % For each dimension: 
                    for DimIndex = 1:length(Dim),
                        % Extract dimension data
                        d=Dim(DimIndex);
                        Range = Ranges{DimIndex};
                        % check if already in list
                        [AlreadyInList, ArrayIndex] = ismember(d,Ds);
                        if AlreadyInList,
                            % Verify that ranges are compatible
                            if xor( any(isnan(StudyTableRanges{ArrayIndex})) , any(isnan(Range)))
                                error ('Conflict between continuous ranges and indexed ranges while trying to verify table ranges from transition %i in dimension %s', t, char(Dim));
                            end
                            StartInd = isnan(StudyTableRanges{ArrayIndex}(1))+1;
                            if any(StudyTableRanges{ArrayIndex}(StartInd:end) ~= Range(StartInd:end)),
                                error ('Conflict between ranges in the same study while trying to verify table ranges from transition %i in dimension %s', t, char(Dim));
                            end
                        else
                            error('Different dimension sets for transitions in the same study were detected in transition %i in dimension %s this is not supported by the system', t, char(Dim))
                        end
                    end
                    
                    Ds = Dim;
                    StudyTableRanges = Ranges;
                    FirstTransition = false;
                end
            else
                % This happens for a regression study
                Ds{1}='Dummy';
                StudyTableRanges{1}=[-inf inf];
            end            
        end % of study transition loop
         
        % Initialize the model dimension list
        D=Ds;
        UnifiedRanges=StudyTableRanges;
 
        % Traverse all the relevant model transitions 
        for t = RelevantTransInd,
            % Check if the transition contains categorical data
            if isstruct(Transitions(t).Probability),
                % Extract table info
                [ Data, Dim, Sizes, Ranges] = TableInfo( Transitions(t).Probability );
                
                % For each dimension 
                for DimIndex = 1:length(Dim),
                    % Exclude time
                    if ~strcmp(Dim(DimIndex),'Time'),
 
                        % Extract dimension data
                        d=Dim(DimIndex);
                        Range = Ranges{DimIndex};
 
                                                   
                        % Check if the dimension is in the list
                        [AlreadyInList, ArrayIndex] = ismember(d,D);
                        if AlreadyInList,
                            % Verify that ranges are compatible
                            if xor( any(isnan(UnifiedRanges{ArrayIndex})) , any(isnan(Range)))
                                error ('Conflict between continuous ranges and indexed ranges while trying to merge table range from transition %i in dimension %s', t, Dim);
                            end
                            % Unify range bounds (splitters). note that 
                            % union returns a sorted list. Since NaN should
                            % remain at the start and it is not sorted 
                            % differently
                            StartInd = isnan(UnifiedRanges{ArrayIndex}(1)) +1;
                            NewRange=union(UnifiedRanges{ArrayIndex}(StartInd:end),Range(StartInd:end));
                            UnifiedRanges{ArrayIndex}(StartInd:StartInd+length(NewRange)-1) = NewRange;
                        else
                            % Add to the list of Unified dimensions
                            ArrayIndex = length(D)+1;
                            D(ArrayIndex)=d;
                            UnifiedRanges{ArrayIndex} = Range;
                        end
                    else
                        error ('The time dimension has been detected in transition %i',t);
                    end
                end % Dimensions loop
            end
        end % study and relevant model transitions loop
        
        
        % Collect the Study table size
        StudyTableSize=[];
        for i=1:length(Ds),
            StudyTableSize(i)=length(StudyTableRanges{i})-1;
        end

        
        % Create an IndexLookup table to help traversing the table later on
        EmptyStudyTable = sym(MyInitMat(StudyTableSize));
        StudyTableInd = MyInitMat(StudyTableSize);
        StudyTableInd(:) = 1:prod(StudyTableSize);
        

        % If the dimension list is empty
        if isempty(D),
            D{1}='Dummy';
            UnifiedRanges{1}=[-inf inf];
        end
 
 
        % Collect the unified dimension table size
        UnifiedTableSize=[];
        for i=1:length(D),
            UnifiedTableSize(i)=length(UnifiedRanges{i})-1;
        end
        
        % Create an IndexLookup table to help traversing the table later on
        EmptyUnifiedTable = sym(MyInitMat(UnifiedTableSize));
        UnifiedTableInd = MyInitMat(UnifiedTableSize);
        UnifiedTableInd(:) = 1:prod(UnifiedTableSize);
 
        
        %%%% End Calculate the unified Study-Model table dimension set %%%%
        
 
        %%%%%%%%%%%%%%%%%%% Model Preprocessing %%%%%%%%%%%%%%%%%%%
 
        % Traverse all model transitions associated with the study
        for t = RelevantTransInd,
            % Create a new multidimensional table. The data in this table 
            % will be associated with the unified dimension list
            Transitions(t).Data=EmptyUnifiedTable;
 
            % If the table holds an expression. Note that this code will
            % take care of the case of a single coefficient as well
            if ~isstruct(Transitions(t).Probability),
                % Parse the expression string to find the names of the 
                % parameters and of the coefficients
                [ SubProbExpr , NameList] = ParseMarkovTerm( char(Transitions(t).Probability), []);
                % Extract parameters by removing coefficients
                UsedParameters = setdiff(NameList,Coefficients);
                
                % Traverse all table cells
                for c=UnifiedTableInd(:)',
                    % Extract the expression that may later be changed
                    ProbExpr = Transitions(t).Probability;
                    % For each parameter previously extracted
                    for ParamID = 1:length(UsedParameters),
                        Param = UsedParameters(ParamID);
                        % Check if in the constant list
                        [IsInConstantList IndexInConstantList] = ismember (Param,Constants);
                        [IsInUnifiedDimList IndexInUnifiedDimList] = ismember (Param,D);
 
                        % A consistency check
                        if IsInConstantList && IsInUnifiedDimList,
                           error('The parameter %s is used as a constant and as a dimension in a table', Param); 
                        elseif IsInConstantList,
                            % In case of a constant, use its value
                            ParamValue = ConstantValues(IndexInConstantList);
                        elseif IsCategoricalStudy,
                            % This code is still under development and
                            % intended to mean substitute covariates
                            % according to population information.
                            if ~FUTURE_IMPLEMENTATION,
                                error('The code below represents future implementation intended to support mean substitution of covariates and the concept is still under development and not been successfully tested. To use this code, the variable FUTURE_IMPLEMENTATION should be set to true before running this code')
                            end
                            % In case of a categorical study where the
                            % parameter is in the unified list or defined 
                            % by population, mean substitution will happen.
                            if IsInUnifiedDimList,
                                % If it is in the unified list:
                                % Extract cell index for each dimension
                                CellIndex = MyInd2Sub(UnifiedTableSize,c);
                                % Collect the cell range for each
                                % dimension/parameter and build a vector
                                CellRange = [UnifiedRanges{IndexInUnifiedDimList}(CellIndex(IndexInUnifiedDimList)), UnifiedRanges{IndexInUnifiedDimList}(CellIndex(IndexInUnifiedDimList)+1)];
                            else % it is defined by population
                                % extract the range from the population set
                                CellRange = [-inf inf];
                            end
                            % Using the previously defined range, calculate
                            % the mean value for this parameter by
                            % accessing the correct population
                            ParamValue = MeanCalcPop(Param,CellRange,Pop(Studies(StudyID).PopID));
                        else
                            % In case of a regression study, where the 
                            % parameter is not in the constant list, no
                            % substitution will occur
                            continue                            
                        end
                        % Substitute the parameter with the computed value
                        ProbExpr = subs(ProbExpr, Param, sym(ParamValue));
                    end % Parameter loop
                    % Place the expression in the table
                    Transitions(t).Data(c)=ProbExpr;
                end % Traversing cells in the table
            else % This means that the table holds categorical data
                
                %%%% Collapse Dimensions to create Source Table %%%%
                % Retrieve Probability table data
                [ProbData, ProbDim, ProbSizes, ProbRanges] = TableInfo( Transitions(t).Probability );
                
                %just copy the Probability table and its characteristics
                SourceTab=sym(Transitions(t).Probability.Data);
                SrcData = ProbData;
                SrcDim = ProbDim;
                SrcSizes = ProbSizes;
                SrcRanges = ProbRanges;
                % Init the index to 0 in the size of the sizes vector
                SrcIndex = ProbSizes*0;                   
 
                %%%%%%%% Expand to Fit Unified Dimensions %%%%%%%%%
 
                % Check source Table data
                if all(size(SourceTab)==1),
                    % If the Source table is truncated to 1 cell, than
                    % expanding the table, means just copying this cell to 
                    % the new data table
                    Transitions(t).Data(:) = SourceTab;
                else
                    % Otherwise, cell by cell complicated copying is 
                    % required during expansion    
                
                    % Traverse destination table cells
                    for c=UnifiedTableInd(:)',
                        % Extract cell index as a subscript vector
                        DestIndex = MyInd2Sub(UnifiedTableSize,c);
 
                        % For each dimension in D / DestIndex 
                        for i=1:length(D),
                            % Collect the dimension cell range high limit
                            CellRangeHigh = UnifiedRanges{i}(DestIndex(i)+1);
 
                            % Try to locate this dimension in the source 
                            % table
                            Dim=D(i);
                            [IsInSrcDimList IndexInSrcDimList] = ismember(Dim,SrcDim);
 
                            % If the dimension is in the source table
                            if IsInSrcDimList,
 
                                SrcRange = SrcRanges{IndexInSrcDimList};
                                % if indexed array, replace NaN with -inf
                                if isnan(SrcRange(1)),
                                    SrcRange(1)=-inf;
                                end
 
                                RangeIndexInSrcDimList = find ( CellRangeHigh > SrcRange(1:end-1) & CellRangeHigh <= SrcRange(2:end));
 
                                % This is an assertion check
                                if isempty (RangeIndexInSrcDimList),
                                    error('No source dimension is found for dimension %s while pre-processing transition %s in the model', char(Dim) ,char(RelevantTransNames{find(t==RelevantTransInd)}))
                                end
                                SrcIndex(IndexInSrcDimList) = RangeIndexInSrcDimList;
                            end
                        end % i in 1:length(D)
                        % Access the source table using the source index 
                        % vector Note that SrcIndex is a subscript index 
                        % vector that is converted to a single index into 
                        % a flat vector
                        FlatIndex = MySub2Ind( SrcSizes, SrcIndex );
                        SrcCellValue = sym(SourceTab(FlatIndex));
                        % now copy the correct cell from the source to the
                        % destination cell
                        Transitions(t).Data(c) = SrcCellValue;                
                    end % Traversing cells in the table
 
                end % If about SourceTab size
                
                %%%%%%%% End of Expand to Fit Unified Dimensions %%%%%%%%%
 
            end % If categorical or expression data
        end % Relevant model transitions for the study loop
 
        %%%%%%%%%%%%%%% End of Model Preprocessing %%%%%%%%%%%%%%%
 
        %%%%%%%%%%%%%%%% Study Preprocessing %%%%%%%%%%%%%%%%
        
        % Traverse all Transitions in the study 
        for t = StudyTransInd,
            % If transition is not from Null and the information in the 
            % study is not regression information 
 
            if Transitions(t).From ~=0 && isempty(Transitions(t).RegCoefVec),
 
                if ~isstruct(Transitions(t).Probability),
                                      
                    % If the data is a constant then create a new source
                    % table containing the probability data
                    SourceTab = struct('Data',[Transitions(t).Probability],'Dim',[Dummy,-inf,inf]);
                    
                    if ~all(size(Transitions(t).Probability) == [1 1]),
                        error('Non table information is not a constant in transition #%i',t );
                    end
 
                    if length(Ds)~=1 || ~ ismember('Dummy',Ds),
                        error ('Assertion error - Ds dimensions are not Dummy where the table is a constant')
                    end
                    
                else
                    % If the data is a table then just copy the table to
                    % the Source Variable
                    SourceTab = Transitions(t).Probability;
                end
 
                % Create the mapping index table. It will be the same size
                % as defined by Ds and each cell will start with an empty
                % matrix that will later be populated with cell index 
                % values from D corresponding to the cell in Ds
                Transitions(t).MappingIndexTable = cell(size(EmptyStudyTable));
 
                % Create the mapping cell prevalence table that provides
                % the prevalence of the model cell corresponding to the
                % study cell
                Transitions(t).MappingCellPrevalenceTable = EmptyUnifiedTable;
 
                
                % Retrieve source table data
                [SrcData, SrcDim, SrcSizes, SrcRanges] = TableInfo( SourceTab );
 
                % Find the 'Time' index in the source table
                [IsTimeInSrcTable IndexTimeInSrcTable] = ismember ('Time',SrcDim);
 
                % init the new table
                DsWithTime = Ds;
                StudyTableWithTimeSize = StudyTableSize;
                StudyWithTimeRanges = StudyTableRanges;
                
                if ~IsTimeInSrcTable 
                    % If time is not the transition just display an
                    % error that that time is expected and not present
                    % ignore the cases where initial count or the
                    % probability is defined implicitly in the study
                    if Transitions(t).To ~=0 && isstruct(Transitions(t).Probability),
                        error('Categorical data, that signifies an observed count lacks time dimension in a study in transition %i',t)
                    else
                        TimeSize = 1;
                        TimeRange = [NaN Studies(StudyID).StudyLength];
                    end
                else
                    % If time is define, add time to the dimensions
 
                    % Before than check that it is logical to find time in
                    % this transition
                    if Transitions(t).To ==0,
                        error('The Time dimension is detected in a transition signifying an initial number in transition %i',t)
                    end
                    
                    % Extract the Time dimension Size
                    TimeSize = SrcSizes(IndexTimeInSrcTable);
                    TimeRangeCell = SrcRanges(IndexTimeInSrcTable);
                    TimeRange = Cell2mat(TimeRangeCell);
 
                    % A consistency check note that Time range is a cell with a
                    % vector inside it and therefore the double indexing
                    if ~isnan(TimeRange(1)),
                        error ('The time dimension is indexed and therefore should start with NaN in transition %i',t)
                    end
                    % Another consistency check to see times are chronological.
                    % Again, note the unusual indexing to support the vector in
                    % the cell structure
                    if ~all(sort(TimeRange(2:end))==TimeRange(2:end)),
                        error ('The time indices should be sorted in ascending order to represent chronological order in transition %i',t)
                    end
 
                    % Another consistency check to see if the study and the
                    % last time specified in the table are the same
                    if TimeRange(end) ~= Studies(StudyID).StudyLength,
                        error ('The time specified in the table for transition #%i should be the same as the one specified in the Study', t)
                    end
 
                    % Handle the case that only the Time dimension was
                    % defined by the user
                    if length(SrcDim) == 1,
                        SrcDim(1) = {'Dummy'};
                        SrcSizes(1)= 1; 
                        SrcRanges(1)= {[-inf inf]};
                        SrcDim(2) = {'Time'};
                        SrcSizes(2)= TimeSize; 
                        SrcRanges(2)= TimeRangeCell;
                        IndexTimeInSrcTable = 2;
                    end
                    
                    % Add the 'Time' dimension at the end of the study list
                    % Note the fact that this is the last dimension as this
                    % will simplify accessing it later on.
                    DsWithTime(end+1)={'Time'};
 
                    % Create the same to handle time, with the correct size
                    StudyTableWithTimeSize(length(DsWithTime)) = TimeSize;
                    % Handle the range list
                    StudyWithTimeRanges(length(DsWithTime)) = TimeRangeCell;
                end
 
 
                % Create an IndexLookup table to help traverse the table later
                EmptyStudyWithTimeTable = sym(MyInitMat(StudyTableWithTimeSize));
                StudyTableWithTimeInd = MyInitMat(StudyTableWithTimeSize);
                StudyTableWithTimeInd(:) = 1:prod(StudyTableWithTimeSize);
 
                % Create a new multidimensional table. The data in this 
                % table will be associated with the study dimension list
                Transitions(t).Data=EmptyStudyWithTimeTable;
                % Record table definition info as well
                Transitions(t).StudyDim=DsWithTime;
                Transitions(t).StudySize=StudyTableWithTimeSize;
                Transitions(t).StudyRange=StudyWithTimeRanges;
                Transitions(t).TimeSize=TimeSize;
                Transitions(t).TimeRange=TimeRange;
 
                % Create the dimension mapping index vector
                DimMapIndexVec=zeros(1,length(Ds));
                
                % Traverse all dimensions in the source table
                for i = 1:length(SrcDim),
                    % Check where this dimension is in Ds+Time and record
                    % the index in the dimension mapping index vector
                    Dim = SrcDim(i);
                    [IsInDsWithTimeList IndexInDsWithTimeList] = ismember(Dim,DsWithTime);
                    if ~IsInDsWithTimeList,
                        error ('Assertion error - All source dimensions must be found in Ds+Time');
                    end
                    DimMapIndexVec(IndexInDsWithTimeList)=i;
                end
 
                % Traverse destination table cells
                for c=StudyTableWithTimeInd(:)',
 
                    % Extract cell index as a subscript vector
                    DestIndex = MyInd2Sub(StudyTableWithTimeSize,c);
 
                    % Extract the Source cell index by reshuffling the 
                    % vector to be arranged according to the dimension 
                    % mapping vector 
                    SrcIndex = DestIndex(DimMapIndexVec);
                    
                    % Convert the subscript index back to regular form
                    SrcIndexFlat = MySub2Ind(SrcSizes, SrcIndex);
                    
                    % Copy the cells from source to destination. This will
                    % essentially will reorder dimensions so to be consistent
                    Transitions(t).Data(c) = SourceTab.Data(SrcIndexFlat);
                    
                end % of traversing StudyTable cells
                
 
                % Reset the ranges
                StudyCellRange = [];
                ModelCellRange = [];
 
                % Traverse destination table cells
                for c=UnifiedTableInd(:)',
                    % Extract cell index as a subscript vector
                    ModelIndex = MyInd2Sub(UnifiedTableSize,c);
                    
                    % Reset the Study index vector
                    StudyIndexSub=zeros(1,length(Ds));
                    
                    % For each dimention in D 
                    for i=1:length(D),
                        % Collect the dimension cell range limits
                        ModelCellRangeLow = UnifiedRanges{i}(ModelIndex(i));
                        ModelCellRangeHigh = UnifiedRanges{i}(ModelIndex(i)+1);
 
                        % Collect the ranges in a matrix
                        ModelCellRange(i,:) = [ModelCellRangeLow ModelCellRangeHigh];
 
                        % Try to locate this dimension in the source table
                        Dim=D(i);
                        [IsInStudyDimList IndexInStudyDimList] = ismember(Dim,Ds);
 
                        % If the dimension is in the source table
                        if IsInStudyDimList,
 
                            StudyRange = StudyTableRanges{IndexInStudyDimList};
                            % if indexed array, replace NaN with -inf
                            StudyRangeCheck=StudyRange;
                            if isnan(StudyRange(1)),
                                StudyRangeCheck(1)=-inf;
                            end
 
                            RangeIndexInStudyDimList = find ( ModelCellRangeHigh > StudyRangeCheck(1:end-1) & ModelCellRangeHigh <= StudyRangeCheck(2:end));
 
                            % This is an assertion check
                            if isempty (RangeIndexInStudyDimList),
                                error('No source Data is found for dimension %s for value %g while pre-processing transition %i',char(Dim),ModelCellRangeHigh ,t)
                            end
                            StudyIndexSub(IndexInStudyDimList)=RangeIndexInStudyDimList;
 
                            % Organize the ranges in a matrix corresponding
                            % to the dimensions of the study table
                            StudyCellRange(IndexInStudyDimList,:) = StudyRange(RangeIndexInStudyDimList: (RangeIndexInStudyDimList+1));
                        end
                    end % i in 1:length(D)
 
                    % Update the Mapping index vector with the new study
                    % index after this is converted to flat form. Note
                    % that the model flat index is added to the proper
                    % study cell vector
                    StudyFlatIndex = MySub2Ind(StudyTableSize, StudyIndexSub);
                    Transitions(t).MappingIndexTable{StudyFlatIndex} = [ Transitions(t).MappingIndexTable{StudyFlatIndex} c ];
                    
                    
                    if isstruct(Transitions(t).Probability),
                        % If the probability was given as a table, then
                        % the study cell prevalence is not 1 and needs 
                        % calculation 
 
                        % A consistency check 
                        if any((StudyCellRange(:,1)-StudyCellRange(:,2))==0) ,
                            error('Assertion error - the extracted study cell range is invalid, this implies a source dimension that is not represented in the model table');
                        end
                                        
                        % Access the source table using the flat study 
                        % index. Note that SrcIndex is a subscript index 
                        % vector that is converted to a single index into 
                        % a flat vector
                        StudyCellValue = Transitions(t).Data(StudyFlatIndex);
 
                        % Extract the prevalence value associated with the
                        % the study dimensions without the time dimension
                        StudyPrevalence = sym(PrevalenceCalcPop( Ds, StudyCellRange, Pop(Studies(StudyID).PopID)));
                    else
                        % Otherwise, if the transition probability was
                        % provided as a constant, then the study
                        % prevalence is 1 and there is no need to calculate
                        % the prevalence. Actually this way of coding
                        % avoids potential problems with the code
                        StudyPrevalence = 1;
                        % The data is in one cell in this case
                        StudyCellValue = sym(Transitions(t).Data(1));
                    end
 
                    % Extract the prevalence value associated with the
                    % destination table without the time dimension
                    ModelPrevalence = sym(PrevalenceCalcPop( D, ModelCellRange, Pop(Studies(StudyID).PopID)));
                    
                    % calculate the prevalence between the cells
                    Prevalence = ModelPrevalence/StudyPrevalence;
 
                    % Now store the correct prevalence for mapping the 
                    % model cell to the study cell. This is stored in the
                    % mapping cell prevalence table for this transition
                    Transitions(t).MappingCellPrevalenceTable(c) = Prevalence;
                    
                end % Traversing cells in the table
                   
            end  % if From~=0...
            
        end % of Study Transitions Loop
        
        %%%%%%%%%%%%%%%% End of Study Preprocessing %%%%%%%%%%%%%%%%
        
        if DebugVec(1), 
            % Just some debug print outs. 
            disp ' The Model Regular States Probability Matrix R is: '
            disp (R)
            disp ' The Model Event States Probability Matrix E is: '
            disp (E)
            disp ' The Modified Event States Probability Matrix A is: '
            disp (A)
            disp ' The Model Multiple Event States Probability Matrix A^NumOfEvents is: '
            disp (MultA)            
            disp ' The Model Probability Matrix P is: '
            disp (P)
            disp ' The K matrix is: '
            disp (K)
            disp ' The K1 matrix is: '
            disp (K1)
 
            for i=1:length(Yt)
                disp (sprintf(' The Model probability matrix Pt(1) = (P^Time) for Time=%i is:',i))
                disp (Pt{i})
                disp (sprintf(' The study probability matrix Yt(1) = K*(P^Time)*K1 for Time=%i is:',i))
                disp (Yt{i})
            end
 
            disp ' The list of model transition names that are relevant to the study transitions'
            disp(RelevantTransNames)
 
            disp 'The following lines will display unified model table dimensions and ranges'
            for i=1:length(D)
                disp (D{i})
                disp (UnifiedRanges{i})
            end
 
            disp 'The following lines will display study table dimensions and ranges'
            for i=1:length(Ds)
                disp (Ds{i})
                disp (StudyTableRanges{i})
            end
 
            for t = RelevantTransInd,
                % Some debug prints
                disp (sprintf('After preprocessing the Model transition with index %i and named as %s Now holds the following Preprocessed unified model table:',t, char(RelevantTransNames{find(t==RelevantTransInd)})))
                disp (Transitions(t).Data)
            end
 
            if IsRegressionStudy,
                disp (sprintf ('The Regression Study transition with input index %i holds the following data:',t))
                disp(Transitions(t).Probability)
            else
                % If this is a categorical study
                for t = StudyTransInd,
                    % Some more debug prints
                    disp (sprintf ('After preprocessing the Categorical Study transition with input index %i: Now holds the following Preprocessed table data:',t))
                    disp(Transitions(t).Data)
                    disp (sprintf ('The transition holds the following index cell mapping to the model with the following prevalence values:'))
                    for c = StudyTableInd,
                        FlatIndicesToPrint = [Transitions(t).MappingIndexTable{c}];
                        PrevToPrint = Transitions(t).MappingCellPrevalenceTable(FlatIndicesToPrint);
                        disp (sprintf('FlatCellIndices = %s, Prevalence values = %s',mat2str(FlatIndicesToPrint),mat2str(double(PrevToPrint))))
                    end    
                end
            end
            disp 'End of Study Pre-Processing'
            disp '--------------------------------------'
            disp 'Starting to calculate the Likelihood Term for the Study'
        end % Debug prints for DebugVec(1)
 
        
        %%%%%%%%%%%%%%%% Calculate Study Likelihood %%%%%%%%%%%%%%%%
        
        % Reset the Study Log likelihood term
        StudyLogLikelihood = sym(0);
        
        if IsRegressionStudy,
            %%%%%%% Build Log Likelihood for Regression Studies %%%%%%%
            % Traverse all Transitions in the study 
            for t = StudyTransInd,
                if DebugVec(2), 
                    % Print debug information for building Likelihood
                    disp (sprintf('Processing Transition #%i:',t))
                end
                % If the transition holds data - not retrospective
                if ~isempty(Transitions(t).Probability) || ~isempty(Transitions(t).RegParamVec),
                    % For ease, collect the useful data in these variables
                    CovMat = sym(Transitions(t).CovarianceMatrix);
                    ParamVec = sym(Transitions(t).RegParamVec);
                    CoefVec = sym(Transitions(t).RegCoefVec);
                    nn = length(ParamVec);
                    % Several consistency Checks
                    if length(CoefVec)~=nn,
                       error ('size of coefficient vector and parameter vector are not the same in transition %i',t);
                    end
                    if any(size(CovMat)~=nn),
                       error ('size of the covariance matrix is not the same as the coefficient and parameter vectors in transition %i',t);
                    end
                                        
                    % Construct System Reserved Temporary Coefficients
                    SysResRegCoefVec=sym(zeros(1,nn));
                    % Loop to construct all vector members
                    for i=1:nn,
                        % Build the coefficient name
                        SymbVarName = sprintf('SysCoefT%i_%i_%i',Transitions(t).From,Transitions(t).To,i);
                        % Make it symbolic via an evaluation command
                        CommandStr = sprintf('syms %s ;',SymbVarName);
                        eval(CommandStr);
                        % Add it to the coefficient vector via a command
                        CommandStr = sprintf('%s;',SymbVarName);
                        SysResRegCoefVec(i) = sym(eval(CommandStr));
                    end % i=1:nn loop
 
                    % construct the system temporary coefficients for the
                    % model table.
 
                    % Create an empty table for the model
                    TempModelCoefTable = EmptyUnifiedTable;
                    % Create the temporary individual count table
                    TempIndividualCountTable = EmptyUnifiedTable;
                    % Loop to construct Temporary study term accumulation 
                    % vector members that each is a table
                    TempStudyTermAccTable={};
                    for i=1:(2*nn+1),
                        TempStudyTermAccTable{i} = EmptyUnifiedTable;                    
                    end
 
                    SysCoefficientList = [];
                    % Now populate its cells with the coefficients
                    for c=UnifiedTableInd(:)',
                        % Build coefficient name
                        SymbVarName = sprintf('SysModelC%i',c);
                        % Make it symbolic via a command
                        CommandStr = sprintf('syms %s ;',SymbVarName);
                        eval(CommandStr);
                        % Set the Table cell to the coefficient
                        TempModelCoefTable(c) = eval(SymbVarName);
                        SysCoefficientList{end+1} = SymbVarName;
                    end % Traversing cells in the table                       
                    
                   
                    % Handle Population sets based on Distribution functions
                    if Studies(StudyID).PopID == 0,
                        PopulationDim = {'Dummy'};
                        PopulationData = zeros(1);                        
                    elseif isempty(Pop(Studies(StudyID).PopID).Data),
                        % Populations based on distribution functions are
                        % not supported
                        error('Population sets based on Distribution functions are not yet supported')                        
                    else
                        PopulationDim = {Pop(Studies(StudyID).PopID).Columns.Dim};
                        PopulationData = Pop(Studies(StudyID).PopID).Data;
                    end
 
                    % If this point was reached, this means the population 
                    % is defined by data
 
                    %%%%% Prepare the study correlation expression %%%%%
 
                    % Construct the study correlation term for 
                    % this cell. This is only the base construction that
                    % will later be changed by introducing population 
                    % individual values into the parameters
 
                    % Examine the regression coefficient type
                    switch Transitions(t).RegFunctionType
                        case 1, % Exponential form
                            Func = exp(-SysResRegCoefVec*ParamVec.');
                        case 2, % Linear form
                            Func = SysResRegCoefVec*ParamVec.';
                        case 3, % UKPDS form
                            if ParamVec(nn) ~= Studies(StudyID).StudyLength,
                                error('In transition %i, The last member of the parameter vector in UKPDS form should be equal to the study length',t);
                            end
                            Func = 1-exp(-(1-SysResRegCoefVec(nn)^ParamVec(nn))/(1-SysResRegCoefVec(nn))*prod(SysResRegCoefVec(1:(nn-1)).^ParamVec(1:(nn-1))));
                        case 11, % One Minus Exponential form
                            Func = 1 - exp(-SysResRegCoefVec*ParamVec.');
                        otherwise,
                        error ('Unrecognized Regression function type in transitions %i',t);
                    end
 
                    % This generally defines the study correlation term Ls 
                    Ls = sym(Func);
                    
                    % Copy Ls to StudyCorrelationNoConstant to allow 
                    % construction of the Term by parameter substitution of
                    % the constants.
                    StudyCorrelationNoConstant = Ls;
                    
                    % Prepare the Symbolic vector that later will be
                    % exchanged with replacement values
                    ParamReplacementVec = sym([]);
                    ParamReplacementIndex = [];
                    for i=1:nn,
                        Param = ParamVec(i);
 
                        % Check if this is an actual constant number
 
                        SymbolicVar = findsym(sym(Param));
 
                        if isempty(SymbolicVar),
                            % This means the data is probably a number
                            % and a run time error will be generated 
                            % and caught to give a meaningful user feedback
                            try
                                % Try symbolic to number conversion
                                ParamValue = double(Param);
                                % note that there is no need to perform
                                % substitution in this case as the 
                                % number will remain a number
                            catch
                                % Write a meaningful error
                                error ('The parameter %s defined in the regression parameter vector in transition %i is invalid',Param,t);
                            end
                        else % Meaning a symbolic variable exists
                            % Check that there is only one parameter in
                            % the parameter vector
                            if ~isempty(findstr(',',SymbolicVar)), 
                               error('The regression parameter vector member %i is composed of more than one parameter: "%s"',i, char(Param)); 
                            end
 
                            [IsInConstantList IndexInConstantList] = ismember (char(SymbolicVar),Constants);
                            [IsInPopulationSet DimIndexInPopulationSet] = ismember (char(SymbolicVar),PopulationDim);
                            [IsInCoefficients IndexInCoefficients] = ismember (char(SymbolicVar),Coefficients);
 
                            % A consistency check
                            if IsInConstantList && IsInPopulationSet,
                                error('The parameter %s is defined both as a constant and part of a population set', char(SymbolicVar)); 
                            elseif IsInConstantList,
                                % Convert the constant to the number
                                ParamValue = ConstantValues(IndexInConstantList);
                                % immediately replace this value in the study
                                % correlation expression
                                StudyCorrelationNoConstant = subs(StudyCorrelationNoConstant, SymbolicVar, sym(ParamValue));
                            elseif IsInPopulationSet, 
                                % In case of data in the population set,
                                % then record the symbolic parameter and
                                % the index in the replacement vectors
                                ParamReplacementVec(end+1) = sym(SymbolicVar);
                                ParamReplacementIndex(end+1) = DimIndexInPopulationSet;
                            elseif ~IsInCoefficients,
                                error('Non coefficient parameter detected')    
                            end
                        end
                    end % of loop i=1:nn
 
                    % Create the perturbation vectors matrix
                    % init the perturbation matrix
                    PerturbationVectorMat = [ zeros(1,nn) ; eye(nn)*h ; -eye(nn)*h] ;
                    
                    % Loop through possible perturbations of the
                    % system coefficients and create the perturbed
                    % vector
                    LsPerturbedVector = sym([]);
                    for PerturbationIndex = 1:(2*nn+1),
                        % Substitute coefficient values into the study
                        % correlation 
                        SubsLs = subs (StudyCorrelationNoConstant, SysResRegCoefVec, CoefVec+PerturbationVectorMat(PerturbationIndex,:));
                        LsPerturbedVector(PerturbationIndex) = SubsLs;
                    end % of PerturbationIndex loop
     
                    %%%% End of prepare the study correlation expression %%%%

                    %%%% Modify correlation term using population %%%%
                    
                    % Iterate through individuals in the data set
                    for IndividualIterator = 1:size(PopulationData,1),
                        IndividualRecord = PopulationData(IndividualIterator,:);
                        ContinueToNextIndividual = false;
 
 
                        %%%%%%%% Construct Lm %%%%%%%%
                        
                        % Construct Lm the Model correlation Term
 
                        % If the dimension is not a 'dummy' then 
                        % the dimension need substitution from the
                        % population set.
                        % First locate the cell index of the data in the 
                        % unified table, for this traverse the table
                        % dimensions
                        CellIndex = MyInitMat(length(D))';
 
                        for IndexToIndex = 1:length(D),
                            CurrentDim = D(IndexToIndex);
                            if strcmp('Dummy',CurrentDim),
                                % For the dummy dimension the cell
                                % index is 1
                                CellIndex(IndexToIndex) = 1;
                            else
                                % For any other dimension, the index
                                % actually requires calculation by
                                % accessing the population data
                                [IsInPopulationSet DimIndexInPopulationSet] = ismember (CurrentDim,PopulationDim);
                                if ~IsInPopulationSet, 
                                    error ('Dimension %s is not in population set %i defined in study ID %i', char(CurrentDim), Studies(StudyID).PopID, StudyID)
                                end
 
                                RangeVector = UnifiedRanges{IndexToIndex};
                                CellVal = IndividualRecord(DimIndexInPopulationSet);
                                if isnan (CellVal),
                                    ContinueToNextIndividual = true;
                                    break;
                                end;
                                if isnan(RangeVector(1)),
                                    % This means discrete range, avoid NaN
                                    CurrentDimIndex = find(RangeVector(2:end) == CellVal);
                                else
                                    % This means a continuous range
                                    CurrentDimIndex = find(RangeVector(1:end-1) < CellVal & CellVal <= RangeVector(2:end));
                                end    
                                % Check if out of bounds
                                if CurrentDimIndex == 0 || CurrentDimIndex > UnifiedTableSize(IndexToIndex)
                                    error ('Out of bounds value %g, found in the population set for dimension %s in population set %i in study ID %i',CellVal, CurrentDim, Studies(StudyID).PopID, StudyID)
                                end
 
                                CellIndex(IndexToIndex) = CurrentDimIndex;
                            end
                        end  % For IndexToIndex
                        % Check if to skip to the next individual
                        if ContinueToNextIndividual,
                            continue;
                        end
                        % correct the index size to be 2 in case of 1
                        % dimension to allow proper indexing
                        if length(CellIndex) == 1,
                           CellIndex = [CellIndex 1]; 
                        end
                        
                        % Knowing the cell index, use it to calculate the
                        % transition probability expression
                        FlatIndex = MySub2Ind( UnifiedTableSize , CellIndex );
 
                        % Set the model correlation term Lm
                        Lm = sym(TempModelCoefTable(FlatIndex));
                        
                        %%%%%%%% End Construct Lm %%%%%%%%%%%%%%
                   
                        %%%%% Substitute individual data in Ls %%%%
 
                        % Decide if any substitution is required
                        if length(ParamReplacementIndex) > 0
                            % Create the value substitution vector by using the
                            % previously calculated index vector of the
                            % parameteres
                            ValueReplacementVec = IndividualRecord(ParamReplacementIndex);
 
                            % If any of this parameteres are missing, i.e. NaN
                            % is detects in the vector, then skip to the next
                            % individual
                            if any(isnan(ParamValue)),
                                continue;
                            end
                            
                            % Perform the replacement to construct Ls
                            SubsLsPerturbedVector = subs (LsPerturbedVector, ParamReplacementVec, ValueReplacementVec);
                        else
                            % No susbstitution from population is needed
                            SubsLsPerturbedVector = LsPerturbedVector;
                        end
                        
                        % Convert Ls to a number with double precision to
                        % avoid complicated expressions. This should no
                        % affect overall precision as population data and
                        % regression coefficients, ussually are supplied to
                        % machine precision.
                        SubsLsPerturbedVectorNumeric = double(SubsLsPerturbedVector);
                        
                        % Loop through possible perturbations of the
                        % system coefficients 
                        for PerturbationIndex = 1:(2*nn+1),
                            % Accumulate the study term in the correct
                            % table cell according to the perturbation and
                            % the flat Index previously calcualted
                            TempStudyTermAccTable{PerturbationIndex}(FlatIndex) = TempStudyTermAccTable{PerturbationIndex}(FlatIndex) + SubsLsPerturbedVectorNumeric(PerturbationIndex);
                        end % of PerturbationIndex loop
                        
                        %%%% End of Substitute individual data in Ls %%%%%
                        
                        % Increase the individual count
                        TempIndividualCountTable(FlatIndex) = TempIndividualCountTable(FlatIndex) +1 ;
                        
                    end % of for IndividualIterator

                    %%%% End of Modify correlation term using population %%%%
                    
                    %%%% Correlation equation and Jacobian solution %%%%
                    
                    % Define the size of the vector as the table size 
                    mm = length(UnifiedTableInd(:));
                    
                    ExpectedSol = zeros(1,mm);
                    ForwardPerturbationSol = zeros(mm,nn);
                    BackwardPerturbationSol = zeros(mm,nn);
                    
                    % Define the output vector
                    CoefficientsUsed = sym([]);
                    % Traverse all table cells
                    for c=UnifiedTableInd(:)',
                        % Rewrite the table as a vector
                        CoefficientsUsed(c) = TempModelCoefTable (c);
                        % Check that there was enough information in the
                        % population to calculate this cell
                        if TempIndividualCountTable(c)==0,
                            error ('The population does not provide sufficient information to construct the study term for cell #%i' , c);
                        end
 
                        
                        % Loop through all perturbations
                        for PerturbationIndex = 1:(2*nn+1),
                            CorrelationSolution = TempStudyTermAccTable{PerturbationIndex}(c) / TempIndividualCountTable(c);
 
                            if PerturbationIndex==1
                                ExpectedSol(c) = CorrelationSolution;
                            elseif PerturbationIndex <= nn+1,
                                ForwardPerturbationSol(c,PerturbationIndex-1) = CorrelationSolution;
                            else
                                BackwardPerturbationSol(c,PerturbationIndex-nn-1) = CorrelationSolution;
                            end
                        end % of PerturbationIndex loop
                    end % Traversing cells in the table
 
                    Jacobian = (ForwardPerturbationSol - BackwardPerturbationSol) / (2*h);
                    
                    % For convenience calculate the difference between the 
                    % coefficients and the expected values
                    CoefDiff = CoefficientsUsed - ExpectedSol;
                    
                    % Calculate the new covariance matrix by using the
                    % calculated Jabocian and the old covariance matrix
                    NewCovMat = Jacobian * CovMat * Jacobian';
                    
                    %%%% End of Correlation equation and Jacobian solution %%%%
                   
                    %%%% Construct Lt - translated likelihood term %%%%
 
                    % Build the un-substituted translated regression 
                    % likelihood term
                    Lu = (-sym(mm)/sym(2))*sym(log(2*pi)) - 0.5*log(det(NewCovMat)) -sym(0.5)*CoefDiff*inv(NewCovMat)*(CoefDiff).';
 
                    % Extract the probability expression from the Y
                    % matrix for the end time of the study
                    TimeEndProbMat = Yt{Studies(StudyID).StudyLength};
                    % calculate indices to access matrix
                    MatFromInd = StudyMatIndMap(Transitions(t).From);
                    MatToInd = StudyMatIndMap(Transitions(t).To);
                    % extract the expression from the matrix
                    ProbExprIn = TimeEndProbMat(MatFromInd, MatToInd);
                    
                    % Substitute its coefficients by the actual model
                    % coefficients. Start by traversing the used coefficient 
                    % vector
                    CoefVectorToSubstitute = sym(zeros(1,mm));
                    for c=UnifiedTableInd(:)',
                        CoefToProcess = CoefficientsUsed(c);
                        CoefString = char(CoefToProcess);
                        % Make an assertion check
                        if ~ all( CoefString(1:9) == 'SysModelC' )
                            error ( 'Unrecognized coefficient %s encountered while writing the likelihood expression', CoefString)
                        end
                        CoefIndex = str2double(CoefString(10:end));
 
                        % Extract cell index as a subscript vector
                        CellIndex = MyInd2Sub(UnifiedTableSize,CoefIndex);
                        
                        % Parse the expression and perform substitutions
                        % of the appropriate transitions
                        SubProbExpr = ParseMarkovTerm(ProbExprIn, CellIndex, Transitions);
                        CoefVectorToSubstitute(c) = sym(SubProbExpr);
                    end % Traversing cells in the table
                    
                    % Create the substituted translated likelihood
                    % expression by substituting the temporary coefficients
                    % with the coefficient vector to substitute
                    
                    Lt = subs (Lu,CoefficientsUsed,CoefVectorToSubstitute);
                    
                    %%%% End of Construct Lt - translated likelihood term %%%%
 
                    
                    % Update the General Log Likelihood expression, by 
                    % adding to it the log likelihood Term for this 
                    % transition
                    StudyLogLikelihood = StudyLogLikelihood + Lt;
                    
                    if DebugVec(2), 
                        % Print debug information for building Likelihood
                        disp 'Constructing Regression Study Log Likelihood:'
                        disp 'The new expected value vector:'
                        disp (ExpectedSol)
                        disp 'The Jacobian for these values with respect to original coefficients is:'
                        disp (Jacobian);
                        disp 'The New Covariance Matrix:'
                        disp (NewCovMat)
                        disp 'The Un-substituted translated log likelihood term (Lu):'
                        disp (Lu)
                        disp 'Substituted translated log likelihood term (Lt):'
                        disp (Lt)
                    end                    
                end % If data is meaningful
            end % of Study Transitions Loop
            
            %%%%%%% End of Build Log Likelihood for Regression Studies %%%%%%%
        elseif IsCategoricalStudy,
            %%%%%%% Build Log Likelihood for Categorical Studies %%%%%%%%%%
          
            if DebugVec(2), 
                % Print debug information for building Likelihood
                disp 'Constructing Categorical Study Log Likelihood:'
            end
            % Traverse all Transitions in the study 
            for t = StudyTransInd,
                if DebugVec(2), 
                    % Print debug information for building Likelihood
                    disp (sprintf('Processing Transition #%i:',t))
                end
                % If the transition holds data - not retrospective and it
                % does not represent starting sample size
                if ~isempty(Transitions(t).Probability) && (Transitions(t).To~=0),
 
                    % Init the transition Log Likelihood term
                    TransLogLikelihood = sym(0);
 
                    if isstruct(Transitions(t).Probability),
                        % This means that a table is defined rather than a
                        % constant probability.
                        % Therefore, time should be extracted from the
                        % table
                        
                        % Retrieve the Study table data 
                        StudyData = Transitions(t).Data ;
                        % Record table definition info as well
                        StudyDim = Transitions(t).StudyDim;
                        StudySizes = Transitions(t).StudySize;
                        StudyRanges = Transitions(t).StudyRange;
                    end
                    
                    % Extract the time range in any case
                    TimeSize = Transitions(t).TimeSize;
                    TimeRange = Transitions(t).TimeRange;
                    
                    % Traverse all study table cells without time
                    for c=StudyTableInd(:)',
 
                        % Set the previous probability expression and the
                        % previous observed data count to zero
                        OldModelCellExpressionAccumulator=sym(0);
                        OldObservedCount=sym(0);
 
                        % Extract cell index as a subscript vector
                        StudyCellIndexNoTime = MyInd2Sub(StudyTableSize,c);
 
                        % Init the Study table cell Log Likelihood term
                        StudyCellLogLikelihood = sym(0);
                        
                        % Loop through all outcome times
                        for TimeInd = 1:TimeSize
                            % Skip the NaN to extract the time
                            TheTime = TimeRange(TimeInd+1);
                            
                            % Extract the probability expression from the Y
                            % matrix for the appropriate year of the study
                            TimeProbMat = Yt{TheTime};
                            % calculate indices to access matrix
                            MatFromInd = StudyMatIndMap(Transitions(t).From);
                            MatToInd = StudyMatIndMap(Transitions(t).To);
                            % extract the expression from the matrix
                            ProbExprIn = TimeProbMat(MatFromInd, MatToInd);
                            
                            % Initialize the model cell expression
                            ModelCellExpressionAccumulator = sym(0);
                            
                            % Traverse all model cells associated with the
                            % study cell according to the mapping index
                            % table prepared during preprocessing
                            for ModelCellIndexFlat = Transitions(t).MappingIndexTable{c}
                                % Extract cell index as a subscript vector
                                ModelCellIndex = MyInd2Sub(UnifiedTableSize,ModelCellIndexFlat);
 
                                % If the index is of size 1, then add
                                % another 1 at the end of the vector to fit
                                % the matlab convention that sizes a
                                % constant or a vector as a matrix
                                if length(ModelCellIndex) == 1,
                                    ModelCellIndex = [ModelCellIndex 1];
                                end
                                
                                % Parse the expression and perform substitutions
                                % of the appropriate transitions
                                [SubProbExpr , NameList] = ParseMarkovTerm(ProbExprIn, ModelCellIndex, Transitions);
 
                                % Convert the string to a symbolic expression
                                SubProbExprSym = sym(SubProbExpr);
                                
                                % Extract the model cell prevalence 
                                Prevalence =  Transitions(t).MappingCellPrevalenceTable(ModelCellIndexFlat);
 
                                ModelCellExpressionAccumulator = ModelCellExpressionAccumulator + Prevalence*SubProbExprSym;
                            end
                            
 
                            % Extract the observed count for this time and
                            % cell from the study table. Note that the
                            % Table is constructed such that it has all the
                            % study dimensions are first and time is the
                            % last dimension. Therefore it is easy to
                            % update the subscript index to include the 
                            % time index.
                          
                            FlatStudyIndex = MySub2Ind( [StudyTableSize TimeSize] , [StudyCellIndexNoTime TimeInd] );
                            ObservedCount = sym(Transitions(t).Data(FlatStudyIndex));
 
                            % Update the Cell log likelihood Term 
                            
                            % Check if there is a need to process the cell
                            CountDiff = ObservedCount - OldObservedCount 
                            if CountDiff == 0,
								disp (sprintf('Warning: Observed Count difference is 0, skipping cell calculation for Study Cell flat index %i, Time = %i', c, TheTime))
								TimeStepLogLikelihood = sym(0)
                            else,
								TimeStepLogLikelihood = CountDiff * log(ModelCellExpressionAccumulator - OldModelCellExpressionAccumulator);
							end
                            StudyCellLogLikelihood = StudyCellLogLikelihood + TimeStepLogLikelihood;
                            
                            if DebugVec(2), 
                                % Print debug information for building Log Likelihood
                                disp (sprintf('Study Time Step - Log Likelihood term for Time = %i:',TheTime))
                                disp (TimeStepLogLikelihood)
                            end                    
 
                            
                            % Update the Previous Observed count and the
                            % previous probability expressions to the
                            % current values to serve in the next iteration
                            % or when the loop is exited
                            OldModelCellExpressionAccumulator = ModelCellExpressionAccumulator;
                            OldObservedCount = ObservedCount;                            
                        end %PointInTime Loop
                        
                        % Find the transition that corresponds to the
                        % sample size for this transition. This transition
                        % belongs to this study and start from the same
                        % state, however it ends in a NULL state. 
                        
                        StudyCheck = [Transitions.StudyID] == StudyID;
                        FromCheck = [Transitions.From] == Transitions(t).From;
                        ToCheck = [Transitions.To] == 0;
 
                        % The study index should satisfy all three indices
                        SampleTransIndex = find(StudyCheck & FromCheck & ToCheck);
 
                        % Data consistency checks
                        if length(SampleTransIndex)>1,
                            error('Input data has more than one transition defining sample size for study #%i - "%s" from state #%i to state NULL ',StudyID, Studies(StudyID).Name, Transitions(t).From);                    
                        end
 
                        if isempty(SampleTransIndex),
                            error('no sample size transition is defined for study #%i - "%s" from state #%i to state NULL ',StudyID, Studies(StudyID).Name, Transitions(t).From);                    
                        end
 
                        SampleSize = sym(Transitions(SampleTransIndex).Data(c));
                        
 
                        % Update the Cell likelihood
                        SampleCountDiff = (SampleSize - OldObservedCount)
                        % Check if there is a need to process the cell
                        if SampleCountDiff == 0,
							disp (sprintf('Warning: Observed Sample Count difference is 0, skipping sample calculation for Study Cell flat index %i', c))
							SampleSizeLogLikelihood = sym(0)
                        else,
							SampleSizeLogLikelihood = SampleCountDiff * log(sym(1) - OldModelCellExpressionAccumulator);
						end
                        
                        StudyCellLogLikelihood = StudyCellLogLikelihood+SampleSizeLogLikelihood;
                        TransLogLikelihood = TransLogLikelihood+StudyCellLogLikelihood;
                        FunctionListStudyCellLogLikelihood{end+1} = StudyCellLogLikelihood;
                        
                        if DebugVec(2), 
                            % Print debug information for building Likelihood
                            disp 'Sample Size Log Likelihood term:'
                            disp (SampleSizeLogLikelihood)
                            disp 'Unified Table Cell Log Likelihood term:'
                            disp (StudyCellLogLikelihood)
                        end                    
 
                    end % Traversing cells in the table
                    
                    % Update the overall Log Likelihood term
                    StudyLogLikelihood = StudyLogLikelihood+TransLogLikelihood;
                    FunctionListTransLogLikelihood{end+1} = StudyLogLikelihood;
                    
                    if DebugVec(2), 
                        % Print debug information for building Log Likelihood
                        disp (sprintf('Transition Log Likelihood Term for transition #%i:',t))
                        disp (TransLogLikelihood)
                    end
                    
                end % If Meaningful Transition
            end % of Study Transitions Loop
        
            %%%%%%% End of Build Log Likelihood for Categorical Studies %%%%%%%
        end
 
        % Update the General Log Likelihood
        L = L + StudyLogLikelihood;
        FunctionListStudyLogLikelihood{end+1} = StudyLogLikelihood;
 
 
        disp (sprintf('The Log Likelihood term for study %i - "%s" is therefore:', StudyID,Studies(StudyID).Name))
        disp (StudyLogLikelihood)
 
        disp 'End of study output'
        disp '--------------------------------------'
        
    end % Studies loop
    
    
end % of Sub-Process loop
 
disp '--------------------------------------'
disp '--------------------------------------'
disp 'The Log Likelihood expression for the Entire Project is:'
disp (L)
 
 
 
if DebugVec(2), 
    % Show the time required for calculation
    disp 'Calculation time was:'
    toc
end
disp '--------------------------------------'
disp '--------------------------------------'
 
 
 
% Perform a consistency check to see if all studies have been used
        
UnusedStudyIDs = find(StudiesUsed==0);
% Remove the model from the list of studies before performing the check
UnusedStudyIDs = setdiff(UnusedStudyIDs,ModelID);
if ~isempty(UnusedStudyIDs)
   error ('There were Studies that were not used - possibly due to incompatibility with the sub-process structure of the model. These studies ID numbers are: [%s]',num2str(UnusedStudyIDs));
end
