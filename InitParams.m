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
% You should have received a copy of the GNU General Public License along
% with the Indirect Estimation Prototype.  If not, see
% <http://www.gnu.org/licenses/>.
%
% This script initializes the Params vector and creates associated symbolic
% variables that correspond to the parameters in Param. This file should be
% used after defining constants, covariates, and coefficients. It should be
% called before calling CalcLikelihood.m and before using the symbolic
% counterparts of any of the parameters. 
% This script also defines the symbolic variables 'Time' and 'Dummy'
%
% INPUT:
% Constants - All the Constant names as stings in a cell vector.
% Covariates - All the Covariate names as stings in a cell vector.
% Coefficients - All the coefficient names as stings in a cell vector.
%
% OUTPUT:
% Params - The union of Covariates Coefficients and Constants. 
% (Symbolic Variables) - A set of symbolic variables with the same names
%                        defined in param. In addition the Symbolic 'Dummy'
%                        and 'Time' are defined
% (Printout of Copyright and License) - This file also prints a the
%                                       Copyright and License notice.

% Print Copyright notice
disp(' The Indirect Estimation Prototype (Version 0.39)')
disp(' Copyright (C) 2007,2009 The Regents of the University of Michigan')
disp(' Initially developed by Deanna Isaman, Jacob Barhak')
disp(' This program comes with ABSOLUTELY NO WARRANTY; for details see the file')
disp(' Copying.txt. This is free software, and you are welcome to redistribute')
disp(' it under certain conditions; see the file Copying.txt for details. ')

% Define the system reserve symbolic parameters Time and Dummy
syms Time Dummy

Params = union (Covariates, Coefficients);
Params = union (Params, Constants);

% Create symbolic variables from the parameters
for i=1:length(Params)
   Str=cell2mat(Params(i));
   Expression =  sprintf ( 'syms %s', Str);
   eval (Expression);
end