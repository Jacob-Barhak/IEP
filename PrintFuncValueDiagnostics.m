function [FuncOut] = PrintFuncValueDiagnostics (X, ExpressionOfX, PrintInputText, PrintOutputText, OnlyIfWarning)
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
%function [FuncOut] = PrintFuncValueDiagnostics (X, ExpressionOfX, PrintInputText, PrintOutputText, OnlyIfWarning)
%This function can be used as a wrapper function to another function that 
%calculates the expression ExpressionOfX with the input value X so that the
%input value X and the Output Value FuncOut are printed. 
%If PrintInputText and PrintOutputText are non empty strings that input and
%output are printed and shown to them.
%The function returns the value of the expression FuncOut
 
FuncOut = ExpressionOfX;
 
PrintStr = '';
 
if ~isempty(PrintInputText),
    PrintStr = strcat(PrintStr , PrintInputText, mat2str(X));
end
if ~isempty(PrintOutputText),
    PrintStr = strcat(PrintStr , PrintOutputText, mat2str(ExpressionOfX));
end
 
if OnlyIfWarning,
    if isempty(lastwarn()),
        return
    else
        % just clear the warning flag
        lastwarn('');        
    end
end
 
disp (PrintStr)
return

