Copyright (C) 2007,2009 The Regents of the University of Michigan
Initially developed by Deanna Isaman, Jacob Barhak

This file is part of the Indirect Estimation Prototype. The Indirect
Estimation Prototype is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

The Indirect Estimation Prototype is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ADDITIONAL CLARIFICATION

The Indirect Estimation Prototype is distributed in the hope that it will 
be useful, but "as is" and WITHOUT ANY WARRANTY of any kind, including any 
warranty that it will not infringe on any property rights of another party 
or the IMPLIED WARRANTIES OF MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. THE UNIVERSITY OF MICHIGAN assumes no responsibilities with  
respect to the use of the Indirect Estimation Prototype.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CONTENTS:
1. INTRODUCTION
2. INSTALLATION
3. HOW TO USE THE SOFTWARE
4. FILES
5. KNOWN BUGS AND FEATURES
6. VERSION HISTORY
7. CONTACT INFORMATION
8. ACKNOWLEDGEMENTS


1. INTRODUCTION:

This software is the Indirect Estimation Prototype.

This software enables estimating Markov model coefficients given a set of 
studies with sufficient data observations and associated population data.
The software implements a method referred to as the Lemonade Method since 
when study data give you lemons, make lemonade.

The system performs the following tasks: 
1) Pre-Processing the input data
2) Constructing the Likelihood expression symbolically
3) Optimizing the Likelihood to extract optimal coefficient values

Note that this Matlab prototype version is intended to be migrated to
Python environment. At the time of this release, the Python version of the
Indirect Estimation and Simulation Tool (IEST) is also released and it also
implements this estimation Prototype including GUI.

The Method behind this Prototype is described in part in the following 
Publications:
- Isaman, D.J.M., Herman, W.H., & Brown, M.B. A discrete-state and discrete-
  time model using indirect estimates. Statistics in Medicine 2006 25(march):
  1035-1049. DOI: 10.1002/sim.2241
- D.J.M. Isaman, J. Barhak, W. Ye: Indirect Estimation of a Discrete-State 
  Discrete-time model using Secondary Data Analysis of Regression Data. 
  Statistics in Medicine Volume 28, Number 16, Pages 2095 – 2115. 
  DOI:10.1002/sim.3599 

For more information about the software see the following web site:
http://www.med.umich.edu/mdrtc/cores/DiseaseModel



2. INSTALLATION:

The software was tested on Matlab Version 7.8.0.347 (R2009a) 64-bit (glnxa64) 
on a Linux machine. Previous versions were designed and tested on Matlab 
Version 7.2.0.232 (R2006a) on a PC computer with Windows XP Professional SP3. 
The software and uses the Symbolic toolbox and the statistics toolbox within 
Matlab. 

To install the software copy all its files to the same directory. 



3. HOW TO USE THE SOFTWARE:

To verify that this software works properly, run the script TestScript.m
from within the Matlab Environment. The results generated in 
TestingOutput.txt should be similar to the expected outcomes in
ExamplesForIndirectEstimation.pdf . 

To create an example of your own:
1. Create an input script similar to ProjectDataExample*.m
2. Run the script in the same directory as all the other files.



4. FILES:

The software contains the following files:

~~~~~~~~~~~~~~ Major Scripts ~~~~~~~~~~~~~~
The following files contain the most major scripts to be run by the user:

CalcLikelihood.m : A main script that generates the Likelihood expression
OptimizeLikelihood.m : The main script to optimizing likelihood expression
InitParams.m : A script that initializes the parameters
LikelihoodGraphs.m : A script to visually represent the results

~~~~~~~~~~~~~~ Documentation ~~~~~~~~~~~~~~
The following documents are provided:
README.txt : This file that you are now reading
ExamplesForIndirectEstimation.pdf : A document providing the test examples
ExamplesForIndirectEstimation.doc : The original word version from which ExamplesForIndirectEstimation.pdf was created.
DeveloperGuide.pdf : Documents the estimation algorithm & the Python code
DeveloperGuide.doc : The original word version from which DeveloperGuide.pdf was created.
License.txt : The GPL license text

~~~~~~~~~~~~~~ Supporting Functions ~~~~~~~~~~~~~~
The following files contain functions used by the above scripts

MeanCalcPop.m : A function to calculate bounded mean of a population
PrevalenceCalcPop.m : A function to compute prevalence in a population
ParseMarkovTerm.m : A function that primarily parses a Markov term
PrepareForNumericSolution.m : A function preparing for a numeric solver
MyCDF.m : A function calculating a Cumulative distribution Function
MyInd2Sub.m : A function to replace Matlab’s Ind2Sub function
MySub2Ind.m : A function to replace Matlab’s MySub2Ind function
MyLog.m : A function to encapsulate maltabs log operation
MyDeal.m : A function to replace Matlab’s deal function 
MyInitMat.m : A function to replace Matlab’s Zeros and Ones
PrintFuncValueDiagnostics.m : A function to add printouts for a function
TableInfo.m : A function that parsed a table 
FindSubProcessStates.m : A function used to deal with sub-process data

~~~~~~~~~~~~~~ Sample inputs ~~~~~~~~~~~~~~
# The following files contain scripts to demonstrate software capabilities
TestScript.m : A test script that runs all the below examples
ProjectDataExample1.m : A script defining example 1
ProjectDataExample2.m : A script defining example 2
ProjectDataExample3.m : A script defining example 3
ProjectDataExample4.m : A script defining example 4
ProjectDataExample5.m : A script defining example 5
ProjectDataExample6.m : A script defining example 6
ProjectDataExample7.m : A script defining example 7
ProjectDataExample7_OnlyOneStudy.m  : A script, example 7 with only 1 study
ProjectDataExample8.m : A script defining example 8
ProjectDataExample9.m : A script defining example 9 - 2 years
ProjectDataExample9_OneYear.m : A script defining example 9 - 1 year
ProjectDataExample10_OneYear.m : A script defining example 10 - 2 year
ProjectDataExample10_TwoYears.m : A script defining example 10 - 3 years
ProjectDataExample10_ThreeYears.m : A script defining example 10 - 1 years
ProjectDataExample11_OneYear.m : A script defining example 11 - 1 year
ProjectDataExample11_TwoYears.m : A script defining example 11 - 2 years
ProjectDataExample11_ThreeYears.m : A script defining example 11 - 3 years
ProjectDataExample12_OneYear.m : A script defining example 12 - 1 year
ProjectDataExample12_TwoYears.m : A script defining example 12 - 2 year
ProjectDataExample12_ThreeYears.m : A script defining example 12 - 3 years
ProjectDataExample13.m : A script defining example 13
ProjectDataExample14.m : A script defining example 14
ProjectDataExample15.m : A script defining example 15
ProjectDataExample16.m : A script defining example 16
ProjectDataExample17.m : A script defining example 17
ProjectDataExample17_OneYear.m : A script defining example 17 – 1 year
ProjectDataExample18.m : A script defining example 18
ProjectDataExample18_OneYear.m : A script defining example 18 – 1 year
ProjectDataExample19.m : A script defining example 19
ProjectDataExample19_OneYear.m : A script defining example 19 – 1 year
ProjectDataExample20.m : A script defining example 20
ProjectDataExample20_OneYear.m : A script defining example 20 – 1 year
ProjectDataExample21.m : A script defining example 21
ProjectDataExample22.m : A script defining example 22
ProjectDataExample23.m : A script defining example 23
ProjectDataExample24.m : A script defining example 24
ProjectDataExample25.m : A script defining example 25
ProjectDataExample26.m : A script defining example 26
ProjectDataExample27.m : A script defining example 27

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



5. KNOWN BUGS AND FEATURES:

This prototype demonstrates the feasibility of the system. It has been 
tested using the provided example sets and final results were accurate 
to at least the one but last presented digit of the output for all the 
above examples.

The system contains disabled code that is considered future implementation
and still requires testing and may require modification. This code is 
currently disabled and will generate an error.

When using regression studies, model tables should not contain the same 
unknown coefficient in two different cells. It is the user’s responsibility
to take care of this, since the system may successfully generate an answer
in this case by using cell-wise optimization rather than coefficient-wise
optimization that are mathematically not the same.



6. VERSION HISTORY

Version 0.39 - 09-Dec-2009:
	- Now works on Matlab R2009a on a Linux machine
	- NaN prevalence is now ignored if Observed count difference is zero
	- Test suite has changed for examples 12,25
	- Documentation was updated

Version 0.38a - 29-Apr-2009:
    - No changes to code, only this readme file was updated for release
 
Version 0.38 - 3-Feb-2009:
    - 'Standard Error' was properly replaced with 'Variance'
    - Tables dimension range and sort order were corrected
    - likelihood is also kept in a vector of per-study likelihood
    - Typo correction in ExamplesForIndirectEstimation.pdf
    - DeveloperGuide.pdf added to fit the Python IEST software developer guide
    - TestScript.m now calculates/displays errors from the ideal results
    - The Lemonade Method was mentioned
 
Version 0.37:
    - Internal version never released see Version 0.38 for accumulated changes

Version 0.36 - 18-Dec-2007 :
    - First release




7. CONTACT INFORMATION:

Jacob Barhak Ph.D.
Application Systems Analyst / Programmer Senior
Department of Biostatistics
School of Public Health
University of Michigan
Room M4317D SPH II
1420 Washington Heights
Ann Arbor, Michigan 48109-2029
Phone: +1 (734) 647-5686
Fax: +1 (734) 763-2215
Email: jbarhak@umich.edu
http://www-personal.umich.edu/~jbarhak/


Wen Ye, Ph.D 
Research Assistant Professor
wye@umich.edu 
Phone: 734-615-9051
Fax: 734-763-2215 
Address:
Department of Biostatistics
School of Public Health
University of Michigan
M4537 SPH II 
1420 Washington Heights 
Ann Arbor, Michigan 48109-2029


8. ACKNOWLEDGEMENTS

Special thanks to the following (in no particular order) for their help in
various stages of this work:
Honghong Zhou
Morton B. Brown
William H. Herman
Fredric Isaman
Donghee Lee 
Shari Messinger Cayetano
Michael Brandle

Sponsors:
This work was supported by the National Institutes of Health through the 
Biostatistics Core of the Michigan Diabetes Research and Training Center 
under Grant (NIH: P60-DK20572), and through grant R21-DK075077 "Chronic 
Disease Modeling for Clinical Research Innovations" by the National 
Institutes of Health (EY03083).


