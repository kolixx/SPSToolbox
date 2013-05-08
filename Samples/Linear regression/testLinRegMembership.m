%
%  This script checks that the prediction error minimizer estimate does
%  always belong to the confidence set.
%

%  Copyright 2013 Sándor Kolumbán (kolumban@aut.bme.hu)
%
%  The program is distributed under the terms of the GNU General Public License.
%
%  This file is part of SPSToolbox
%
%  SPSToolbox is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  SPSToolbox is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with SPSToolbox.  If not, see <http://www.gnu.org/licenses/>.

clc

% clear the console and the workspace
clc

% the nominal model
theta= [1;...
        0.3;...
        2];
% generate a random regressor matrix
% the number of samples
N = 100;
% the dimension of theta
n = length(theta);
% the random regressors
X = rand(n, N);
% the SPS setup
[sps85] = GenerateSPSSetup(3, 20, N);

% generate a noise sequence for the trial
E = (rand(N,1)-0.5)*0.1;

% calculate the noisy output of the system
[YN] = X.'*theta + E;

% get the LS estimate
thetaLS = (X*X.')\X*YN;

display('Result of membership test (should be 1):')
IsRegressionParamPartOfSPSConfidenceSet(thetaLS, sps85, YN, X)






