%
%  This script checks that the prediction error minimizer estimate does
%  always belong to the confidence set.
%

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






