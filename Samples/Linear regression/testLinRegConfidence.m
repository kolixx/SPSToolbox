%
% This script is an example that illustrates the correct behaviour of the
% SPS hypothesis test
%
%
% Confidence set membership query is perfomed for a large number of
% realizations. For a confidece set of confidence 95%, the membership of
% the nominal model should be true for 95% of the trials. Different
% condidence levels are checked.
%
% The fist part of the script is time consuming, it calculates a lot of
% indicator values. See more detailed description inside the script.
%
% The second part visualizes the results.

%
%  Copyright 2013 Sándor Kolumbán (kolumban@aut.bme.hu)
%
%  The program is distributed under the terms of the GNU General Public License.
%
%  This file is part of SPSToolbox
%
%  SPSToolboxis free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  Foobar is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

%% Claculation of indicator values.

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

% the number of independ trials for different confidence sets
numTest = 15000;

% these arrays hold 1 or 0 depending wheter the nominal model belongs to
% the corresponding confidence region or not
indicators95 = zeros(1, numTest);
indicators90 = zeros(1, numTest);
indicators85 = zeros(1, numTest);
indicators80 = zeros(1, numTest);
indicators75 = zeros(1, numTest);
indicators70 = zeros(1, numTest);

% run the trials
for testInstance = 1:numTest
    
    % write out the iteration counter to see some progress
    testInstance
    
    % generate a noise sequence for the trial
    E = (rand(N,1)-0.5)*0.1;
    % calculate the noisy output of the system
    YN = X.'*theta + E;

    % generate sps setups for the different confidence levels
    [sps95] = GenerateSPSSetup(1, 20, N);
    [sps90] = GenerateSPSSetup(1, 10, N);
    [sps85] = GenerateSPSSetup(3, 20, N);
    [sps80] = GenerateSPSSetup(1, 5, N);
    [sps75] = GenerateSPSSetup(1, 4, N);
    [sps70] = GenerateSPSSetup(3, 10, N);
    
    % ask the membership of the nominal model in the confidence regions
    % defined by the different sps setups
    indicators95(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps95, YN, X);
    indicators90(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps90, YN, X);
    indicators85(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps85, YN, X);
    indicators80(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps80, YN, X);
    indicators75(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps75, YN, X);
    indicators70(testInstance) = IsRegressionParamPartOfSPSConfidenceSet(theta, sps70, YN, X);

end

save workspace

%% Visualize the results. 

load workspace

idx = [0.7 0.75 0.8 0.85 0.9 0.95];
vals = [mean(indicators70) ...
        mean(indicators75) ...
        mean(indicators80) ...
        mean(indicators85) ...
        mean(indicators90) ...
        mean(indicators95)];
    
variances = [var(indicators70) ...
            var(indicators75) ...
            var(indicators80) ...
            var(indicators85) ...
            var(indicators90) ...
            var(indicators95)]; 
%What is expected is that the red signs are exactly on top of the blue ones
close all
hold on
plot(idx, idx,'bo', idx, vals, 'rx')
xlim([0.65, 1])
ylim([0.65, 1])
hold off
        
