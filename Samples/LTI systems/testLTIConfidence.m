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
%  The second part visualizes the results.
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

% the input signal will be one for the first K samples and -1 for the
% second K samples
K = 50;
U = ones(K,1);
U = [U; -U];

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


% the definition of the nominal model
model.A = [1];
model.F = [1 (0.9+0.8) 0.9*0.8];
model.B = [1];
model.C = [1 0];
model.D = [1 -0.1];

% run the trials
for testInstance = 1:numTest
    
    % write out the iteration counter to see some progress
    testInstance
    
    % generate a noise sequence for the trial
    N = (rand(length(U),1)-0.5)*0.1;
    % calculate the noisy output of the system
    [YN] = SimulateBJSample(model, U, N);

    % generate sps setups for the different confidence levels
    [sps95] = GenerateSPSSetup(1, 20, length(N));
    [sps90] = GenerateSPSSetup(1, 10, length(N));
    [sps85] = GenerateSPSSetup(3, 20, length(N));
    [sps80] = GenerateSPSSetup(1, 5, length(N));
    [sps75] = GenerateSPSSetup(1, 4, length(N));
    [sps70] = GenerateSPSSetup(3, 10, length(N));
    
    % ask the membership of the nominal model in the confidence regions
    % defined by the different sps setups
    indicators95(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps95, YN, U);
    indicators90(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps90, YN, U);
    indicators85(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps85, YN, U);
    indicators80(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps80, YN, U);
    indicators75(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps75, YN, U);
    indicators70(testInstance) = IsModelPartOfSPSConfidenceSet(model, sps70, YN, U);

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
        
