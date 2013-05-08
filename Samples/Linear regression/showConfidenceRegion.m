%
%  Explores a confidence set in two dimensions and visualizes it.
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


% clear the console
clc
close all

% the nominal model
theta= [1;...
        2];
% generate a random regressor matrix
% the number of samples
N = 100;
% the dimension of theta
n = length(theta);
% the random regressors
% - the saved version is loaded: X = rand(n, N);

% generate an SPS setup for 0.85 confidence level
% - the saved version is loaded: sps85 = GenerateSPSSetup(3, 20, N);

% generate a noise sequence
% - the saved version is loaded: E = randn(N,1)*0.03;
%save('confRegionRandomness.mat', 'E', 'sps85', 'X')
load('confRegionRandomness.mat', 'E', 'sps85', 'X')

% simulate the system
[YN] = X.'*theta+E;

thetaLS = (X*X.')\X*YN;

display('Is the LS estimate is inside the confidence region (should be 1)?')
IsRegressionParamPartOfSPSConfidenceSet(thetaLS, sps85, YN, X)

%% generate points inside a polygone
upperPoints = [0.985 0.985 0.994 1.01;...
               2.014 2.015  2.016  2];
lowerPoints = [0.985 0.985 1     1.01 1.01; ...
               2.014 2.005 1.99  1.985 2];
stepsize = [0.0005; 0.0005];



[TR] = GenerateTestRegion(upperPoints, lowerPoints, stepsize);
% visualize the test region
figure(2), plot(TR(1,:), TR(2,:), '.', lowerPoints(1,:), lowerPoints(2,:), '-', upperPoints(1,:), upperPoints(2,:), '-');



%% test each individual model if it is inside the confidence set or not
nconf = size(TR,2);

Z = zeros(1, nconf);
percent = 0;
display('Percantage of evaluated models on the grid:')
for k=1:nconf
    if (k/nconf*100 >= percent)
        display(sprintf('%d %', percent));
        percent = percent + 10;
    end
    thetak = TR(:, k);
   
    Z(k) = IsRegressionParamPartOfSPSConfidenceSet(thetak, sps85, YN, X);
end

%% save those models that belong to the confidence set
FB = zeros(2, sum(sum(Z)));
pos = 1;
for k=1:nconf
    if (Z(k) == 1)
        FB(:, pos) = TR(:, k);
        pos = pos + 1;
    end
end

figure(3) 
plot(TR(1,:), TR(2,:), 'b.', ... % plot the investigated region
    lowerPoints(1,:), lowerPoints(2,:), 'b-', upperPoints(1,:), upperPoints(2,:), 'b-', ... % plot the region boundaries
    FB(1,:), FB(2,:), 'r.', ... % plot the confidence set
    thetaLS(1), thetaLS(2),'gx'... % plot the pem estimate
    );
