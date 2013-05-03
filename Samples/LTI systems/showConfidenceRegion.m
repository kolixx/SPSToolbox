%
%  Explores a confidence set in two dimensions and visualizes it.
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

% clear the console
clc
close all

% aply a step and a negative step as input
K = 25;
U = ones(K,1);
U = [U; -U];

% generate a nominal model with two unknown paramters to be estimated (f_1
% and b_1)
model.A = [1];
model.F = [1 -0.6];
model.B = [1];
model.C = [1];
model.D = [1];

% generate an SPS setup for 0.85 confidence level
% - the saved version is loaded: sps85 = GenerateSPSSetup(3, 20, length(U));

% generate a noise sequence
% - the saved version is loaded: N = randn(length(U),1)*0.03;
%save('confRegionRandomness.mat', 'N', 'sps85')
load('confRegionRandomness.mat', 'N', 'sps85')

% simulate the system
[YN] = SimulateBJSample(model, U, N);

% the noisy output
figure(1), plot(YN)

data = iddata(YN, U, 1);
sys = bj(data, [1 0 0 1 1], 'InitialCondition','zero');

 % create the pem model
modelPEM.A = sys.A;
modelPEM.F = sys.F;
modelPEM.B = sys.B(2);
modelPEM.C = sys.C;
modelPEM.D = sys.D;

display('Is the PEM estimate is inside the confidence region (should be 1)?')
IsModelPartOfSPSConfidenceSet(modelPEM, sps85, YN, U)

%% generate points inside a polygone
upperPoints = [-0.602 -0.60 -0.593 -0.588;...
               1      1.01  1.026  1.03 ];
lowerPoints = [-0.602 -0.598 -0.59 -0.588;...
                1      1     1.018 1.03];
stepsize = [0.0001; 0.0005];



[TR] = GenerateTestRegion(upperPoints, lowerPoints, stepsize);
% visualize the test region
figure(2), plot(TR(1,:), TR(2,:), '.', lowerPoints(1,:), lowerPoints(2,:), '-', upperPoints(1,:), upperPoints(2,:), '-');



%% test each individual model if it is inside the confidence set or not
Fs = TR(1,:);
Bs = TR(2,:);
nconf = length(Fs);

Z = zeros(1, nconf);
percent = 0;
display('Percantage of evaluated models on the grid:')
for k=1:nconf
    if (k/nconf*100 >= percent)
        display(sprintf('%d %', percent));
        percent = percent + 10;
    end
    mkl.A = [1];
    mkl.F = [1 Fs(k)];
    mkl.B = [Bs(k)];
    mkl.C = [1];
    mkl.D = [1];
   
    Z(k) = IsModelPartOfSPSConfidenceSet(mkl, sps85, YN, U);
end

%% save those models that belong to the confidence set
FB = zeros(2, sum(sum(Z)));
pos = 1;
for k=1:nconf
    if (Z(k) == 1)
        FB(:, pos) = [Fs(k); Bs(k)];
        pos = pos + 1;
    end
end

figure(3) 
plot(TR(1,:), TR(2,:), 'b.', ... % plot the investigated region
    lowerPoints(1,:), lowerPoints(2,:), 'b-', upperPoints(1,:), upperPoints(2,:), 'b-', ... % plot the region boundaries
    FB(1,:), FB(2,:), 'r.', ... % plot the confidence set
    modelPEM.F(2), modelPEM.B(1),'gx'... % plot the pem estimate
    );
