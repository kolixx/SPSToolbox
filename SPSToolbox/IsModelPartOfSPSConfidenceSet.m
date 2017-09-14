function [inconf, Z] = IsModelPartOfSPSConfidenceSet(model, sps, Y, U, performanceMeasure)
%
%  Gets if a given model is part of the sps confidence set or not.
%
%  Input arguments:
%     - the model that describes the system structure
%     - sps: the sps setup
%     - Y, U: column vectors of the same size containing the output and the
%     input respectively
%     - performanceMeasure: a function handle to evaluate DP peformance
%
%
%  Output arguments:
%     - inconf: boolean indicator of the membership to the confidence set
%     - Z: the different Z values that were used to determine the
%     membership
%
%

%  Copyright 2013 Sándor Kolumbán (s.kolumban@tue.nl)
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

    % the matrix containing the sign perturbed noise ralizatons
    NN = zeros(length(Y), sps.m);
    
    % get the noise realization corresponding to the given model
    NN(:,1) = CalculateBJNoiseRealization(model, Y, U);
    
    % generate the sign perturbed noise samples
    for k=2:sps.m
        NN(:,k) = NN(:,1).*sps.Signs(:,k);
    end
    
%     % get the perturbed outputs
%     ------ this is only needed, if the gradients are calculated based on
%            Y and U (and to give a more didactic presentation of the SPS method)
%            but the values are only needed if the gradients are calculated
%            using CalculateBJGradientYU
%     ------
%     YY = zeros(length(Y), sps.m);
%     YY(:,1) = Y;
%     % this should be the nominal response
%     Y1 = SimulateBJSample(model, U, NN(:,1));
%     for k=2:sps.m
%         YY(:,k) = SimulateBJSample(model, U, NN(:,k));
%     end
%     ------
    
    % calcualte Z_k value for the given ralization, and put it into a
    % vector
    Z = zeros(1, sps.m);
    
    for k=1:sps.m
        %GetSPSdefinedPerformance(model,  NN(:,k), U);
        Z(k) = feval(performanceMeasure, model,  NN(:,k), U);
    end
    
    % get the order of Z(1)
    % the number of indeces that are bigger than the Z(1)
    r = DetermineRank(Z, sps.TieOrder);
    
    % is it larger or equal to q, return true
    inconf = r >= sps.q;
end







