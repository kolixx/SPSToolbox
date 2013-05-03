function [inconf] = IsModelPartOfSPSConfidenceSet(model, sps, Y, U)
%
%  Gets if a given model is part of the sps confidence set or not.
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

    % the matrix containing the sign perturbed noise ralizatons
    NN = zeros(length(Y), sps.m);
    
    % get the noise realization corresponding to the given model
    NN(:,1) = CalculateBJNoiseRealization(model, Y, U);
    
    % generate the sign perturbed noise samples
    for k=2:sps.m
        NN(:,k) = NN(:,1).*sps.Signs(:,k);
    end
    
    % get the perturbed outputs
    YY = zeros(length(Y), sps.m);
    YY(:,1) = Y;
    % this should be the nominal response
    Y1 = SimulateBJSample(model, U, NN(:,1));
    for k=2:sps.m
        YY(:,k) = SimulateBJSample(model, U, NN(:,k));
    end
    
    % calculate the psi matrix for each realization
    % calcualte Z_k value for the given ralization, and put it into a
    % vector
    Z = zeros(1, sps.m);
    
    for k=1:sps.m
        [Psi] = CalculateBJGradient(model, YY(:,k), U);
        
        % sum l=1..n alpha it psi t Nt
        Sum = Psi*NN(:,k);
         
        [~,R] = qr(Psi.');
        Cov = R.'*R;
        
        Z(k) = Sum.'/Cov* Sum;
    end
    
    % get the order of Z(1)
    % the number of indeces that are bigger than the Z(1)
    r = DetermineRank(Z, sps.TieOrder);
    
    % is it larger or equal to q, return true
    inconf = r >= sps.q;
end







