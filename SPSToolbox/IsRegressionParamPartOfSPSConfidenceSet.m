function [inconf, Z] = IsRegressionParamPartOfSPSConfidenceSet(theta, sps, Y, X)
%
%  Gets if a given parameter vector theta is part of the sps confidence set or not.
%
%  The supposed model is that Y = X.' theta + error
%
%  Input parameters:
%    - theta: a column vector of size n
%    - sps: an sps setup correspnding to the same data length as X
%    - Y: column vector of size N
%    - X: an n x N matrix containing the regressor vectors as columns
%
%  Ouptut parameter:
%    - inconf: true if the given model was part of the confidence set
%     - Z: the different Z values that were used to determine the
%     membership

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
    NN(:,1) = Y - X.'*theta;
    
    % generate the sign perturbed noise samples
    for k=2:sps.m
        NN(:,k) = NN(:,1).*sps.Signs(:,k);
    end
    
    % calcualte Z_k value for the given ralization, and put it into a
    % vector
    Z = zeros(1, sps.m);
    
    [~,R] = qr(X.');
    Cov = R.'*R;
        
    for k=1:sps.m
        % sum l=1..n alpha it psi t Nt
        Sum = X*NN(:,k);
        Z(k) = Sum.'/Cov* Sum;
    end
    
    % get the order of Z(1)
    % the number of indeces that are bigger than the Z(1)
    r = DetermineRank(Z, sps.TieOrder);
    
    % is it larger or equal to q, return true
    inconf = r >= sps.q;
end