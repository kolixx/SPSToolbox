function [Z] = GetSPSdefinedPerformanceForBJModelNU(model, N, U)
%
%   Calculates the performance measure corresponding to Box-Jenkins models
%   defined by the [1], that is the norm of the gradient of the
%   least-squares cost function.
%
%   [1] Csáji, B. Cs.; Campi, M. C.; Weyer, E.: Sign-Perturbed Sums (SPS): A Method for Constructing Exact Finite-Sample Confidence Regions for General Linear Systems, 51st IEEE Conference on Decision and Control (CDC 2012), Maui, Hawaii, 2012, pp. 7321–7326
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
%  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

        
    %------ Calculate the gradinets using Y and U
    %       Calculation based on the noise samples is faster  
    %------
    %[Psi] = CalculateBJGradientYU(model, Y, U);
    %------
    [Psi] = CalculateBJGradientNU(model, N, U);

    % sum l=1..n alpha it psi t Nt
    Sum = Psi*N;

    [~,R] = qr(Psi.');
    Cov = R.'*R;

    Z = Sum.'/Cov* Sum;
end