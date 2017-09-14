function [N] = CalculateBJNoiseRealization(model, Y, U)
%
%  Based on the given box jenkins model, the measured outputs and the
%  inputs, the routine returns the noise realization used to generate the
%  measurements, assuming zero initial conditions.
%
%  Input arguments:
%    - model: a struct containing the coefficients A,B,C,D,F of a
%    BoxJenkings model.
%    - Y: column vector of length K containing the outputs N[1] ... N[K]
%    - U: column vector of length K containing the inputs U[1] ... U[K]
%
%  Output arguments:
%    - N: column vector of length K containing the noise samples N[1] ... N[K]

%  Copyright 2013 S�ndor Kolumb�n (s.kolumban@tue.nl)
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

    % check the input
    % deg: AF > B, AD > C
    % monic: A,F,D,C
    if (length(model.A)+length(model.F) < length(model.B))
        error('degree error');
    end
    if (length(model.A)+length(model.D) < length(model.C))
        error('degree error');
    end
    if (model.A(1) ~= 1 || model.F(1) ~= 1 || model.D(1) ~= 1 || model.C(1) ~= 1)
        error('monic error')
    end
    
    [un,um] = size(U);
    if (um~=1)
        error('Dimension issue');
    end
    [yn,ym] = size(Y);
    if (ym~=1 || un ~= yn)
        error('Dimension issue');
    end
    
    % yt = B/(AF) U + C/(AD) N = Y1 + Y2
    
    AF = conv(model.A, model.F);
    AD = conv(model.A,model.D);
    
    B = model.B;
    B = [zeros(1, length(AF)-length(B)) B];
    C = model.C;
    if (length(C)~= length(AD))
        % assumption 2
        error('Noise model not properly invertible');
    end
    if (B(1)~= 0)
        % assumption 2
        error('The first element of the inpulse response should be zero.')
    end
    if (max(abs(roots(C))) >=1)
        % assumption 2
        error('No stable inverse for noise model');
    end
    
    % the predicted output wothout noise
    Y1 = SimulateSystem(B, AF, U);
    
    % the output of the noise system
    No = Y-Y1;
    N = SimulateSystem(AD, C, No);
end