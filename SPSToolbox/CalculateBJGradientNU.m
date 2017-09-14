function [Psi] = CalculateBJGradientNU(model, N, U)
%
%   Calculates the gradient of the prediction errors or noise ralizaions
%   from the output and input and model data. These column vectors are
%   stacked into a nparam x nsample matrix.
%
%   The order of the first dimension is the following: a1, ..., b0, ...,
%   f1, ..., c1, ...., d1, ....
%
%   Input parameters: 
%      - model: the model structure describing the model
%      - N, U: the noise realization entering the noise model and the input
%      of the nominal model colum vectors (dimension N x 1).
%
%   Output parameters:
%      - Psi: an n_theta times N dimension matrix containing the
%      derivatives of the prediction errors with respect to the model
%      parameters as columns.
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

    nparam = length(model.A)-1 +...
             length(model.B) + ...
             length(model.F)-1 + ...
             length(model.C)-1 + ...
             length(model.D)-1;
         
    Psi = zeros(nparam, length(N));
    
    % the helper variable is D/(CF)u(t-k), when it is calculated it should
    % memorized
    % this variable shows the bigest index for which the helper variable is
    % calculated already
    helperVariableMaxK = 0;
    % the helper variable might be needed in the calculation of the A and F
    % coefficient derivatives
    % reserve space for its values depending on the highest degree
    helperVariable = zeros(length(N), max(length(model.A)-1, length(model.F)-1));
    
    % the last used row of the Psi matrix
    usedindeces = 0;
    
    % add the derivatives related to A
    % A is monic, start only at the second position if there is any
    for k=1:(length(model.A)-1)
        % pad the U with k zeros, and cut k terms at the end
        Upadded = [zeros(k, 1); U(1:length(U)-k)];
        % pad the N with k zeros, and cut k terms at the end
        Npadded = [zeros(k, 1); N(1:length(N)-k)];
        
        % calculate the helper variable
        CF = conv(model.C, model.F);
        D = [model.D zeros(1, length(CF)-length(model.D))];
        hv = SimulateSystem(D, CF, Upadded);
        helperVariableMaxK = max(helperVariableMaxK, k);
        helperVariable(:, k) = hv;
        
        B = model.B;
        A = model.A;
        [A,B] = paddatend(A,B);
        One = 1;
        One = [One zeros(1, length(A)-length(One))];
        Psi(usedindeces+1,:) = SimulateSystem(B, A, hv).' ...
                   + ...
                   SimulateSystem(One, A, Npadded).';
       
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to B
    % B is not monic, start only at the first position if there is any
    for k=((0:(length(model.B)-1)) + length(model.F)- length(model.B))
        
        hv = [];
        
        if (k<=helperVariableMaxK)
            % if the helper variable is calculated, use it
            hv = helperVariable(:, k);
        else
            % pad the U with k zeros, and cut k terms at the end
            Upadded = [zeros(k, 1); U(1:length(N)-k)];
        
            CF = conv(model.C, model.F);
            D = [model.D zeros(1, length(CF)-length(model.D))];
            hv = SimulateSystem(D, CF, Upadded);
            helperVariableMaxK = max(helperVariableMaxK, k);
            helperVariable(:, k) = hv;
        end
        
        Psi(usedindeces+1,:) = hv.';
        
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to C
    % C is monic, start only at the second position if there is any
    for k=1:(length(model.C)-1)
        % pad the N with k zeros, and cut k terms at the end
        Npadded = [zeros(k, 1); N(1:length(N)-k)];
        
        C = model.C;
        One = 1;
        One = [One zeros(1, length(C)-length(One))];
                
        Psi(usedindeces+1,:) = -SimulateSystem(One, C, Npadded).';
        
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to D
    % D is monic, start only at the second position if there is any
    for k=1:(length(model.D)-1)
        % pad the N with k zeros, and cut k terms at the end
        Npadded = [zeros(k, 1); N(1:length(N)-k)];
        
        D = model.D;
        One = 1;
        One = [One zeros(1, length(D)-length(One))];
                
        Psi(usedindeces+1,:) = SimulateSystem(One, D, Npadded).';
        
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to F
    % F is monic, start only at the second position if there is any
    for k=1:(length(model.F)-1)
                
        hv = [];
        
        if (k<=helperVariableMaxK)
            % if the helper variable is calculated, use it
            hv = helperVariable(:, k);
        else
            % pad the U with k zeros, and cut k terms at the end
            Upadded = [zeros(k, 1); U(1:length(N)-k)];
        
            CF = conv(model.C, model.F);
            D = [model.D zeros(1, length(CF)-length(model.D))];
            hv = SimulateSystem(D, CF, Upadded);
            helperVariableMaxK = max(helperVariableMaxK, k);
            helperVariable(:, k) = hv;
        end
        
        B = model.B;
        F = model.F;
        [B,F] = paddatend(B,F);
        Psi(usedindeces+1,:) = SimulateSystem(B, F, hv).';
        
        usedindeces = usedindeces + 1;
    end
end

function [A,B] = paddatend(A,B)
    if (length(A) < length(B))
        A = [A zeros(1, length(B)-length(A))];
    else
        B = [B zeros(1, length(A)-length(B))];
    end
end





