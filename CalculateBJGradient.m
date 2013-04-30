function [Psi] = CalculateBJGradient(model, Y, U)
%
%   Calculates the gradient of the prediction errors or noise ralizaions
%   from the output and input and model data. These column vectors are
%   stacked into a nparam x nsample matrix.
%
%   The order of the first dimension is the following: a1, ..., b0, ...,
%   f1, ..., c1, ...., d1, ....
%

    nparam = length(model.A)-1 +...
             length(model.B) + ...
             length(model.F)-1 + ...
             length(model.C)-1 + ...
             length(model.D)-1;
         
    Psi = zeros(nparam, length(Y));
    
    % the last used row of the Psi matrix
    usedindeces = 0;
    
    % add the derivatives related to A
    % A is monic, start only at the second position if there is any
    for k=1:(length(model.A)-1)
        % pad the Y with k zeros, and cut k terms at the end
        Ypadded = [zeros(k, 1); Y(1:length(Y)-k)];
        Psi(usedindeces+1,:) = SimulateSystem(model.D, model.C, Ypadded).';
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to B
    % B is not monic, start only at the first position if there is any
    for k=((0:(length(model.B)-1)) + length(model.F)- length(model.B))
        % pad the U with k zeros, and cut k terms at the end
        Upadded = [zeros(k, 1); U(1:length(Y)-k)];
        
        CF = conv(model.C, model.F);
        D = [model.D zeros(1, length(CF)-length(model.D))];
        
        Psi(usedindeces+1,:) = SimulateSystem(D, CF, Upadded).';
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to C
    % C is monic, start only at the second position if there is any
    for k=1:(length(model.C)-1)
        % pad the Y with k zeros, and cut k terms at the end
        Ypadded = [zeros(k, 1); Y(1:length(Y)-k)];
        % pad the U with k zeros, and cut k terms at the end
        Upadded = [zeros(k, 1); U(1:length(Y)-k)];
        
        CCF = conv(conv(model.C, model.C), model.F);
        DB = conv(model.D, model.B);
        DB = [DB zeros(1, length(CCF)-length(DB))];
        
        CC = conv(model.C, model.C);
        DA = conv(model.D, model.A);
        [DA,CC] = paddatend(DA,CC);
        
        Psi(usedindeces+1,:) = SimulateSystem(DB, CCF, Upadded).'...
            - ...
            SimulateSystem(DA, CC, Ypadded).';
        
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to D
    % D is monic, start only at the second position if there is any
    for k=1:(length(model.D)-1)
        % pad the Y with k zeros, and cut k terms at the end
        Ypadded = [zeros(k, 1); Y(1:length(Y)-k)];
        % pad the U with k zeros, and cut k terms at the end
        Upadded = [zeros(k, 1); U(1:length(Y)-k)];
        
        A = model.A;
        C = model.C;
        [A,C] = paddatend(A,C);
        
        B = model.B;
        CF = conv(model.C, model.F);
        [B,CF] = paddatend(B,CF);
        
        Psi(usedindeces+1,:) = - SimulateSystem(B, CF, Upadded).'...
            + ...
            SimulateSystem(A, C, Ypadded).';
        
        usedindeces = usedindeces + 1;
    end
    
    % add the derivatives related to F
    % F is monic, start only at the second position if there is any
    for k=1:(length(model.F)-1)
        % pad the U with k zeros, and cut k terms at the end
        Upadded = [zeros(k, 1); U(1:length(Y)-k)];
        
        DB = conv(model.D, model.B);
        CFF = conv(conv(model.C, model.F), model.F);
        [DB,CFF] = paddatend(DB,CFF);
        Psi(usedindeces+1,:) = SimulateSystem(DB, CFF, Upadded).';
        
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




