function [Y] = SimulateSystem(B, A, U)
%
%   Simulates a B/A linear system with input U from initial starting
%   condition.
%
    if (length(B)~=length(A))
        error('Nominator and denominator should have the same length');
    end
    
    if (A(1)~= 1)
        error('The denominator should have leading coefficient one')
    end
    
    Y = zeros(length(U), 1);
    
    % n-1, the length of the needed padding
    n = length(A);
    
    Ubuffer = zeros(n, 1);
    Ybuffer = zeros(n, 1);
    
    Atilde = A;
    Atilde(1) = 0;
    
    for k=1:length(U)
        % shift and fill the buffers
        Ubuffer(2:n) = Ubuffer(1:(n-1));
        Ubuffer(1) = U(k);
        Ybuffer(2:n) = Ybuffer(1:(n-1));
        
        Y(k) = B*Ubuffer - Atilde*Ybuffer;
        Ybuffer(1) = Y(k);
    end
end