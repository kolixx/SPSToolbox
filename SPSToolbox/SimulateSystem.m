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
    
    Aorig = A;
    Borig = B;
    
    % Calculate the output of the system using hankel matrices. Not as
    % memory efficient but faster.
    
    % create a state space model for the system
    n = length(A)-1;
    ssA = zeros(n, n);
    ssB = zeros(n, 1);
    ssC = zeros(1, n);
    ssD = 0;
   
    % if the system has nonzero inpulse response at zero, do a polynomial
    % division
    if (B(1) ~= 0)
        [ssD, B] = deconv(B, A);
    end
    ssC = B(2:length(B));
    ssB = [1; zeros(n-1, 1)];
    ssA(1,:) = -A(2:length(A));
    for k=2:n
        ssA(k, k-1) = 1;
    end
    
    N = length(U);
    c = zeros(N,1);
    c(1) = ssD;
    if (n ~= 0)
        % this is not defined for static systems
        for k=2:N
            c(k) = ssC*ssA^(k-2)*ssB;
        end
    end
    H = toeplitz(c, zeros(1, N));
    Y = H*U;
    
%     -----------    
%     Calculate the output of the system using the standard sliding mode
%     simulation. Very slow.
%
%     A = Aorig;
%     B = Borig;
%     Y = zeros(length(U), 1);
%     % n-1, the length of the needed padding
%     n = length(A);
%     
%     Ubuffer = zeros(n, 1);
%     Ybuffer = zeros(n, 1);
%     
%     Atilde = A;
%     Atilde(1) = 0;
%     
%     for k=1:length(U)
%         % shift and fill the buffers
%         Ubuffer(2:n) = Ubuffer(1:(n-1));
%         Ubuffer(1) = U(k);
%         Ybuffer(2:n) = Ybuffer(1:(n-1));
%         
%         Y(k) = B*Ubuffer - Atilde*Ybuffer;
%         Ybuffer(1) = Y(k);
%     end
%     ------------

end




