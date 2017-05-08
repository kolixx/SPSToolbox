function [Z] = GetSPSdefinedPerformanceForBJModelNU(model, N, U)
        
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