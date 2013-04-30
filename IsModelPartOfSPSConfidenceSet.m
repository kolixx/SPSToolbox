function [inconf] = IsModelPartOfSPSConfidenceSet(model, sps, Y, U)
%
%  Gets if a given model is part of the sps confidence set or not.
%

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
         
        [Q,R] = qr(Psi.');
        Cov = R.'*R;
        
        Z(k) = Sum.' *inv(Cov)* Sum;
    end
    
    % get the order of Z(1)
    % order the values of Z in descending orders
    [Zfirstround, idx] = sort(Z,'descend');
    
    % resolv ties around Z(1)
    z1pos = find(idx == 1, 1);
    tieIndeces = find(Zfirstround == Z(1));
    indecesToSort = idx(tieIndeces);
    [~, idxOrder] = sort(sps.TieOrder(indecesToSort));
    idx(indecesToSort) = indecesToSort(idxOrder);
    
    % the number of indeces that are bigger than the Z(1)
    r = z1pos - 1;
    
    % is it larger or equal to q, return true
    inconf = r >= sps.q;
end







