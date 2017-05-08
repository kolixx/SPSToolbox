function [alphaL, alphaU] = FindConfidenceRegionComponentBoundaryInDirection(theta, sps, Y, X, direction, tolerance, membershipTestFunc)
%
%
%
%
%

    if (~feval(membershipTestFunc, theta, sps, Y, X))
        error('theta not part of the confidence region. Boundary is searched only from within the confidence set.');
    end
    
    % apply binary search for the boundary point up to tolerance
    alphaL = 0;
    alphaU = 1;
    % find an upperbound in the direction
    while(feval(membershipTestFunc, theta+alphaU*direction, sps, Y, X))
        alphaU = 2*alphaU;
    end
    
    while (abs(alphaU-alphaL) > tolerance)
       % check the half point
       thetaHalf = theta+(alphaL + alphaU)/2*direction;
       if (feval(membershipTestFunc, thetaHalf, sps, Y, X))
           % inside the region, move the inside lower boundary to the half
           % point
           alphaL = (alphaL + alphaU)/2;
       else
           % half outside the set, move the external bound closer
           alphaU = (alphaL + alphaU)/2;
       end
    end
end