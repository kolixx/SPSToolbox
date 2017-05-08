function [contour] = Approximate2DConfidenceRegionContour(thetaStart, sps, Y, X, membershipTestFunc)
%
%
%
%
    if (length(thetaStart)~= 2)
        error('Only the contour of 2D sets can be approximated with this routine.');
    end
    
    if (~feval(membershipTestFunc, thetaStart, sps, Y, X))
        error('Initial point is not inside the confidence region');
    end
    
    % find the boundary of the confidence region in a random direction
    dir = randn(2,1);
    dir = [1; 0];
    dir = dir/norm(dir);
    
    [alphaL1, alphaU1] = FindConfidenceRegionComponentBoundaryInDirection(thetaStart, sps, Y, X, dir, 10^-10, membershipTestFunc);
    [~, alphaU2] = FindConfidenceRegionComponentBoundaryInDirection(thetaStart, sps, Y, X, -dir, 10^-10, membershipTestFunc);
    
    % the contour is going to be of resolution diameter/100 in a random
    % direction
    stepsize = (alphaU1+alphaU2)/100;
    
    % start the contour from an outside point next to the boundary
    contour = thetaStart+alphaU1*dir;
    % set the initial direction perpendicular to the direction from the
    % central start point
    % such that the outside of the interval is to the left
    dir = rotPhi(-pi/2)*dir;
    
    diff = [];
    
    while ((size(contour, 2) < 100 || norm(contour(:,1)-contour(:,end)) > 10*stepsize) && size(contour,2) < 8000)
        diff  = [diff norm(contour(:,1)-contour(:,end))];
        [newpoint, dir] = findNextPoint(contour(:,end), dir, stepsize, sps, Y, X, membershipTestFunc);
        contour = [contour newpoint];
    end
    
    %close the region
    contour = [contour contour(:,1)];
end

function [newpoint, newdir] = findNextPoint(theta, direction, stepsize, sps, Y, X, membershipTestFunc)
%
%  Tracing the contour, numerically problematic.
%

    thetaX = [ 1.005555523856740;    0.839053506328663];
    if (theta(1) < 1.0055)
    end
   
    % find the angle of the next step using binary search
    backtrackAngle = 10;
    angU = pi-backtrackAngle/360*2*pi;
    angL = -angU;
    
    if (feval(membershipTestFunc, theta+rotPhi(angU)*direction*stepsize, sps, Y, X) ...
        ||...
        ~feval(membershipTestFunc, theta+rotPhi(angL)*direction*stepsize, sps, Y, X))
        error('not proper directionality')
    end
    
    while (angU-angL > 10^-5)
        
        midAng = (angU+angL)/2;
        
        if (feval(membershipTestFunc, theta+rotPhi(midAng)*direction*stepsize, sps, Y, X))
            angL = midAng;
        else
            angU = midAng;
        end
    end
    
    if (angU > pi/4)
    end

    % perform the step, for an outer point
    newpoint = theta+rotPhi(angU)*direction*stepsize;
    
    newdir = newpoint-theta;
    newdir = newdir/norm(newdir);
end

function [R] = rotPhi(phi)

    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];

end











