function [contour] = Approximate2DConvexConfidenceRegionContour(thetaStart, sps, Y, X, membershipTestFunc, maxRadius)
%
%  Approximates the connected component around a given start point of the
%  confidence set. The algorithm relies on the assumption that the
%  component in starconvex around the given point.
%
%  Input arguments:
%  - thetaStart: a parameter vector that is assumed to be part of the
%  confidence region.
%  - sps: the data perturbation setup belonging to the confidence region
%  - membershipTestFunc: the test function of the model structure
%  - maxRadius: the maximum distance in the parameter space that should be
%  searched.
%
%  Output arguments:
%  - contour: a 2xK matrix with each column corresponding to a point on the
%  contour. The first and last points are the same closing the contour.
%

%  Copyright 2017 Sándor Kolumbán (s.kolumban@tue.nl)
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

    if (nargin < 6)
        maxRadius = inf;
    end

    if (length(thetaStart)~= 2)
        error('Only the contour of 2D sets can be approximated with this routine.');
    end
    
    if (~feval(membershipTestFunc, thetaStart, sps, Y, X))
        error('Initial point is not inside the confidence region');
    end
    
    ang = 0:0.01:2*pi;
    contour = zeros(2, length(ang));
    
    for k=1:length(ang)
        dir = rotPhi(ang(k))*[1; 0];
        [~, alphaU] = FindConfidenceRegionComponentBoundaryInDirection(thetaStart, sps, Y, X, dir, 10^-10, membershipTestFunc, maxRadius);
        contour(:, k) = thetaStart+alphaU*dir;
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
    backtrackAngle = 90;
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











