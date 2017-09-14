function [alphaL, alphaU] = FindConfidenceRegionComponentBoundaryInDirection(theta, sps, Y, X, direction, tolerance, membershipTestFunc, maxRadius)
%
%  Starting from an initial point theta the routine finds an upper and
%  lower approximation of the confidence region boundary in the given
%  directon, such that:
%
%  - theta+alphaL*direction is part of the confidence region
%  - theta+alphaU*direction is not part of the confidence region
%  - both alphaL and alphaU are positive
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

    if (nargin < 8)
        maxRadius = inf;
    end

    if (~feval(membershipTestFunc, theta, sps, Y, X))
        error('theta not part of the confidence region. Boundary is searched only from within the confidence set.');
    end
    
    % apply binary search for the boundary point up to tolerance
    alphaL = 0;
    alphaU = 1;
    % find an upperbound in the direction
    while(norm(alphaU*direction) < maxRadius^2 && feval(membershipTestFunc, theta+alphaU*direction, sps, Y, X))
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