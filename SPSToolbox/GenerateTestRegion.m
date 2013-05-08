function [TR] = GenerateTestRegion(upperPoints, lowerPoints, stepsize)
%
%  Generates points inside a convex polygone with the given stepsize for
%  each direction.
%

%  Copyright 2013 Sándor Kolumbán (kolumban@aut.bme.hu)
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

    corners  = [upperPoints, lowerPoints];
    minx = min(corners(1, :));
    maxx = max(corners(1, :));
    
    miny = min(corners(2, :));
    maxy = max(corners(2, :));
    
    %x = minx:stepsize(1):maxx;
    %y = miny:stepsize(2):maxy;
    x = maxx:-stepsize(1):minx;
    y = maxy:-stepsize(2):miny;
    
    Z = zeros(length(x), length(y));
    
    for ix = 1:length(x)
        for iy=1:length(y)
            p = [x(ix); y(iy)];
            Z(ix, iy) = isPointInsidePolygone(p, upperPoints, lowerPoints);
        end
    end
    
    pos = 1;
    l = sum(sum(1-Z));
    TR = zeros(2, l);
    for ix = 1:length(x)
        for iy=1:length(y)
            if (Z(ix, iy) == 1)
                TR(:, pos) = [x(ix);y(iy)];
                pos = pos + 1;
            end
        end
    end
    
    TR = TR(:,1:(pos-1));
end

function [isInside] = isPointInsidePolygone(p, upperPoints, lowerPoints)

    isInside = 1;
    slope_old = inf;
    for segidx = 1:(size(upperPoints,2)-1)
        dir = upperPoints(:, segidx+1)-upperPoints(:, segidx);
        slope = dir(2)/dir(1);
        if (slope > slope_old)
            upperPoints(:, segidx)
            upperPoints(:, segidx+1)
            error('non convex domain, upper');
        end
        slope_old = slope;
        
        hatp = upperPoints(2, segidx)+slope*(p(1) - upperPoints(1, segidx));
        if (hatp < p(2))
           isInside = 0;
        end
    end
    
    slope_old = -inf;

    for segidx = 1:(size(lowerPoints,2)-1)
        dir = lowerPoints(:, segidx+1)-lowerPoints(:, segidx);
        slope = dir(2)/dir(1);
        
        if (slope < slope_old)
            lowerPoints(:, segidx)
            lowerPoints(:, segidx+1)
            error('non convex domain, lower');
        end
        slope_old = slope;
        hatp = lowerPoints(2, segidx)+slope*(p(1) - lowerPoints(1, segidx));
        if (hatp > p(2))
            isInside = 0;
        end
    end
end