function [r] = DetermineRank(Z, perm)
%
%  Gets the number of values in the vector Z that are greater than Z(1).
%  Ties are solved according to the given permutation.
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
%  along with SPSToolbox.  If not, see <http://www.gnu.org/licenses/>.

    % get the order of Z(1)
    % order the values of Z in descending orders
    [Zfirstround, idx] = sort(Z,'descend');
    
    % resolv ties around Z(1)
    z1pos = find(idx == 1, 1);
    % where are the tie positions after the first round in Zfirstround
    tieIndeces = find(Zfirstround == Z(1));
    % the indeces of the Z values that are on these tie positions
    indecesToSort = idx(tieIndeces);
    [~, idxOrder] = sort(perm(indecesToSort), 'descend');
    idx(tieIndeces) = idx(tieIndeces(idxOrder));
    z1pos = find(idx == 1, 1);
    
    % the number of indeces that are bigger than the Z(1)
    r = z1pos - 1;

end