function [sps] = GenerateSPSSetup(q, m, N)
%
%  This routine generates an sps setup for confidence level 1-q/m for a set
%  of measurements consising og N samples.
% 
%  Output arguments
%  - sps.Signs: an N x m sign matrix. The first column is the all one
%  column.
%  - sps.TieOrder: to resolve possible tie situations, a random order is also
%  given.
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

    sps.Signs = round(rand(N, m-1))*2-1;
    sps.Signs = [ones(N, 1) sps.Signs];
    
    sps.TieOrder = randperm(m);
    
    sps.q = q;
    sps.m = m;
    sps.Confidence = 1-q/m;
end