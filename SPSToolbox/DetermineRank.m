function [r] = DetermineRank(Z, perm)
%
%  Gets the number of values in the vector Z that are greater than Z(1).
%  Ties are solved according to the given permutation.
%

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