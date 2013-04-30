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
    sps.Signs = round(rand(N, m-1))*2-1;
    sps.Signs = [ones(N, 1) sps.Signs];
    
    sps.TieOrder = randperm(m);
    
    sps.q = q;
    sps.m = m;
    sps.Confidence = 1-q/m;
end