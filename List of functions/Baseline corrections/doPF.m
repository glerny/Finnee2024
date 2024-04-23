%% Description
% Polynomial fitting
% Unpublished results
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, w] = doPF(XY, n)
% input: XY 2xn array,
%        n order iof the polynomial
% output:
%        z basline model
%        bslPts false if peak pts, true in baseline pts
%% CHANGES 05/02/2024
N = size(XY, 1);
w = true(N, 1);
iterMax = 10;
L = N;

for ii = 1:iterMax
    
    p = polyfit(XY(w, 1), XY(w, 2), n);
    z = polyval(p, XY(:,1));
    w = ~isoutlier(XY(:,2) - z);
    
    %%% ORIGINAL
    
    if sum(w) ==  L, break; end
    if sum(w) <= 5, break; end
    
    L = sum(w);
    % Check the baseline points and stop the loop if no more points are
    % detected as baseline
    
end

