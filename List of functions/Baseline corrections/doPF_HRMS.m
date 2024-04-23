%% Description
% Polynomial fitting
% Unpublished results
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, BslPts, p] = doPF_HRMS(XY, k, n)
maxIter = 2;
maxIterO = 10;
cXY = XY;
cXY(isnan(XY(:,2)), :) = [];
p0 = zeros(n, 1);
options = optimset('Display','none');
[p,fval,exitflag] = fminsearch(@objectivefcn1, p0, options);
count = 1;

while exitflag == 0  & count < maxIter
    [p,fval,exitflag] = fminsearch(@objectivefcn1, p, options);
    count = count +1;

end

z = polyval(p, XY(:,1));

% Center 
XY(:, 3) = XY(:, 2) - z; XY(:, 4) = 1;
count = 1;
while 1
    L = (1:size(XY, 1))';  
    L(XY(:, 4) == 0) = [];
    IO = isoutlier(XY(XY(:, 4) ==1, 3));
    if ~any(IO), break; end
    XY(L(IO), 4) = 0; 
    count = count +1;if count > maxIterO, break; end
end
z = z +median(XY(XY(:, 4)==1, 3));
BslPts = XY(:, 4) == 0;

    function [f, bsl] = objectivefcn1(a)
        bsl = polyval(a, cXY(:,1));
        prof = cXY(:, 2)-bsl;
        id = prof < 0;
        prof(id) = prof(id).^k;
        f = ChrMoment([cXY(:, 1), prof], 1);
    end
end