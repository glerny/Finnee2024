function [mzAxis, nXY, stats] = mkMasterMZ(MS, split, n, Lim, ratio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    Lim = [min(MS(:, 1)), max(MS(:, 1))];
    split = 3;
    n   = 3;
    ratio = 1;

elseif nargin == 2

    Lim = [min(MS(:, 1)), max(MS(:, 1))];
    n   = 3;
    ratio = 1;

elseif nargin == 3
    Lim = [min(MS(:, 1)), max(MS(:, 1))];
    ratio = 1;

elseif nargin == 4
    ratio = 1;
end

W2C = zeros([size(MS, 1)+2 5]);
W2C(2:end-1, 1) = MS(:,1);
W2C(1:end-2, 2) = MS(:,1);
W2C(2:end-1, 3) = MS(:,2);
W2C(:,4) = circshift(W2C(:,3), -1);
W2C(:,5) = circshift(W2C(:,3), +1);
W2C(W2C(:,3) == 0 | W2C(:,4) == 0 | W2C(:,5) == 0, :) = [];
XY = [W2C(:,1), W2C(:,2)-W2C(:,1)];
MStep = round(height(XY)/split);
nXY = [];

myCut(1, :) = [1, MStep];
while myCut(end, 2) < height(XY)
    myCut(end+1, :) = [myCut(end, 2)+1, min(height(XY), myCut(end, 2)+1+MStep)];

end

mzAxis = Lim(1, 1);
for ii = 1:height(myCut)
    cXY = [XY(myCut(ii, 1):myCut(ii, 2), 1), XY(myCut(ii, 1):myCut(ii, 2), 2)];

    [p{ii}, S{ii}, mu{ii}] = polyfit(cXY(:, 1), cXY(:, 2), n);
    cXY(:, 3) = polyval(p{ii}, cXY(:,1), S{ii}, mu{ii});
    cXY(:, 4) = sqrt((cXY(:, 2) - cXY(:, 3)).^2);
    cXY(cXY(:, 4) > mean(cXY(:, 4)) + 5*std(cXY(:, 4)), :) = [];

    [p{ii}, S{ii}, mu{ii}] = polyfit(cXY(:,1), cXY(:,2), n);
    cXY(:, 3) = polyval(p{ii}, cXY(:,1), S{ii}, mu{ii});
    cXY(:, 4) = sqrt((cXY(:, 2) - cXY(:, 3)).^2);

    R = corrcoef(polyval(p{ii}, cXY(:,1), S{ii}, mu{ii}), cXY(:,2));
    r2{ii} = R(1,2);

    int_old = 0;
    while mzAxis(end, 1) < XY(myCut(ii, 2), 1)
        int = polyval(p{ii}, mzAxis(end, 1), S{ii}, mu{ii})/ratio;
        mzAxis(end+1, 1) = mzAxis(end, 1) + int; %#ok<AGROW>
        if int_old > int
            error("")
        else
            int_old = int;
        end
    end

    nXY = [nXY; cXY];
end

while mzAxis(end, 1) < Lim(2)
    int = polyval(p{ii}, mzAxis(end, 1), S{ii}, mu{ii})/ratio;
    mzAxis(end+1, 1) = mzAxis(end, 1) + int; %#ok<AGROW>
    if int_old > int
        error("")
    else
        int_old = int;
    end
end

myCut = [XY(myCut(:, 1) ,1), XY(myCut(:, 2) ,1)];
stats.Cuts = myCut;
stats.r2 = r2;
stats.p  = p;
stats.S  = S;
stats.mu = mu;

end


