function splitedProfiles = splitProfiles(XY, minsize)
XY(XY(:,2) == 0, 2) = nan;
IdNan = isnan(XY(:,2));
splitedProfiles = {};
nbSplit = 1;
IdS = find(~IdNan, 1, 'first');
if isempty(IdS)
    splitedProfiles{1} = XY;
    return
end


while 1
    IdE = find(IdNan(IdS:end), 1, 'first');
    if isempty(IdE)
        IdE = size(XY, 1);
        if IdE - IdS > minsize
            splitedProfiles{nbSplit} = XY(IdS:IdE, :);
        end
        break

    else
        IdE = IdE + IdS - 2;

    end

    if IdE - IdS > minsize
        splitedProfiles{nbSplit} = XY(IdS:IdE, :);
        nbSplit = nbSplit+1;
    end

    IdS = find(~IdNan(IdE+1:end), 1, 'first');
    if isempty(IdS)
        break
    else
        IdS = IdE + IdS;
    end


end