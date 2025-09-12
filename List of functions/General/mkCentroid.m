function [MS_ctr, MS_prf] = mkCentroid(MS_prf)
wdw = 3;
spkSz = 3;
iStart = 1;
count = 1;
doStop = false;
MS_ctr = [];
MS_prf = spikesRemoval(MS_prf, spkSz);
newMS = {};

while 1
    IdS = find(MS_prf(iStart:end, 2) > 0, 1, 'first');
    if isempty(IdS)
        break

    end
    IdS = max(1, IdS+iStart-2);
    IdE = find(MS_prf(IdS+1:end, 2) == 0, 1, 'first');
    if isempty(IdE)
        doStop = true;
        IdE = height(MS_prf);
        
    else
        IdE = IdE+IdS;
        iStart = IdE +1;
    end

    smallMS = MS_prf(IdS:IdE, :);
    LM = LocalMaxima(smallMS, wdw, 0); 
    if height(LM) == 1
        M = ChrMoment(smallMS, 3);
        M(4) = count;
        MS_ctr = [MS_ctr; M];
        newMS{count} = smallMS;
        count = count + 1;

    elseif height(LM) == 2
        IMax = findCloser(smallMS(:, 1), LM(1, 1));
        IStart = 1;
        IStop = findCloser(smallMS(:, 1), LM(2, 1));
        IdStart = find(smallMS(IStart:IMax, 2) == min(smallMS(IStart:IMax, 2)), 1, 'last') + IStart -1;
        IdEnd = find(smallMS(IMax:IStop, 2) == min(smallMS(IMax:IStop, 2)), 1, 'last') + IMax -1;
        M = ChrMoment(smallMS(IdStart:IdEnd, :), 3);
        M(4) = count;
        MS_ctr = [MS_ctr; M];

        IMax = findCloser(smallMS(:, 1), LM(2, 1));
        IStart = findCloser(smallMS(:, 1), LM(1, 1));
        IStop = height(smallMS);
        IdStart = find(smallMS(IStart:IMax, 2) == min(smallMS(IStart:IMax, 2)), 1, 'last') + IStart -1;
        IdEnd = find(smallMS(IMax:IStop, 2) == min(smallMS(IMax:IStop, 2)), 1, 'last') + IMax -1;
        M = ChrMoment(smallMS(IdStart:IdEnd, :), 3);
        M(4) = count;
        MS_ctr = [MS_ctr; M];

        newMS{count} = smallMS;
        count = count + 1;
    
    elseif height(LM) > 2

        IMax = findCloser(smallMS(:, 1), LM(1, 1));
        IStart = 1;
        IStop = findCloser(smallMS(:, 1), LM(2, 1));
        IdStart = find(smallMS(IStart:IMax, 2) == min(smallMS(IStart:IMax, 2)), 1, 'last') + IStart -1;
        IdEnd = find(smallMS(IMax:IStop, 2) == min(smallMS(IMax:IStop, 2)), 1, 'last') + IMax -1;
        M = ChrMoment(smallMS(IdStart:IdEnd, :), 3);
        M(4) = count;
        MS_ctr = [MS_ctr; M];

        for ii = 2:height(LM)-1
            IMax = findCloser(smallMS(:, 1), LM(ii, 1));
            IStart = findCloser(smallMS(:, 1), LM(ii-1, 1));
            IStop =  findCloser(smallMS(:, 1), LM(ii+1, 1));
            IdStart = find(smallMS(IStart:IMax, 2) == min(smallMS(IStart:IMax, 2)), 1, 'last') + IStart -1;
            IdEnd = find(smallMS(IMax:IStop, 2) == min(smallMS(IMax:IStop, 2)), 1, 'last') + IMax -1;
            M = ChrMoment(smallMS(IdStart:IdEnd, :), 3);
            M(4) = count;
            MS_ctr = [MS_ctr; M];

        end

        IMax = findCloser(smallMS(:, 1), LM(2, 1));
        IStart = findCloser(smallMS(:, 1), LM(1, 1));
        IStop = height(smallMS);
        IdStart = find(smallMS(IStart:IMax, 2) == min(smallMS(IStart:IMax, 2)), 1, 'last') + IStart -1;
        IdEnd = find(smallMS(IMax:IStop, 2) == min(smallMS(IMax:IStop, 2)), 1, 'last') + IMax -1;
        M = ChrMoment(smallMS(IdStart:IdEnd, :), 3);
        M(4) = count;
        MS_ctr = [MS_ctr; M];

        newMS{count} = smallMS;
        count = count + 1;

    end

    if doStop, break, end

end

MS_prf = [];
for ii = 1:numel(newMS)
    MS_prf = [MS_prf; newMS{ii}];
end

end