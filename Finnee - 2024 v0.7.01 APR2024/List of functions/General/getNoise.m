function noise = getNoise(XY, options)

if nargin == 1
    options.DegPolynomial   = 2;
    options.Window          = 6;

end

noiseVector = [];
for ii = 1:size(XY, 1)-options.Window+1
    xy = XY(ii:ii+options.Window-1, :);
    p = polyfit(xy(:,1), xy(:,2), options.DegPolynomial);
    noiseVector(ii) = max(xy(:,2)-polyval(p, xy(:,1))) - min(xy(:,2)-polyval(p, xy(:,1)));

end

if mean(noiseVector) > 2*std(noiseVector)
    noise = mean(noiseVector) + 2*std(noiseVector);

else
    while 1
        io = isoutlier(noiseVector);
        if ~any(io)
            break
        end
        noiseVector(io) = [];

    end
    if mean(noiseVector) > 2*std(noiseVector)
        noise = mean(noiseVector) + 2*std(noiseVector);

    else
        io = isoutlier(noiseVector, "ThresholdFactor", 2);
        noiseVector(io) = [];
        noise = mean(noiseVector) + 2*std(noiseVector);
    end
end
