function VectorOut = OutliersRemoval_robust(VectorIn)

option.firstAlfa  = 3;
option.secondAlfa = 2;
option.test4outlier = 2;
VectorOut = VectorIn;

if mean(VectorOut) <= option.test4outlier*std(VectorOut)
    while 1
        Ido = isoutlier(VectorOut, "ThresholdFactor", option.firstAlfa);

        if sum(Ido) == 0
            break

        else
            VectorOut(Ido) = [];

        end
    end

    if mean(VectorOut) <= option.test4outlier*std(VectorIn)
        Ido = isoutlier(VectorOut, "ThresholdFactor", option.secondAlfa);
        VectorOut(Ido) = [];

    end
end