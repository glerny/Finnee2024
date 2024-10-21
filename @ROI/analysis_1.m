%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: (19/10/2023)
%       1. Description
%       2. Control when 2 peaks merged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChromPeak = analysis_1(obj, Options)

%% INITIALISATION
tStart = tic; 

if nargin == 1
    Options.minMZ           = 3;
    Options.minTime         = 5;
    Options.WindowMZ        = 2;
    Options.WindowTime      = 5;
    Options.Min_nnz         = 10;
    Options.sgfd            = 'average';
    Options.StepSize        = 1.1; % in step increment tm or mz
    Options.CritRes         = 0.8;

end

variable2_names_types = [["ID2ROI", "double", "%.0f", ""];  ...
    ["PeakNbr", "double", "%.0f", ""];  ...
    ["3D_nnz", "double", "%.0f", ""]; ...
    ["3D_CV", "double", "%.2f", ""]; ...
    ["3D_Kurto", "double", "%.2f", ""]; ...
    ["3D_Qual1", "double", "%.2f", ""]; ...
    ["3D_Qual2", "double", "%.2f", ""]; ...
    ["3D_LMD", "double", "%.2f", ""]; ...
    ["chrom_Time_IA", "double", "%.2f", "min"]; ...
    ["chrom_intensity_apex", "double", "%.0f", ""]; ...
    ["chrom_area", "double", "%.2f", "Arb. units"]; ...
    ["chrom_centroid", "double", "%.4f", ""]; ...
    ["chrom_variance", "double", "%.2e", ""];...
    ["chrom_ptsPerSigm", "double", "%.2f", ""];...
    ["ms_Accu_mass_1", "double", "%.4f", ""]; ...
    ["ms_std_Accu_mass_1", "double", "%.4f", ""]; ...
    ["ms_Accu_mass_2", "double", "%.4f", ""];...
    ["ms_std_Accu_mass_2", "double", "%.4f", ""];...
    ["ms_variance", "double", "%.4f", ""];...
    ["ms_ptsPerSigm", "double", "%.2f", ""]; ...
    ["ROILims", "cell", "%.2f", ""]];

ChromPeak = table('Size', [0, size(variable2_names_types,1)],...
    'VariableNames', variable2_names_types(:,1),...
    'VariableTypes', variable2_names_types(:,2));
ChromPeak.Properties.VariableUnits = variable2_names_types(:,4);
ChromPeak.Properties.VariableDescriptions = variable2_names_types(:,3);
ChromPeak = addprop(ChromPeak,{'Options4Creations','DateOfCreation','cmpTime'}, ...
              {'table','table','table'});

ChromPeak.Properties.CustomProperties.Options4Creations = Options;
ChromPeak.Properties.CustomProperties.DateOfCreation = datetime;

if any(size(obj.Data) < [2*Options.WindowMZ+1 2*Options.WindowTime+1])
    ChromPeak.Properties.CustomProperties.cmpTime = toc(tStart);
    return
end

%% STEP 1. Smoothed XY to Data
XY =  obj.Data; XY(~isfinite(XY)) = 0;
switch Options.sgfd
    case 'low'
        h = sgsdf_2d(-1:1, -1:1, 1, 1, 0);
        Data = filter2(h, XY, 'same');

    case 'average'
        h = sgsdf_2d(-2:2, -2:2, 2, 2, 0);
        Data = filter2(h, XY, 'same');


    case 'high'
        h = sgsdf_2d(-5:5, -5:5, 2, 2, 0);
        Data = filter2(h, XY, 'same');

    case 'none'
        Data = XY;

    otherwise 
        error('Options for sgfd are ''none'', ''low'', ''average'' or ''high''')

end

%% STEP 2. Find local maxima in Data
LocMax = [];
for ii = 2:size(Data, 1)-1
    for jj = 2:size(Data, 2)-1
        Int = [max(1, ii-Options.WindowMZ) min(ii+Options.WindowMZ, size(Data, 1)) ...
            max(1, jj-Options.WindowTime) min(jj+Options.WindowTime, size(Data, 2))];
        if Data(ii, jj) == max(Data(Int(1):Int(2), Int(3):Int(4)), [], 'all', 'omitnan') && ...
                Data(ii, jj) > 0 && ...
                all(XY(Int(1):Int(2), Int(3):Int(4)) > 0, 'all')
            LocMax(end+1, :) = [XY(ii, jj), ii, jj];

        end
    end
end

if isempty(LocMax)
    ChromPeak.Properties.CustomProperties.cmpTime = toc(tStart);
    return

end
LocMax = sortrows(LocMax, -1);

%% STEP 3. get and initialise hyperlimits
Hsurf = zeros([size(XY) + 2, size(LocMax, 1)]);
Hsurf(1, :, :) = -1;
Hsurf(end, :, :) = -1;
Hsurf(:, 1, :) = -1;
Hsurf(:, end, :) = -1;
Hlim = [];

for ii = 1:size(LocMax, 1)
    Hlim(ii, :) = [0 LocMax(ii, 2)+1 LocMax(ii, 3)+1 1 1 1 0 0 0];
    Mask = mkMask(Hlim(ii, 2:6), Hsurf(:, :, ii));
    if any(Hsurf(:, :, :).*Mask ~= 0, 'all')
        error("test")
    end
    Hlim(ii, 1) = sum(XY.*Mask(2:end-1, 2:end-1), 'all');

    list3D = find(1:size(Hsurf, 3) ~= ii);
    for jj = 1:numel(list3D)
        Hsurf(:, :, list3D(jj)) = Hsurf(:, :, list3D(jj)) + ii*Mask;
    end

end

%% STEP 4. Optimise the "best" surface under peaks
while 1
    AddVol = zeros(3, size(LocMax,1));

    for ii = 1:size(LocMax, 1)
        xini = Hlim(ii, 2:6);

        if Hlim(ii, 7) == 0
            x = xini; x(3) = x(3) + Options.StepSize;
            Mask = mkMask(x, Hsurf(:, :, ii));
            if any(Hsurf(:, :, ii).*Mask ~= 0, 'all')
                lVal =  Hsurf(Hsurf(:, :, ii).*Mask ~= 0);
                Hlim(ii, 7) = lVal(1);
                AddVol(ii, 1) = 0;
                
            else
                AddVol(ii, 1) = sum(XY.*Mask(2:end-1, 2:end-1), 'all') - Hlim(ii, 1);

            end
        end
        
        if Hlim(ii, 8) == 0
            x = xini; x(4) = x(4) + Options.StepSize;
            Mask = mkMask(x, Hsurf(:, :, ii));
            if any(Hsurf(:, :, ii).*Mask ~= 0, 'all')
                lVal = Hsurf(Hsurf(:, :, ii).*Mask ~= 0);
                Hlim(ii, 8) = lVal(1);
                AddVol(ii, 2) = 0;
                
            else
                AddVol(ii, 2) = sum(XY.*Mask(2:end-1, 2:end-1), 'all') - Hlim(ii, 1);

            end
        end
        
        if Hlim(ii, 9) == 0
            x = xini; x(5) = x(5) + Options.StepSize;
            Mask = mkMask(x, Hsurf(:, :, ii));
            if any(Hsurf(:, :, ii).*Mask ~= 0, 'all')
                lVal = Hsurf(Hsurf(:, :, ii).*Mask ~= 0);
                Hlim(ii, 9) = lVal(1);
                AddVol(ii, 3) = 0;

            else
                AddVol(ii, 3) = sum(XY.*Mask(2:end-1, 2:end-1), 'all') - Hlim(ii, 1);

            end
        end

    end

    if max(AddVol, [], 'all') <= obj.BlankNoise
        break
    end
    
    [lin, col] = find(AddVol == max(AddVol, [], 'all')); lin = lin(1); col = col(1);
    Hlim(lin, col+3) =  Hlim(lin, col+3) + 1;
    Mask = mkMask(Hlim(lin, 2:6), Hsurf(:, :, lin)); 
    Hlim(lin, 1) = sum(XY.*Mask(2:end-1, 2:end-1), 'all');
    
    Hsurf(Hsurf == lin) = 0;
    list3D = find(1:size(Hsurf, 3) ~= lin);
    for jj = 1:numel(list3D)
        Hsurf(:, :, list3D(jj)) = Hsurf(:, :, list3D(jj)) + lin*Mask;
    end
end

%% STEP 5. Calculate the FOMs
ic = 1;
for ii = 1:size(Hlim, 1)
    Mask = mkMask(Hlim(ii, 2:6), Hsurf(:, :, ii)); 
    myCut = Mask(2:end-1, 2:end-1);
    xy = XY.*myCut;

    if nnz(xy) >= Options.Min_nnz

        %%% get descriptor of the 3D peak
        vector = xy(:);
        vector(vector == 0) = [];
        ratio = sort(vector);
        ratio(:, 2)  = nan(size(vector));
        for kk = 1:numel(vector)
            ratio(kk, 2) = sum(ratio(1:kk, 1))/sum(ratio(kk:end, 1));
            if ratio(kk, 2) > 1; break; end

        end
        ratio(isnan(ratio(:, 2)), :) = [];

        if size(ratio, 1) == 1
            Heq =  max(vector);

        else
            Heq = interp1(ratio(:, 2), ratio(:,1), 1);

        end

        locMax = LocalMaxima3D(xy, [1 1], 'SmoothLow');
        if isempty(locMax)
            LMD = 0;
        else

            LMD = size(locMax, 1)/numel(vector);
        end

        Descriptors = [numel(vector), std(vector)/mean(vector), ...
            mean(vector)/max(vector), Heq/max(vector), ...
            sum(vector>Heq)/sum(vector<Heq), LMD];
         %%% END descriptor of the 3D peak

        y  = obj.AxisY.Data;
        x  = obj.AxisX.Data;
        lmS = LocalMaxima( [y, sum(xy, 2)],  Options.WindowMZ, 0);

        % Check if kkep peak
        if nnz(xy) >= Options.Min_nnz & ...
                ~isempty(lmS) & ...
                sum(sum(xy, 2) ~= 0, 'all') >= Options.minMZ & ...
                sum(sum(xy, 1) ~= 0, 'all') >= Options.minTime

            %%% MEASURE THE MS PARAMETERS (PROJECTION OF 3D PEAK TO M/Z
            %%% AXIS - PART ms1

            
            lmS = lmS(lmS(:, 2) == max(lmS(:, 2)), :);

            CM_MS = ChrMoment([y, sum(xy, 2)], 3);
            MSPeak = [lmS(1, 1), lmS(1, 2), CM_MS(1), CM_MS(2), ...
                CM_MS(3), sqrt(CM_MS(3))/((y(end)-y(1))/numel(y))];
            %%% END OF PART ms1

            %%% MAKE THE XIP FOR THE TARGET MS PEAK
            XIP = zeros(numel(x), 8);
            XIP(:,1) = x;
            for jj = 1:numel(x)
                c_xy =  [y, xy(:, jj)];

                if nnz(c_xy(:, 2)) >= Options.minMZ
                    lmxy = LocalMaxima(c_xy(:, 1:2),  Options.WindowMZ, 0);
                    if ~isempty(lmxy)
                        lmxy = sortrows(lmxy, -2);
                        if ~isempty(lmxy)
                            CM = ChrMoment(c_xy(:, 1:2), 3);
                            XIP(jj, 2:7) = [lmxy(1, 2), CM(1), lmxy(1, 1), CM(2), CM(3), size(lmxy, 1)];

                        end
                    end
                end
            end

            % keep one trailing and ending zeros
            IdS = find(XIP(:, 2) > 0, 1, 'first');
            if isempty(IdS)
                IdS = size(XIP, 1)-1;

            end
            IdS = max(1, IdS-1);
            IdE = find(XIP(:, 2) > 0, 1, 'last');
            if isempty(IdE)
                IdE = size(XIP, 1);

            end
            IdE = min(size(XIP, 1), IdE+1);
            XIP = XIP(IdS:IdE, :);
            % XIP IS MADE!

            lmP = LocalMaxima( [x, sum(xy, 1)'],  Options.WindowTime, 0);
            if ~isempty(lmP) & size(XIP, 1) > Options.minTime

                lmP = lmP(lmP(:, 2) == max(lmP(:, 2)), :);
                Id2M = XIP(:, 2) ~= 0;

                CM_PF = ChrMoment(XIP(:, [1 3]), 3);
                accu = sum(XIP(Id2M,3).*XIP(Id2M,4))/sum(XIP(Id2M,3));
                accu(2) = std(XIP(Id2M,4)-accu(1));
                accu(3) = sum(XIP(Id2M,3).*XIP(Id2M,5))/sum(XIP(Id2M,3));
                accu(4) = std(XIP(Id2M,5)-accu(3));
                FLim = Hlim(ii, 2:6);

                ChromPeak(ic, :) = ...
                    {nan, ic, Descriptors(1), Descriptors(2), Descriptors(3),...
                     Descriptors(4),  Descriptors(5),  Descriptors(6),...
                     lmP(1), lmP(2), CM_PF(1), CM_PF(2), CM_PF(3), ...
                     sqrt(CM_PF(3))/((x(end)-x(1))/numel(x)), ...
                     accu(1), accu(2), accu(3), accu(4), MSPeak(5), MSPeak(6), ...
                     {FLim}};
                if any(table2array(ChromPeak(ic, 2:20)) <0)
                    disp("pp")
                end
                ic = ic + 1;

            end
        end
    end
end

MinResolution = nan(size(ChromPeak, 1), 2);
for ii = 1:size(ChromPeak, 1)
    RS = abs((ChromPeak.chrom_centroid - ChromPeak.chrom_centroid(ii)))./...
        (2*(sqrt(ChromPeak.chrom_variance) + sqrt(ChromPeak.chrom_variance(ii))));
    RS(:, 2) = abs((ChromPeak.ms_Accu_mass_2 - ChromPeak.ms_Accu_mass_2(ii)))./...
        (2*(sqrt(ChromPeak.ms_variance) + sqrt(ChromPeak.ms_variance(ii))));
    RS(:, 3) = sqrt(RS(:, 1).^2 + RS(:, 2).^2);
    IdRs = find(RS(:, 3) == min(RS(RS(:, 3) ~= 0, 3)));

    if ~isempty(IdRs)
        MinResolution(ii, :) = [RS(IdRs, 3), IdRs];
    end
end

if numel(MinResolution) > 2 & any(isnan(MinResolution(:, 1)))
    disp("pp")
end

if any(MinResolution(:, 1) <= Options.CritRes) & numel(MinResolution) > 2
    newOptions = Options;
    newOptions.WindowMZ   = newOptions.WindowMZ+1;
    newOptions.WindowTime = newOptions.WindowTime+1;
    ChromPeak = obj.analysis_1(newOptions);

else
    ChromPeak.CriticalResolution = MinResolution;

end

    function Mask = mkMask(x, XYF)
        Mask = zeros(size(XYF));
        for iM = 1:size(XYF, 1)
            for jM = 1:size(XYF, 2)

                if jM <= x(2)
                    Mask(iM, jM) = isPointInEllipse(iM, jM, x(1), x(2), x(3), x(4));

                else
                    Mask(iM, jM) = isPointInEllipse(iM, jM, x(1), x(2), x(3), x(5));

                end
            end
        end
    end
end

