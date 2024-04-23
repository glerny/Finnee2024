%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [myROIAnalysis, ExportData] = analysis_2(obj, Options)


if nargin == 1
    Options.sigMZ        = 1;
    Options.minTime      = 5;
    Options.WindowTime   = 4;
    Options.WindowMZ     = 2;
    Options.sgfd         = 'average';
    Options.critRes      = 1;
    Options.DynRange     = 0.05;
    Options.Orders       = 1; 

end
myTarget = obj.AditionalInformation.target;
ExportData = {};

variable2_names_types = [["chrom_apex", "double", "%.4f", "m/z"]; ...
    ["chrom_intensity_apex", "double", "%.0f", ""]; ...
    ["chrom_area", "double", "%.2f", "Arb. units"]; ...
    ["chrom_centroid", "double", "%.4f", ""]; ...
    ["chrom_variance", "double", "%.2e", ""];...
    ["chrom_assym", "double", "%.2e", ""];...
    ["Peak_start", "double", "%.2e", ""];...
    ["Peak_end", "double", "%.2e", ""];...
    ["Accu_mass", "double", "%.4f", ""]; ...
    ["ms_variance", "double", "%.4f", ""]; ...
    ["ms_assym", "double", "%.4f", ""]; ...
    ["ms_Peak_start", "double", "%.2e", ""];...
    ["ms_Peak_end", "double", "%.2e", ""];...
    ["Signal", "double", "%.4f", ""]; ...
    ["Noise", "double", "%.4f", ""]; ...
    ["n2n", "double", "%.4f", ""]];
ChromPeak = table('Size', [0, size(variable2_names_types,1)],...
    'VariableNames', variable2_names_types(:,1),...
    'VariableTypes', variable2_names_types(:,2));
ChromPeak.Properties.VariableUnits = variable2_names_types(:,4);
ChromPeak.Properties.VariableDescriptions = variable2_names_types(:,3);
IsEmpty = table('Size', [0, size(variable2_names_types,1)],...
    'VariableNames', variable2_names_types(:,1),...
    'VariableTypes', variable2_names_types(:,2));
IsEmpty.Properties.VariableUnits = variable2_names_types(:,4);
IsEmpty.Properties.VariableDescriptions = variable2_names_types(:,3);
IsEmpty(1, :) = {NaN, NaN, NaN, NaN,  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
NoBasCorr = true;

if isempty(obj.Data)
    myROIAnalysis = IsEmpty;
    ExportData = {};
    return;

end

%% STEP 1. Smoothed XY to Data
XY =  obj.Data; XY(~isfinite(XY)) = 0;
ExportData.X = obj.AxisX.Data;
ExportData.Y = obj.AxisY.Data;
ExportData.OriXY = XY;

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
ExportData.SmoXY = Data;

XYZ = zeros(size(Data));
Y_new = sort(unique([obj.AxisY.Data; ...
    [myTarget.Targetmz - Options.sigMZ*myTarget.PeakStdDev_MS; myTarget.Targetmz + Options.sigMZ*myTarget.PeakStdDev_MS]]));
[X,Y] = meshgrid(obj.AxisX.Data, obj.AxisY.Data);
[Xq,Yq] = meshgrid(obj.AxisX.Data, Y_new);
XYZ = interp2(X, Y, Data, Xq, Yq);
XYZ(XYZ <0) = 0;
Itgt = ...
    [find(Y_new == myTarget.Targetmz - Options.sigMZ*myTarget.PeakStdDev_MS), ...
    find(Y_new == myTarget.Targetmz + Options.sigMZ*myTarget.PeakStdDev_MS)];
XYChr = [obj.AxisX.Data, (trapz(Y_new(Itgt(1):Itgt(2)), XYZ(Itgt(1):Itgt(2), :)))'];
XYChr(isnan(XYChr)) = 0;

%% Baseline correction
if all(XYChr(:, 2) ~= 0)
    [z, w] = doPF(XYChr, Options.Orders);
    Noise = mean(XYChr(w, 2)) + 1*std(XYChr(w, 2));
    XYChr(:, 2) = XYChr(:, 2) - z;
    XYChr(XYChr(:, 2) < 0, 2) = 0;
    NoBasCorr = false;

end
lmChr = LocalMaxima(XYChr, Options.WindowTime, Options.DynRange*max(XYChr(:, 2)));
if isempty(lmChr)
    myROIAnalysis = IsEmpty;
    ExportData = {};
    return;

end

lmChr = unique([XYChr(1, :); lmChr; XYChr(end, :)], "rows");
for ii = 1:size(lmChr,1 )
    lmChr(ii, 3) = findCloser(lmChr(ii, 1), XYChr(:, 1));

end

PEAKS = [];
for ii = 2:size(lmChr,1 )-1
    IdL = find(XYChr(lmChr(ii-1, 3):lmChr(ii, 3), 2) == ...
        min(XYChr(lmChr(ii-1, 3):lmChr(ii, 3), 2)), 1, "last") ...
        + lmChr(ii-1, 3) - 1;
    IdL(2) = find(XYChr(lmChr(ii, 3):lmChr(ii+1, 3), 2) == ...
        min(XYChr(lmChr(ii, 3):lmChr(ii+1, 3), 2)), 1, "first") ...
        + lmChr(ii, 3) - 1;
    if IdL(2) - IdL(1) > Options.minTime
        PEAKS(end+1, :) = [IdL, (XYChr(IdL, 1))', lmChr(ii, 1:2), ChrMoment(XYChr(IdL(1):IdL(2), :), 3)];
    end

end

if isempty(PEAKS)
    myROIAnalysis = IsEmpty;
    ExportData = {};
    return;

end

while 1
    if size(PEAKS, 1) > 1
        
        PEAKS = sortrows(PEAKS, 1);
        % 1. check for the critical resolution
        for jj = 1:size(PEAKS, 1)
            RS = abs(PEAKS(jj, 8) - PEAKS(:, 8))./(2*(sqrt(PEAKS(jj, 9)) + sqrt(PEAKS(:, 9))));
            IdX = find(RS == min(RS(RS~=0)));
            PEAKS(jj, 10) = IdX;
            PEAKS(jj, 11) = RS(IdX);

        end

    else
        break

    end

    if min(PEAKS(PEAKS(:, 11) ~= 0, 11)) <=  Options.critRes
        Id = find(PEAKS(:, 11) == min(PEAKS(PEAKS(:, 11) ~= 0, 11)), 1, 'first');
        Id(2) = PEAKS(Id, 10);
        IdX = [min(PEAKS(Id, 1)) max(PEAKS(Id, 2))];
        mLM = PEAKS(Id, 5:6); mLM = mLM(mLM(:, 2) == max(mLM(:, 2)), :);
        PEAKS(Id, :) = [];
        PEAKS(end+1, 1:9) = [IdX, (XYChr(IdX, 1))', mLM(1), mLM(2), ChrMoment(XYChr(IdX(1):IdX(2), :), 3)];

    else
        break

    end
end

% GET THE BACKGROUND NOISE
if NoBasCorr
   Noise = nan;

end

% GET THE MASS SPECTRA
MassSpectra = {};
for ii = 1:size(PEAKS, 1)
    MassSpectra{ii}(:, 1) = obj.AxisY.Data;
    MassSpectra{ii}(:, 2) = mean(XY(:, PEAKS(ii, 1):PEAKS(ii, 2)), 2);
end

% CALCUL AND RECORD
for ii = 1:size(PEAKS, 1)
    Vector = XYChr(PEAKS(ii, 1):PEAKS(ii, 2), :);
    %MCr  = ChrMoment(Vector);
    lmMS = LocalMaxima(MassSpectra{ii}, Options.WindowMZ, Options.DynRange*max(MassSpectra{ii}(:, 2)));

    if isempty(lmMS)
        [M, IM] = max(MassSpectra{ii}(:, 2));
        lmMS =  [MassSpectra{ii}(IM, 1), M];
        IdMS_start = 1; IdMS_end = numel(MassSpectra{ii}(:, 2));

    else
        lmMS = unique([MassSpectra{ii}(1, :); lmMS; MassSpectra{ii}(end, :)], "rows");
        IdX = findCloser(myTarget.Targetmz, lmMS(:, 1));

        if IdX == 1 | IdX == size(lmMS, 1)
            [M, IM] = max(MassSpectra{ii}(:, 2));
            lmMS =  [MassSpectra{ii}(IM, 1), M];
            IdMS_start = 1; IdMS_end = numel(MassSpectra{ii}(:, 2));

        else

            IdmmMS = zeros(size(lmMS, 1), 1);
            for jj = 1:size(lmMS, 1)
                IdmmMS(jj) = findCloser(lmMS(jj, 1), MassSpectra{ii}(:, 1));
            end

            [~, IdMS_start] = min(MassSpectra{ii}(IdmmMS(IdX-1):IdmmMS(IdX), 2));
            IdMS_start = IdMS_start + IdmmMS(IdX-1) - 1;
            [~, IdMS_end] = min(MassSpectra{ii}(IdmmMS(IdX):IdmmMS(IdX+1), 2));
            IdMS_end = IdMS_end + IdmmMS(IdX) - 1;
        end

    end

    MMS  = ChrMoment(MassSpectra{ii}(IdMS_start:IdMS_end, :));
    lmMS = lmMS(findCloser(MMS(2), lmMS(:, 1)), :);
    ratio = sort(Vector(:, 2));
    ratio(ratio == 0) = [];
    ratio(:, 2)  = nan(size(ratio));
    for kk = 1:numel(Vector(:, 2))
        ratio(kk, 2) = sum(ratio(1:kk, 1))/sum(ratio(kk:end, 1));
        if ratio(kk, 2) > 1; break; end

    end
    ratio(isnan(ratio(:, 2)), :) = [];

    try
        Heq = interp1(ratio(:, 2), ratio(:,1), 1);
        Sig = sum(Vector(:, 2) > Heq);
        n2n = sum(Vector(:, 2) > Heq)/sum(Vector(:, 2) < Heq);

    catch
        Sig = NaN;
        n2n = NaN;
    end


    ChromPeak(ii, :) = {PEAKS(ii, 5), PEAKS(ii, 6), PEAKS(ii, 7), PEAKS(ii, 8), ...
        PEAKS(ii, 9),  (PEAKS(ii, 8) - PEAKS(ii, 5))/ PEAKS(ii, 8), ...
        PEAKS(ii, 3),  PEAKS(ii, 4) , ...
        MMS(2), MMS(3), (MMS(2) - lmMS(1))/MMS(2),...
        MassSpectra{ii}(IdMS_start, 1), MassSpectra{ii}(IdMS_end, 1), ...
        Sig, Noise, n2n};
end

ExportData.TargetChr = XYChr;
ExportData.TargetMS = MassSpectra;
myROIAnalysis = ChromPeak;
end

