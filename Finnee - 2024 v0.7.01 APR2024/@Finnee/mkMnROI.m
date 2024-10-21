%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2021 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myROIs = mkMnROI(obj, dts, limits, thresholds, FullFolder)
% TODO: help and description

Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
if nargin == 3
    mzmin  = 0;
    chmin  = 0;
    Intmin = 0;
    FullFolder = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'ROIs');

elseif nargin == 4
    FullFolder = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'ROIs');
    mzmin  = thresholds(1);
    chmin  = thresholds(2);
    Intmin = thresholds(3);

else
    mzmin  = thresholds(1);
    chmin  = thresholds(2);
    Intmin = thresholds(3);

end


% definition
infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;

fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
fclose(fidRead);
Profiles = array2table(Profiles);
ColumnNames = {}; ColumnUnits = {};
for ii = 1:numel(infoDataset.Profiles.column)
    ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
    ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
end
Profiles.Properties.VariableNames = ColumnNames;
Profiles.Properties.VariableUnits = ColumnUnits;

fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Spectra.dat'), 'rb');
Spectra = fread(fidRead, [infoDataset.Spectra.size], "double");
fclose(fidRead);
Spectra = array2table(Spectra);
ColumnNames = {}; ColumnUnits = {};
for ii = 1:numel(infoDataset.Spectra.column)
    ColumnNames{ii} = infoDataset.Spectra.column{ii}.name;
    ColumnUnits{ii} = infoDataset.Spectra.column{ii}.units;
end
Spectra.Properties.VariableNames = ColumnNames;
Spectra.Properties.VariableUnits = ColumnUnits;

infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;

mzAxis = Spectra.("m/z axis");

myROIs = {};
limits = sortrows(limits,'ID','ascend');
limits.Done = false(size(limits.ID));

if ~exist(FullFolder, 'dir')
    mkdir(FullFolder);

else
    rmdir(FullFolder, 's')
    fullFileName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'myROIs.mat');

    if isfile(fullFileName)
        delete(fullFileName)
    
    end
    mkdir(FullFolder);
   
end

myROIs.AxisX       = Profiles.Scan_Time; AxisX = Profiles.Scan_Time;
myROIs.AxisY       = Spectra.("m/z axis"); AxisY = Spectra.("m/z axis");
myROIs.minNoise    = infoDataset.minNoise;
myROIs.NoiseVector = infoDataset.minNoise;
myROIs.myPath     = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset});
myROIs.Path2ROIs  = FullFolder;
myROIs.FinneeFile = obj;
myROIs.Dataset    = dts;
limits.sizeROI    = zeros(size(limits, 1), 2);


indMZ = zeros(size(limits, 1), 2);
Data4ROI{size(limits, 1)} = [];
X{size(limits, 1)} = [];
Y{size(limits, 1)} = [];

for ii = 1:size(limits, 1)
    ix = find(AxisY >= limits.mzMin(ii), 1, 'first'); indMZ(ii, 1) = max(1, ix);
    ix = find(AxisY <= limits.mzMax(ii), 1, 'last');
    indMZ(ii, 2) = min(ix, length(AxisY));
    X{ii} = AxisY(indMZ(ii, 1):indMZ(ii, 2));
end


TLim      = [min(limits.TimeMin), max(limits.TimeMax)];
IdTLim1 = find(AxisX < TLim(1), 1, 'last');
if isempty(IdTLim1), IdTLim1 = 1; end

IdTLim2 = find(AxisX >  TLim(2), 1, 'first');
if isempty(IdTLim2), IdTLim2 = length(AxisX); end

IdTLim = [IdTLim1 IdTLim2];
IsOpenRois = false(size(limits, 1), 1);

for ii = IdTLim(1):IdTLim(2)

    ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);
    fMS = Spectra.("m/z axis");

    if ~isempty(myMS)

        [~, Id1] = intersect(fMS, myMS(:,1));
        if numel(Id1) < numel(myMS(:,1)), error(""); end
        fMS(Id1, 2) = myMS(:,2);
    else

        fMS(:, 2) = 0;
    end

    Id = limits.TimeMin <= AxisX(ii) & limits.TimeMax >= AxisX(ii);
    fId = find(Id);
    for jj = 1:length(fId)
        Data4ROI{fId(jj)}(:, end+1) = fMS(indMZ(fId(jj),1):indMZ(fId(jj),2), 2);
        IsOpenRois(fId(jj)) = true;
    end

    cId = find(~Id & IsOpenRois);
    for jj = 1:length(cId)

        if all(size(Data4ROI{cId(jj)}) > [chmin mzmin]) & max(Data4ROI{cId(jj)}, [], 'all') > Intmin

            fileName = fullfile(FullFolder, ['ROI#', num2str(limits.ID(cId(jj))), '.dat']);

            [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
            fwrite(fidWriteDat, Data4ROI{cId(jj)}(:), "double");
            limits.sizeROI(cId(jj), :) = size(Data4ROI{cId(jj)});
            fclose(fidWriteDat);


            % DESCRIPTIVE VARIABLES (non-null data)
            % 1. IdROI | nnz | CV(std/mean) | Kurto (mean/Ampl) |
            % ROIQual1 | ROIQual1 | LocMaxDen

            vector = Data4ROI{cId(jj)}(:);
            vector(vector == 0) = [];
            ratio = sort(vector);
            ratio(:, 2)  = nan(size(vector));
            for kk = 1:numel(vector)
                ratio(kk, 2) = sum(ratio(1:kk, 1))/sum(ratio(kk:end, 1));
                if ratio(kk, 2) > 1; break; end

            end
            ratio(isnan(ratio(:, 2)), :) = [];

            try
                Heq = interp1(ratio(:, 2), ratio(:,1), 1);
            catch
                Heq = NaN;
            end

            try
                locMax = LocalMaxima3D(Data4ROI{cId(jj)}, [1 1], 'SmoothLow');
                if isempty(locMax)
                    LMD = 0;
                    nlm = 0;
                else

                    LMD = size(locMax, 1)/numel(vector);
                    nlm = size(locMax, 1);
                end
            catch
                LMD = NaN;
                nlm = NaN;

            end


            Descriptors(cId(jj), :) = [limits.ID(cId(jj)), numel(vector), ...
                std(vector)/mean(vector), mean(vector)/max(vector), ...
                Heq/max(vector), sum(vector>Heq)/sum(vector<Heq), LMD, nlm];
            % END OF DESCRIPTIVE VARIABLE

            Data4ROI{cId(jj)} = [];
            IsOpenRois(cId(jj)) = false;
            limits.Done(cId(jj)) = true;
        end
    end
end

%% CLOSE THE ROI STILL OPEN
cId = find(IsOpenRois);
for jj = 1:length(cId)

    if all(size(Data4ROI{cId(jj)}) > [chmin mzmin]) & max(Data4ROI{cId(jj)}, [], 'all') > Intmin
        fileName = fullfile(FullFolder, ['ROI#', num2str(limits.ID(cId(jj))), '.dat']);
        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
        fwrite(fidWriteDat, Data4ROI{cId(jj)}(:), "double");
        limits.sizeROI(cId(jj), :) = size(Data4ROI{cId(jj)});
        fclose(fidWriteDat);

        % DESCRIPTIVE VARIABLES (non-null data)
        % 1. IdROI | nnz | CV(std/mean) | Kurto (mean/Ampl) |
        % ROIQual1 | ROIQual1 | LocMaxDen

        vector = Data4ROI{cId(jj)}(:);
        vector(vector == 0) = [];
        ratio = sort(vector);
        ratio(:, 2)  = nan(size(vector));
        for kk = 1:numel(vector)
            ratio(kk, 2) = sum(ratio(1:kk, 1))/sum(ratio(kk:end, 1));
            if ratio(kk, 2) > 1; break; end

        end
        ratio(isnan(ratio(:, 2)), :) = [];
        Heq = interp1(ratio(:, 2), ratio(:,1), 1);

        locMax = LocalMaxima3D(Data4ROI{cId(jj)}, [1 1], 'SmoothLow');
        if isempty(locMax)
            LMD = 0;
        else

            LMD = size(locMax, 1)/numel(vector);
        end


        Descriptors(cId(jj), :) = [limits.ID(cId(jj)), numel(vector), ...
            std(vector)/mean(vector), mean(vector)/max(vector), ...
            Heq/max(vector), sum(vector>Heq)/sum(vector<Heq), LMD, size(locMax, 1)];
        % END OF DESCRIPTIVE VARIABLE

        Data4ROI{cId(jj)} = [];
        IsOpenRois(cId(jj)) = false;
        limits.Done(cId(jj)) = true;
    end
end
DescriptiveVars = array2table(Descriptors, 'VariableNames',...
    {'IdROI', 'nnz','CV','Kurtosis', 'ROIQual1', 'ROIQual2', ...
    'LocalMaxDens', 'LocalMax'});

myROIs.Targets    = limits;
DescriptiveVars(DescriptiveVars.nnz == 0, :) = [];
myROIs.DescriptiveVars = DescriptiveVars;
Id2Rem = all(myROIs.Targets.sizeROI == [0, 0], 2);
myROIs.Targets(Id2Rem, :) = [];
save(fullfile(myROIs.myPath, 'myROIs.mat'), 'myROIs')

end