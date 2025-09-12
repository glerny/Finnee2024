%% DESCRIPTION
% INTERPOLATE2D 
%
%% INPUT PARAMETERS
% *Compulsory*
% *obj      : The Finnee object
% *dts      : The indice to the dataset that contains the original scans
%             (should be profile scans)
% *newAxis  : The master MZ axis. This parameter can be empty,
%
% *Optional*
% *tLim       : Followed by a 2x1 array of numbers (default [0 inf]). Only
%             records scans between tLim(1) and tLim(2)
% *spikes     : Followed by an integer between 0 and 3 (default 2) (see the
%             method @Finnee\FilterDataset) for additional information.
%             Remove spikes in every MS scans If used, where spikes are any
%             peaks in each MS of length equal or lower that the integer.
%             'spikes' followed by 0 allows to to turn off spikes removal.
% TO BE DONE
%
%% OUTPUT PARAMETERS
% *obj      : The Finnee object. !Important, if not input parameter is used
%             the saved Finnnee object will still be modified.
%
%% EXAMPLES
% myFinnee = myFinnee(
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj, Axis] = Interpolate2D(obj, dts, Axis, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options

narginchk(2, inf)
tStart = tic;

if nargin == 2 & isstruct(dts)
    options = dts;
    dts = options.parentDataset;
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    Axis = options.Axis;
    options.ParallelMe = true;

else
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);

    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    options = checkVarargin(dts, varargin{:});
    if nargin == 2, Axis = {}; end
    options.Axis = Axis;
end

Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
if isempty(Id2Dataset), error(""), end
if obj.Datasets.HasMasterMZ(Id2Dataset), error(""), end
if ~obj.Datasets.isProfileScan{Id2Dataset}, error(""), end
if ~obj.Datasets.IsMS1{Id2Dataset}, error(""), end


% Create new dat file
ii = 1;
while 1
    newDataset = ['Dataset', num2str(ii)];
    if ~any(strcmp(newDataset, obj.Datasets.Name)), break; end
    ii = ii + 1;
end
nds = ii;
obj.Datasets = [obj.Datasets; ...
    {newDataset, '', true, true, true, false, false, ...
    ['Dataset', num2str(dts)], obj.Datasets.Labels{Id2Dataset}, ...
    options, 'Master_mz_axis', {}}];

%% Load or create first mz Axis
% Get the index of the scan with the highest total intensity
fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
fclose(fidRead);
Profiles = array2table(Profiles);
ColumnNames = {}; ColumnUnits = {};
for ii = 1:numel(infoDataset.Profiles.column)
    ColumnNames{ii}     = infoDataset.Profiles.column{ii}.name;
    ColumnLabels{ii}    = infoDataset.Profiles.column{ii}.labelUnits;
    ColumnUnits{ii}     = infoDataset.Profiles.column{ii}.units;
    ColumnFO{ii}        = infoDataset.Profiles.column{1, 1}.formatUnits;
end
Profiles.Properties.VariableNames = ColumnNames;
Profiles.Properties.VariableUnits = ColumnUnits;
Id1 = find(Profiles.Non_Null_Elements_Profile == max(Profiles.Non_Null_Elements_Profile));
Id2 = Profiles.Total_Ions_Profile(Id1) == max(Profiles.Total_Ions_Profile(Id1));
IdAxis = Id1(Id2);

ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdAxis, 1)), '.dat']);
fidRead = fopen(ScanName, 'rb');
myMS = fread(fidRead, ...
    [Profiles.Length_MS_Scan_Profile(IdAxis) Profiles.Column_MS_Scan_Profile(IdAxis)],...
    "double");
fclose(fidRead);

[mzAxis, data4axis, stats] ...
                                            = mkMasterMZ(myMS, options.split, options.orderPolynomial, options.YLim, 1);
AdditionalInformation.FirstMzAxis.Data4Axis = data4axis;
AdditionalInformation.FirstMzAxis.stats     = stats;
AdditionalInformation.FirstMzAxis.n         =  size(data4axis, 1);

myAxis.FirstMzAxis.p      = options.orderPolynomial;
myAxis.FirstMzAxis.stats  = stats;
myAxis.FirstMzAxis.Data   = mzAxis;

doSecondaryMz = false;
myAxis.TimeAxis.Data    = Profiles.Scan_Time;

if isempty(Axis)
    myAxis.TimeAxis.Data    = Profiles.Scan_Time;

else
    if isfield(Axis, 'TimeAxis')
        myAxis.TimeAxis.Data = Axis.TimeAxis;

    end
    if isfield(Axis, 'mzAxis')
        doSecondaryMz = true;
        myAxis.SecondMzAxis.Data = Axis.mzAxis;

    end
end

timeAxis = myAxis.TimeAxis.Data;
mzInterval = [inf 0];
mkdir(fullfile(obj.Path2Fin, newDataset));
mkdir(fullfile(obj.Path2Fin, newDataset, 'Scans'));
newProfiles = [];

if doSecondaryMz
    mzAxis = myAxis.SecondMzAxis.Data;

else
    mzAxis = myAxis.FirstMzAxis.Data;
end
newSpectra = mzAxis; newSpectra(:, 7) = 0;
fAxis = myAxis.FirstMzAxis.Data;

%% Interpolate to MZaxis, average if needed & calcualted base profiles
if ~options.ParallelMe, h = waitbar(0,'processing scans'); end

counter = 1;
ScNu = 0;
CurSeq = zeros(numel(mzAxis), 2);
scanNum = 0;

BNV = []; % Black noise vector
for ii = 1:numel(timeAxis)
    if ~options.ParallelMe, waitbar(ii/numel(timeAxis)); end

    IdS = find (timeAxis(ii) == Profiles.Scan_Time);

    if ~isempty(IdS) % If same axis
        ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
            'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS)), '.dat']);
        fidRead = fopen(ScanName, 'rb');
        myMS = fread(fidRead, ...
            [Profiles.Length_MS_Scan_Profile(IdS) Profiles.Column_MS_Scan_Profile(IdS)],...
            "double");
        fclose(fidRead);

        vq =  interp1(myMS(:,1), myMS(:,2), mzAxis,...
            options.mth4interp1);
        vq(~isfinite(vq)) = 0;
        myMS = [mzAxis, vq];

    else

        % find n-2:n+2 in old axis
        IdS = findCloser(timeAxis(ii), Profiles.Scan_Time);
        if IdS < 2, continue; end
        if IdS > numel(Profiles.Scan_Time)-1; break; end

        if timeAxis(ii) - Profiles.Scan_Time(IdS) < 0
            IdS1 = IdS-1;
            IdS2 = IdS;

        else
            IdS1 = IdS;
            IdS2 = IdS+1;
        end

        ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
            'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS1)), '.dat']);
        fidRead = fopen(ScanName, 'rb');
        myMS = fread(fidRead, ...
            [Profiles.Length_MS_Scan_Profile(IdS1) Profiles.Column_MS_Scan_Profile(IdS1)],...
            "double");
        fclose(fidRead);

        vq =  interp1(myMS(:,1), myMS(:,2), mzAxis,...
            options.mth4interp1);

        ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
            'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS2)), '.dat']);
        fidRead = fopen(ScanName, 'rb');
        myMS = fread(fidRead, ...
            [Profiles.Length_MS_Scan_Profile(IdS2) Profiles.Column_MS_Scan_Profile(IdS2)],...
            "double");
        fclose(fidRead);

        vq(:, 2) =  interp1(myMS(:,1), myMS(:,2), mzAxis,...
            options.mth4interp1);
        vq(~isfinite(vq)) = 0;
     

         vq2D = interp1(Profiles.Scan_Time(IdS1:IdS2), vq', timeAxis(ii));
         vq2D(~isfinite(vq2D)) = 0;
         myMS = [mzAxis, vq2D'];

    end

    % Filter spikes if needed
    if ~isempty(myMS) && options.RemSpks
        spkSz = options.SpkSz;
        myMS    = spikesRemoval(myMS, spkSz );
    end

    mii =  ii;
    mTime = timeAxis(ii);
    ScanNumber = scanNum;
    scanNum = scanNum + 1;

    %Record MZSpectra
    IdNnz = myMS(:, 2) > 0;
    newSpectra(:,2) = newSpectra(:,2) + double(IdNnz);
    IdStart = CurSeq(:, 2) == 0 & IdNnz;
    CurSeq(IdStart, 1) = ScanNumber+1;
    CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

    % Check Length
    IdMatch = newSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
    newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
    newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
    CurSeq(~IdNnz, 2) = 0;
    newSpectra(:,5) = newSpectra(:,5) + myMS(:, 2);
    newSpectra(:,6) = max([newSpectra(:,6), myMS(:, 2)], [], 2);
    
    % reduced trailing zero in excess
    if ~isempty(myMS)
        provMat      = [myMS(2:end, 2); 0];
        provMat(:,2) = myMS(:, 2);
        provMat(:,3) = [0; myMS(1:end-1, 2)];
        myMS         = myMS(sum(provMat, 2) > 0, :);
    end

    if min(myMS(:,1)) < mzInterval(1), mzInterval(1) = min(myMS(:,1)); end
    if max(myMS(:,1)) > mzInterval(2), mzInterval(2) = max(myMS(:,1)); end

    if isempty(myMS)| size(myMS, 1) == 1 | sum(myMS(:, 2)) == 0
        newProfiles(ScanNumber+1, :) = [ScanNumber, ...
            mTime,  nan, 0, sum(myMS(:,2)), 0, ...
            nnz(myMS(:,2)), size(myMS), nan];

    else

        M = ChrMoment(myMS, 2);
        newProfiles(ScanNumber+1, :) = [ScanNumber, ...
            mTime, ...
            nan, max(myMS(:,2)), sum(myMS(:,2)), ...
            min(myMS(myMS(:, 2) ~= 0, 2)), nnz(myMS(:,2)), ...
            size(myMS), M(2)];

    end

    fileName = fullfile(obj.Path2Fin, newDataset, 'Scans', ['Scan#', ...
        num2str(ScanNumber), '.dat']);
    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
    fwrite(fidWriteDat, myMS(:), "double");
    fclose(fidWriteDat);
end
try close(h); catch, end
IdMatch = newSpectra(:,3) < CurSeq(:, 2);
newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
newSpectra(:,7) = newSpectra(:,5)./newSpectra(:,2);

provMat = [newSpectra(2:end, 5); 0];
provMat(:,2) = newSpectra(:, 5);
provMat(:,3) = [0; newSpectra(1:end-1, 2)];
newSpectra = newSpectra(sum(provMat, 2) > 0, :);
short_mzAxis =  newSpectra(:,1);

BNV(isoutlier(BNV)) = [];
minNoise = mean(BNV) + 2*std(BNV);

infoDataset.dateofcreation = datetime("now");
infoDataset.AdditionalInformation = AdditionalInformation;
infoDataset.AdditionalInformation.nonEndingErrors{1} = '';
infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
infoDataset.Label = obj.Datasets.Labels{end};
infoDataset.mzInterval  =  mzInterval;
infoDataset.minNoise    =  minNoise;
infoDataset.IntThres4P8 =  options.Sig2Nois*minNoise;
infoDataset.Options = options;
infoDataset.Profiles.size = size(newProfiles);

infoDataset.Spectra.size = size(newSpectra);
infoDataset.Spectra.column{1}.name = 'm/z axis';
infoDataset.Spectra.column{1}.labelUnits = obj.Datasets.Labels{end}.AxisY.Label;
infoDataset.Spectra.column{1}.units = '';
infoDataset.Spectra.column{1}.formatUnits = obj.Datasets.Labels{end}.AxisY.fo;
infoDataset.Spectra.column{1}.ShowMe = false;
infoDataset.Spectra.column{2}.name = 'Non_Zeros_Elements';
infoDataset.Spectra.column{2}.labelUnits = '';
infoDataset.Spectra.column{2}.units = '';
infoDataset.Spectra.column{2}.formatUnits = '%0.0f';
infoDataset.Spectra.column{2}.ShowMe = false;
infoDataset.Spectra.column{3}.name = 'Max_Continuous_Profile';
infoDataset.Spectra.column{3}.labelUnits = '';
infoDataset.Spectra.column{3}.units = '';
infoDataset.Spectra.column{3}.formatUnits = '%0.0f';
infoDataset.Spectra.column{3}.ShowMe = false;
infoDataset.Spectra.column{4}.name = 'Id_Start_for_Max_Continuous_Profile';
infoDataset.Spectra.column{4}.labelUnits = '';
infoDataset.Spectra.column{4}.units = '';
infoDataset.Spectra.column{4}.formatUnits = '%0.0f';
infoDataset.Spectra.column{4}.ShowMe = false;
infoDataset.Spectra.column{5}.name = 'Total_Intensity_Spectrum';
infoDataset.Spectra.column{5}.labelUnits = obj.Datasets.Labels{end}.AxisZ.Label;
infoDataset.Spectra.column{5}.units =  obj.Datasets.Labels{end}.AxisZ.Label;
infoDataset.Spectra.column{5}.formatUnits =  obj.Datasets.Labels{end}.AxisZ.fo;
infoDataset.Spectra.column{5}.ShowMe = true;
infoDataset.Spectra.column{6}.name = 'Base_Intensity_Spectrum';
infoDataset.Spectra.column{6}.labelUnits = obj.Datasets.Labels{end}.AxisZ.Label;
infoDataset.Spectra.column{6}.units = obj.Datasets.Labels{end}.AxisZ.Units;
infoDataset.Spectra.column{6}.formatUnits = obj.Datasets.Labels{end}.AxisZ.fo;
infoDataset.Spectra.column{6}.ShowMe = true;
infoDataset.Spectra.column{7}.name = 'Average_Intensity';
infoDataset.Spectra.column{7}.labelUnits = obj.Datasets.Labels{end}.AxisZ.Label;
infoDataset.Spectra.column{7}.units = obj.Datasets.Labels{end}.AxisZ.Units;
infoDataset.Spectra.column{7}.formatUnits = obj.Datasets.Labels{end}.AxisZ.fo;
infoDataset.Spectra.column{7}.ShowMe = true;


save(fullfile(obj.Path2Fin, newDataset, 'infoDataset.mat'), 'infoDataset');
%Record Spectra
fileName = fullfile(obj.Path2Fin, newDataset, 'Profiles.dat');
[fidWriteDat, errmsg]  = fopen(fileName, 'wb');
fwrite(fidWriteDat, newProfiles(:), "double");
fclose(fidWriteDat);

fileName = fullfile(obj.Path2Fin, newDataset, 'Spectra.dat');
[fidWriteDat, errmsg]  = fopen(fileName, 'wb');
fwrite(fidWriteDat, newSpectra(:), "double");
fclose(fidWriteDat);

fileName = fullfile(obj.Path2Fin, newDataset, 'MasterMZAxis.dat');
[fidWriteDat, errmsg]  = fopen(fileName, 'wb');
fwrite(fidWriteDat, mzAxis, "double");
fclose(fidWriteDat);

myFinnee = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')

%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECKVARARGIN
    function options = checkVarargin(dts, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.

        % 1- verify the input dataset
        options = {};

        % 2- Defaults optional parameters
        options.RemSpks         = false;
        options.SpkSz           = 0;
        options.mth4interp1     = 'linear';
        options.MstAxis         = false;
        options.YLim            = infoDataset.mzInterval;
        options.orderPolynomial = 3;
        options.split           = 10;
        options.parentDataset   = dts;
        options.AverageScans    = 1;
        options.Spk4noise       = 3;
        options.Sig2Nois        = 3;
        options.MinMinNoise     = [100, 10];
        options.ParallelMe      = false;
        options.ScanRate        = 'average';

        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);

        tgtIx = input('tLim');
        if ~isempty(tgtIx)
            tLim         = varargin{tgtIx +1};
            options.XLim = [min(tLim) max(tLim)];
        end

        tgtIx = input('spikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemSpks = false;
            else
                options.RemSpks = true;
                options.SpkSz  =  spks;
            end
        end

        tgtIx = input('TimeSpikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemTmSpks = false;
            else
                options.RemTmSpks = true;
                options.TmSpkSWd  =  spks;
            end
        end

        tgtIx = input('ScanRate');
        if ~isempty(tgtIx)
            %TODO: test scanRate for 'low', 'normal', 'high'
            options.ScanRate = varargin{tgtIx +1};
        end
    end

    function IdAl = findCloserIndices(Axis_1, Axis_2)
        Axis_1(:, 3) = 0;
        Axis_1(:, 4) = (1:height(Axis_1))';
        Axis_2(:, 3) = 1;
        Axis_2(:, 4) = (1:height(Axis_2))';
        Axis = sortrows([Axis_1; Axis_2], 1);
        Merged = [];

        while 1
            Axis = sortrows(Axis, 1);
            ThisStep = [];
            Id1 = find(Axis(1:end-1, 3) == 0 & Axis(2:end, 3) == 1);
            Axis(Id1, 5:7) = Axis(Id1+1, [1:2, 4]);
            Id2 = find(Axis(2:end, 3) == 0 & Axis(1:end-1, 3) ~= 0);
            Axis(Id2+1, 8:10) = Axis(Id2, [1:2, 4]);

            IdX = Axis(:, 5) ~= 0 & Axis(:, 8) == 0;
            ThisStep = [ThisStep; Axis(IdX, [1 4 5 6 7])];
            Axis(IdX, :) = [];

            IdX = Axis(:, 5) == 0 & Axis(:, 8) ~= 0;
            ThisStep = [ThisStep; Axis(IdX, [1 4 8 9 10])];
            Axis(IdX, :) = [];

            IdX = (Axis(:, 5) ~= 0 & Axis(:, 8) ~= 0) & ...
                (abs(Axis(:, 5) - Axis(:, 1)) <= abs(Axis(:, 8) - Axis(:, 1)));
            ThisStep = [ThisStep; Axis(IdX, [1 4 5 6 7])];
            Axis(IdX, :) = [];

            IdX = (Axis(:, 5) ~= 0 & Axis(:, 8) ~= 0) & ...
                (abs(Axis(:, 5) - Axis(:, 1)) > abs(Axis(:, 8) - Axis(:, 1)));
            ThisStep = [ThisStep; Axis(IdX, [1 4 8 9 10])];
            Axis(IdX, :) = [];
            Axis(:, 5:end) =  [];

            if sum(ThisStep(:, 4)) == 0; break; end

            Merged = [Merged; ThisStep];
        end

        Merged = sortrows(Merged, 1);
        IdAl = Merged(:, [2, 5]);
    end


end