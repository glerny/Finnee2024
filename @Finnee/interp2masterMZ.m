%% DESCRIPTION
% ALIGN2NEWMZ is used to normalise every MS scans in a datasets to a
% common mz axis.
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
% *meth4int'  : Followed by a string to specify an alternative
%             interpolation method: 'nearest', 'next', 'previous',
%             'linear','spline','pchip', or 'cubic'. The default method is
%             'linear'. help interp1 for more information on the
%             interpolation function. Changing this option is not advised.
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

function [obj, Axis] = interp2masterMZ(obj, dts, mzAxis, varargin)

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
    if ~ options.MstAxis, mzAxis = []; end
    options.ParallelMe = true;

else
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);

    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    options = checkVarargin(dts, varargin{:});

    if nargin == 2, mzAxis = []; end
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

%% Load or create master mz Axis
if isempty(mzAxis)

    % Get the index of the scan with the highest total intensity
    fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
    Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
    fclose(fidRead);
    Profiles = array2table(Profiles);
    ColumnNames = {}; ColumnUnits = {};
    for ii = 1:numel(infoDataset.Profiles.column)
        ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
        ColumnLabels{ii} = infoDataset.Profiles.column{ii}.labelUnits;
        ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
        ColumnFO{ii}    = infoDataset.Profiles.column{1, 1}.formatUnits;
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

    [mzAxis, r2, data4axis] = mkMasterMZ(myMS, options.orderPolynomial, options.YLim);
    AdditionalInformation.MasterMZAxis.Data4Axis = data4axis;
    AdditionalInformation.MasterMZAxis.p = options.orderPolynomial;
    AdditionalInformation.MasterMZAxis.r2 = r2;
    AdditionalInformation.MasterMZAxis.n  =  size(data4axis, 1);
    Axis.Data4Axis = data4axis;
    Axis.p = options.orderPolynomial;
    Axis.r2 = r2;

end
Axis.Axis = mzAxis;

mzInterval = [inf 0];
mkdir(fullfile(obj.Path2Fin, newDataset));
mkdir(fullfile(obj.Path2Fin, newDataset, 'Scans'));
newProfiles = [];
newSpectra = mzAxis; newSpectra(:, 7) = 0;

%% Align to MZaxis, average if needed & calcualted base profiles
if ~options.ParallelMe, h = waitbar(0,'processing scans'); end

counter = 1;
skipRecording = false;
ScNu = 0;
CurSeq = zeros(numel(mzAxis), 2);
for ii = 1:numel(Profiles.Scan_Index)
    if ~options.ParallelMe, waitbar(ii/numel(Profiles.Scan_Index)); end

    ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);

    vq =  interp1(myMS(:,1), myMS(:,2), mzAxis,...
        options.mth4interp1);
    vq(~isfinite(vq)) = 0;
    myMS = [mzAxis, vq];

    % Filter spikes if needed
    if ~isempty(myMS) && options.RemSpks
        spkSz = options.SpkSz;
        myMS    = spikesRemoval(myMS, spkSz );
    end

    if options.AverageScans > 1
        if counter == 1
            mMS = myMS;
            mii = ii;
            counter = counter + 1;
            skipRecording = true;
            mTime = Profiles.Scan_Time(ii);

        elseif counter < options.AverageScans
            mMS(:, 2) = mMS(:, 2) +  myMS(:,2);
            mii = mii + ii;
            counter = counter + 1;
            skipRecording = true;
            mTime = mTime + Profiles.Scan_Time(ii);

        else
            mMS(:, 2) = mMS(:, 2) +  myMS(:,2);
            mii = mii + ii;
            mMS(:,2) = mMS(:,2)/counter;
            mii = round(mii/counter);
            mTime = mTime + Profiles.Scan_Time(ii);
            mTime = mTime/counter;
            counter = 1;
            ScanNumber = ScNu;
            ScNu = ScNu + 1;
            skipRecording = false;

        end
    else
            mMS =  myMS;
            mii =  ii;
            mTime = Profiles.Scan_Time(ii);
            ScanNumber = ScNu;
            ScNu = ScNu + 1;

    end

    if skipRecording, continue; end

    % find and remove spikes
    if options.RemSpks
        spkSz  = options.SpkSz;
        mMS = spikesRemoval(mMS, spkSz );
    end

    %Record MZSpectra
    IdNnz = mMS(:, 2) > 0;
    newSpectra(:,2) = newSpectra(:,2) + double(IdNnz);
    IdStart = CurSeq(:, 2) == 0 & IdNnz;
    CurSeq(IdStart, 1) = ScanNumber+1;
    CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

    % Check Length
    IdMatch = newSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
    newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
    newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
    CurSeq(~IdNnz, 2) = 0;
    newSpectra(:,5) = newSpectra(:,5) + mMS(:, 2);
    newSpectra(:,6) = max([newSpectra(:,6), mMS(:, 2)], [], 2);
    
    % reduced trailing zero in excess
    if ~isempty(mMS)
        provMat      = [mMS(2:end, 2); 0];
        provMat(:,2) = mMS(:, 2);
        provMat(:,3) = [0; mMS(1:end-1, 2)];
        mMS         = mMS(sum(provMat, 2) > 0, :);
    end

    if min(mMS(:,1)) < mzInterval(1), mzInterval(1) = min(mMS(:,1)); end
    if max(mMS(:,1)) > mzInterval(2), mzInterval(2) = max(mMS(:,1)); end

    if isempty(mMS)| size(mMS, 1) == 1 | sum(mMS(:, 2)) == 0
        newProfiles(ScanNumber+1, :) = [ScanNumber, ...
            mTime,  nan, 0, sum(mMS(:,2)), 0, ...
            nnz(mMS(:,2)), size(mMS), nan];

    else

        M = ChrMoment(mMS, 2);
        newProfiles(ScanNumber+1, :) = [ScanNumber, ...
            mTime, ...
            nan, max(mMS(:,2)), sum(mMS(:,2)), ...
            min(mMS(mMS(:, 2) ~= 0, 2)), nnz(mMS(:,2)), ...
            size(mMS), M(2)];

    end

    fileName = fullfile(obj.Path2Fin, newDataset, 'Scans', ['Scan#', ...
        num2str(ScanNumber), '.dat']);
    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
    fwrite(fidWriteDat, mMS(:), "double");
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

% Find mz lines with only spikes
IdSpikes = ~isnan(newSpectra(:, 6)) ...
    & newSpectra(:, 3) <= options.Spk4noise;
if sum(IdSpikes) < options.MinMinNoise(1)
    minNoise = options.MinMinNoise(2);
    DoNoise = false;

else
    DoNoise = true;

end

ScNu = 0;
CurSeq = zeros(numel(short_mzAxis), 2);
NoiseVector = [];

for ii = 1:numel(newProfiles(:,1))
    ScanName = fullfile(obj.Path2Fin, newDataset,...
        'Scans', ['Scan#', num2str(newProfiles(ii,1)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [newProfiles(ii, 8) newProfiles(ii, 9)],...
        "double");
    fclose(fidRead);

    if isempty(myMS)
        cMZ = short_mzAxis; cMZ(:, 2) = 0;

    else
        [~, IdI] = intersect(short_mzAxis, myMS(:,1));
        cMZ = short_mzAxis; cMZ(IdI, 2) = myMS(:,2);
    end
    
    % Record Noise
    NoiseVector = [NoiseVector; cMZ(IdSpikes & cMZ(:, 2) > 0, 2)];

end
NoiseVector(isoutlier(NoiseVector)) = [];
NoiseVector(isoutlier(NoiseVector)) = [];
minNoise = round(mean(NoiseVector) + 2*std(NoiseVector));

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
        options.orderPolynomial = 6;
        options.parentDataset   = dts;
        options.AverageScans    = 1;
        options.Spk4noise       = 3;
        options.Sig2Nois        = 3;
        options.MinMinNoise     = [100, 10];
        options.ParallelMe      = false;

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

        tgtIx = input('AverageScans');
        if ~isempty(tgtIx)
            options.AverageScans = varargin{tgtIx +1};
        end
    end

%% Extrapolate master MZ axis
    function [mzAxis, r2, XY] = mkMasterMZ(MS, n, Lim, ratio)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here

        if nargin == 1
            Lim = infoDataset.mzInterval;
            n   = 3;
            ratio = 1;
        elseif nargin == 2
            Lim = infoDataset.mzInterval;
            ratio = 1;
        elseif nargin == 3
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
        [p, S, mu] = polyfit(XY(:,1), XY(:,2), n);
        XY(:, 3) = polyval(p, XY(:,1), S, mu);
        XY(:, 4) = sqrt((XY(:, 2) - XY(:, 3)).^2);
        XY(XY(:, 4) > mean(XY(:, 4)) + 5*std(XY(:, 4)), :) = [];

        [p, S, mu] = polyfit(XY(:,1), XY(:,2), n);
        XY(:, 3) = polyval(p, XY(:,1), S, mu);
        XY(:, 4) = sqrt((XY(:, 2) - XY(:, 3)).^2);

        R = corrcoef(polyval(p, XY(:,1), S, mu), XY(:,2));
        r2 = R(1,2);
        mzAxis = Lim(1, 1);

        int_old = 0;
        while mzAxis(end, 1) < Lim(2)
            int = polyval(p, mzAxis(end, 1), S, mu)/ratio;
            mzAxis(end+1, 1) = mzAxis(end, 1) + int; %#ok<AGROW>
            if int_old > int
                error("")
            else
                int_old = int;
            end
        end
    end


end