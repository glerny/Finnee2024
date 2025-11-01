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

function obj = mkDerivatives(obj, dts,varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options

narginchk(2, inf)
tStart = tic;

if nargin == 2 & isstruct(dts)
    options = dts;
    dts = 2;
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    
    options.ParallelMe = true;

else
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);

    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    options = checkVarargin(dts, varargin{:});
end


Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
if isempty(Id2Dataset), error(""), end
if obj.Datasets.HasMasterMZ(Id2Dataset), error(""), end
if ~obj.Datasets.isProfileScan{Id2Dataset}, error(""), end
if ~obj.Datasets.IsMS1{Id2Dataset}, error(""), end

fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
fclose(fidRead);
Profiles = array2table(Profiles);
ColumnNames = {}; ColumnUnits = {};
for ii = 1:numel(infoDataset.Profiles.column)
        ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
        ColumnLabels{ii} = infoDataset.Profiles.column{ii}.labelUnits;
        ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
        ColumnFO{ii}    = infoDataset.Profiles.column{ii}.formatUnits;
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

fileName_axis = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'MasterMZAxis.dat');
fidRead = fopen(fileName_axis, 'rb');
MasterMz = fread(fidRead, inf, "double");
fclose(fidRead);


%% INITIATE THE MATRIX FOR DERIVATIVE

%% 7. Recording new dataset and saving
% Create new dat file
ii = 1;
while 1
    newDataset = ['Dataset', num2str(ii)];
    if ~any(strcmp(newDataset, obj.Datasets.Name)), break; end
    ii = ii + 1;
end
mzInterval = [inf 0];
mkdir(fullfile(obj.Path2Fin, newDataset));
mkdir(fullfile(obj.Path2Fin, newDataset, 'Scans'));
newProfiles = [];
newSpectra = MasterMz; newSpectra(:, 7) = 0;

%TIP, BPP and more
CurSeq = zeros(numel(newSpectra(:,1)), 2);
iStart = 1; StopMe = false ;

Discard = max(max(options.smoothing), max(options.spreading, [], 'all'))+ 9;
firstLoop = true;

while 1
    iEnd = iStart + Discard + options.sizeMat-1;

    if iEnd > numel(Profiles.Scan_Index)
        iEnd = numel(Profiles.Scan_Index);
        StopMe = true;

    end

    if firstLoop
        cuMatrix = zeros(numel(MasterMz), iEnd-iStart);
        for ii = iStart:iEnd
            ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
            fidRead = fopen(ScanName, 'rb');
            myMS = fread(fidRead, ...
                [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                "double");
            fclose(fidRead);

            if ~isempty(myMS)
                [~, Id1] = intersect(MasterMz, myMS(:,1));
                if numel(Id1) < numel(myMS(:,1)), error(""); end

                cuMatrix(Id1, ii-iStart+1) = myMS(:, 2);

            end
        end
        firstLoop = false;

    else
        cuMatrix = cuMatrix(:, end-(2*Discard+1):end);
        cuMatrix(:, iEnd-iStart) = 0;

        for ii = iStart+(2*Discard+1):iEnd
            ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
            fidRead = fopen(ScanName, 'rb');
            myMS = fread(fidRead, ...
                [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                "double");
            fclose(fidRead);

            if ~isempty(myMS)
                [~, Id1] = intersect(MasterMz, myMS(:,1));
                if numel(Id1) < numel(myMS(:,1)), error(""); end

                cuMatrix(Id1, ii-iStart+1) = myMS(:, 2);

            end
        end

        

    end

    % Filter ions
    if options.FilterIons.do
        FilterOut = false(numel(MasterMz), 1); 

        for ii = 1:numel(options.FilterIons.Ions)
            Dmz = 4*(options.FilterIons.Ions(ii)/options.FilterIons.Dmz);
            IdX = MasterMz >= options.FilterIons.Ions(ii) - Dmz & ...
                MasterMz <= options.FilterIons.Ions(ii) + Dmz;
            FilterOut(IdX) = true;


        end
        cuMatrix(FilterOut, :) = 0;
    end

    % Remove spikes (1 & 2)
    for ii = 3:width(cuMatrix)-2
        IdSp = cuMatrix(:, ii-1) == 0 & cuMatrix(:, ii+1) == 0 |... 
            cuMatrix(:, ii-1) == 0 & cuMatrix(:, ii+2) == 0 |...
            cuMatrix(:, ii-2) == 0 & cuMatrix(:, ii+1) == 0;
        cuMatrix(IdSp, ii) = 0;

    end
    
    cy = options.cycles;
    SG = options.smoothing(1);
    Ga = options.spreading(1);

    if cy == 0
        % [~,g] = sgolay(2, SG);
        % Data2 = filter2(g(:, 1)', cuMatrix, 'same');
        Data2 = smoothdata2(cuMatrix, 'gaussian', {1, SG});

    else
        [~,g] = sgolay(2, SG);
        Data2 = cuMatrix;

        if Ga ~= 0
            Data2 = smoothdata2(Data2, 'gaussian', {1, Ga});
            Data2(isnan(Data2)) = 0;

        end
        Data2 = -filter2(g(:, 3)', Data2, 'same');
        Data2(isnan(Data2)) = 0;

        % %%%%%%%%%%%%%%%%%%%%%%%%PROVISORY&&&&&&&&&&&&&&&&&&&&&&&
        %  Ga = options.spreading(2);
        % 
        % if Ga > 0
        %     Data2 = smoothdata2(Data2, 'gaussian', {3, Ga});
        %     Data2(Data2 < 0) = 0;
        % 
        % end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    Data2(Data2 < 0) = 0;

    for ii = 2:width(Data2)-1 % Remove spikes
        Id2Z = Data2(:, ii-1) == 0 & Data2(:, ii+1) == 0;
        Data2(Id2Z, ii) = 0;

    end

    for jj = 2:cy
        SG = options.smoothing(jj);
        Ga = options.spreading(jj);

        if Ga > 0
            %! TO BE CHANGED TO !1D filter and optimise for speed
           Data2 = smoothdata2(Data2, 'gaussian', {1, Ga});
           Data2(isnan(Data2)) = 0;

        end
        [~,g] = sgolay(2, SG);
        Data2 = -filter2(g(:, 3)', Data2, 'same');
        Data2(isnan(Data2)) = 0;

        % %%%%%%%%%%%%%%%%%%%%%%%%PROVISORY&&&&&&&&&&&&&&&&&&&&&&&
        %  Ga = options.spreading(jj+1);
        % 
        % if Ga > 0
        %     Data2 = smoothdata2(Data2, 'gaussian', {3, Ga});
        %     Data2(Data2 < 0) = 0;
        % 
        % end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Data2(Data2 < 0) = 0;
        for ii = 2:width(Data2)-1
            Id2Z = Data2(:, ii-1) == 0 & Data2(:, ii+1) == 0;
            Data2(Id2Z, ii) = 0;
        end
    end
    
    for ii = iStart+Discard:iEnd-Discard

        MS = [MasterMz, Data2(:, ii-iStart+1)];

        % Filter spikes if needed
        if ~isempty(MS) && options.RemSpks
            spkSz = options.SpkSz;
            MS    = spikesRemoval(MS, spkSz );
        end

        %Record MZSpectra
        IdNnz = MS(:, 2) > 0;
        newSpectra(:,2) = newSpectra(:,2) + double(IdNnz);
        IdStart = CurSeq(:, 2) == 0 & IdNnz;

        CurSeq(IdStart, 1) = ii;
        CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

        % Check Length
        IdMatch = newSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
        newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
        newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
        CurSeq(~IdNnz, 2) = 0;
        newSpectra(:,5) = newSpectra(:,5) + MS(:, 2);
        newSpectra(:,6) = max([newSpectra(:,6), MS(:, 2)], [], 2);

        % reduced trailing zero in excess
        if ~isempty(MS)
            provMat      = [MS(2:end, 2); 0];
            provMat(:,2) = MS(:, 2);
            provMat(:,3) = [0; MS(1:end-1, 2)];
            MS         = MS(sum(provMat, 2) > 0, :);
        end

        if min(MS(:,1)) < mzInterval(1), mzInterval(1) = min(MS(:,1)); end
        if max(MS(:,1)) > mzInterval(2), mzInterval(2) = max(MS(:,1)); end

        if isempty(MS) | size(MS, 1) == 1 | sum(MS(:, 2)) == 0
            newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
                Profiles.Scan_Time(ii), ...
                Profiles.Injection_Time(ii), 0, 0, 0, ...
                nnz(MS(:,2)), size(MS), nan, nan];

        else

            M = ChrMoment(MS);
            [~, IdM] = max(MS(:, 2));
            IdM = min(max(IdM, 3), height(MS(:,1))-2);
            try
                mzApex = sum(MS(IdM-2:IdM+2,1).*MS(IdM-2:IdM+2,2))/sum(MS(IdM-2:IdM+2,2));

            catch
                mzApex = NaN;
            end

            newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
                Profiles.Scan_Time(ii), ...
                Profiles.Injection_Time(ii), max(MS(:,2)), sum(MS(:,2)), ...
                min(MS(MS(:, 2) ~= 0, 2)), nnz(MS(:,2)),...
                size(MS), M(2), mzApex];
        end

        fileName = fullfile(obj.Path2Fin, newDataset, 'Scans', ['Scan#', ...
            num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
        fwrite(fidWriteDat, MS(:), "double");
        fclose(fidWriteDat);
    end

    if StopMe
        break;
    end
    iStart = iEnd - 2*Discard +1;

end

IdMatch = newSpectra(:,3) < CurSeq(:, 2);
newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
newSpectra(:, 7) = newSpectra(:,5)./newSpectra(:,2);

provMat = [newSpectra(2:end, 5); 0];
provMat(:,2) = newSpectra(:, 5);
provMat(:,3) = [0; newSpectra(1:end-1, 2)];
newSpectra = newSpectra(sum(provMat, 2) > 0, :);
short_mzAxis =  newSpectra(:,1);

infoDataset.dateofcreation = datetime("now");
infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
infoDataset.Label = obj.Datasets.Labels{end};
infoDataset.mzInterval =  mzInterval;
infoDataset.Options = options;
infoDataset.Profiles.size = size(newProfiles);
infoDataset.Spectra.size = size(newSpectra);
infoDataset.NoiseVector = [];

infoDataset.Profiles.column{11}.name = 'Apexmz';
infoDataset.Profiles.column{11}.labelUnits = infoDataset.Profiles.column{10}.labelUnits;
infoDataset.Profiles.column{11}.units = infoDataset.Profiles.column{10}.units;
infoDataset.Profiles.column{11}.formatUnits = infoDataset.Profiles.column{10}.formatUnits;
infoDataset.Profiles.column{11}.ShowMe = true;


obj.Datasets = [obj.Datasets; ...
    {newDataset, '', true, true, true, true, false, ...
    ['Dataset', num2str(dts)], obj.Datasets.Labels{Id2Dataset}, ...
    options,  'Derivatives', {}}];

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

fileName_old = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'MasterMZAxis.dat');
fileName_new = fullfile(obj.Path2Fin, newDataset);
copyfile(fileName_old, fileName_new)

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
        options.RemSpks         = true;
        options.SpkSz           = 1;
        options.YLim            = infoDataset.mzInterval;
        % options.smoothing       = [11; 7; 5; 3; 3; 3];
        % options.spreading       = [0; 11; 7; 5; 3; 3];
        options.smoothing       = [3; 3; 3; 3; 3; 3; 3; 3];
        options.spreading       = [61; 41; 31; 21; 11; 9; 7; 5];

        % options.smoothing       = [3; 3; 3; 3; 3; 3; 3; 3];
        % options.spreading       = [1, 91; 1, 51; 1, 43; 1, 35; 1, 27; 1, 19; 1, 11; 1, 3];
        options.cycles          = 8; 
        options.ParallelMe      = false;
        options.sizeMat         = 500;
        options.FilterIons.do   = false;
        options.FilterIons.Dmz  = 50000;
        options.FilterIons.Ions  = [229.1428, 201.1125, 83.0603, 99.5309, ...
            129.0545, 158.9610, 90.5258, 111.0435, 246.1681, 88.0232, 125.9862, ...
            143.9970, 102.0336];

        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);

        tgtIx = input('cycles');
        if ~isempty(tgtIx)
            spreading         = varargin{tgtIx +1};
            options.cycles    = length(spreading);
            options.spreading = spreading;
        end

        tgtIx = input('Dmz');
        if ~isempty(tgtIx)
            options.FilterIons.Dmz = varargin{tgtIx +1};
        end

        tgtIx = input('ionsOut');
        if ~isempty(tgtIx)
            options.FilterIons.do = true;
            options.FilterIons.Ions = varargin{tgtIx +1};
        end
% 
%         tgtIx = input('spikes');
%         if ~isempty(tgtIx)
%             spks = varargin{tgtIx +1};
%             if spks == 0
%                 options.RemSpks = false;
%             else
%                 options.RemSpks = true;
%                 options.SpkSz  =  spks;
%             end
%         end
% 
%         tgtIx = input('TimeSpikes');
%         if ~isempty(tgtIx)
%             spks = varargin{tgtIx +1};
%             if spks == 0
%                 options.RemTmSpks = false;
%             else
%                 options.RemTmSpks = true;
%                 options.TmSpkSWd  =  spks;
%             end
%         end
% 
%         tgtIx = input('AverageScans');
%         if ~isempty(tgtIx)
%             options.AverageScans = varargin{tgtIx +1};
%         end
    end


end