%% DESCRIPTION
% the FILTERDATASET method is used to filter a dataset using a particular
% method and create a new dataset with the filtered scans. FILTERDATASET is
% used in particular to remove the spikes or the background noise. For more
% information see
% <https://github.com/glerny/Finnee2016/wiki/Basic-operation-with-Dataset>
%
%% INPUT PARAMETERS
% *Compulsory
%   _obj_     : The Finnee object
%   _dts_     : The indices to the dataset to correct
%       (i.e. myFinnee.Datasets{dts}
%   _method_  : The method to implement in the format
%   'methodName:param1:param2:...:paramn'
%       + 'RemoveNoise:pm1:pm2:pm3:pm4'
%           Use to remove background noise. For each points, will
%           reconstitued a small matrix, centered around the point of interest
%           and of size 2*pm1+1 in the time dimension and 2*pm2+1 in the
%           m/z dimension. The center point intensity will be set to zero
%           if no values within the matrix has a signal to noise higher than
%           pm3. The function only works in the backgroun noise has been
%           estimated, pm4 defined how the noise is calculated in m/z were
%           no data are available. If pm4 start with % (%5 or %25 for
%           example), missing dat are the mean of the percentil, otherwise
%           it is a set value.
%   - varargin: Possible options
%       + 'spikes:pm1'
%          Only actives with 'RemoveNoise:pm1:pm2:pm3'. Will also performe
%          a spikes removal filter after the noise removal filter.
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = filterDataset(obj, dts, method, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options
narginchk(2, inf)

if nargin == 2 & isstruct(dts)
    options = dts;
    dts = options.parentDataset;
    options.ParallelMe = true;
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    method = options.method;
    MtU = strsplit(method, ':');
    MtU{1} = 'remNoise';
    if length(MtU) < 2
        MtU{2} = 10;
    else
        MtU{2} = str2double(MtU{2});
    end

    if length(MtU) < 3
        MtU{3} = 5;
    else
        MtU{3} = str2double(MtU{3});
    end

    if length(MtU) < 4
        MtU{4} = 10;
    else
        MtU{4} = str2double(MtU{4});
    end

else
    Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
    infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
    infoDataset = infoDataset.infoDataset;
    [options, MtU]   = checkVarargin(method, dts, varargin{:});

    if nargin == 2, mzAxis = []; end
end

tStart = tic;
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
% 2- Load each scan and run the filter algorithm

% Create new dat file
ii = 1;
while 1
    newDataset = ['Dataset', num2str(ii)];
    if ~any(strcmp(newDataset, obj.Datasets.Name)), break; end
    ii = ii + 1;
end
infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;

mzInterval = [inf 0];
mkdir(fullfile(obj.Path2Fin, newDataset));
mkdir(fullfile(obj.Path2Fin, newDataset, 'Scans'));
newProfile = [];

switch MtU{1}
    case 'remNoise'
        %% REMNOISE HERE
        fileName_axis = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'MasterMZAxis.dat');
        fidRead = fopen(fileName_axis, 'rb');
        mzAxis = fread(fidRead, inf, "double");
        fclose(fidRead);
        minNoise = infoDataset.minNoise;
        NoiseVector = infoDataset.NoiseVector;
        [~, IdN] = intersect(mzAxis, Spectra.("m/z axis"));

        newId = [];
        for jj = -MtU{3}-2:MtU{3}+2
            newId = unique([newId; IdN+jj]);
        end
        
        newId( newId < 1 | newId > numel(mzAxis)) = [];

        %%% ATTENTION TO CHECK
        mzAxis(IdN, 2) = minNoise;%Spectra.Noise;
        mzAxis = mzAxis(newId, :);
        mzAxis(mzAxis(:, 2) < minNoise, 2) = minNoise;

        newSpectra = mzAxis(:, 1);
        newSpectra(:, 7) = 0;

        S2NMat  = zeros(size(mzAxis, 1), 2*MtU{2}+1);
        for ii  = 1:MtU{2}+1

            ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
            fidRead = fopen(ScanName, 'rb');
            myMS = fread(fidRead, ...
                [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                "double");
            fclose(fidRead);

            if ~isempty(myMS)
                [~, Id1] = intersect(mzAxis(:,1), myMS(:,1));
                if numel(Id1) < numel(myMS(:,1)), error(""); end
                fMS = zeros(size(mzAxis(:,1)));
                fMS(Id1) = myMS(:,2);
            else
                fMS = zeros(size(mzAxis(:,1)));
            end

            cNoise = mzAxis(:,2);
            Id2Dn = find(NoiseVector(:, 2) <= Profiles.Scan_Time(ii) & NoiseVector(:, 3) >= Profiles.Scan_Time(ii));

            if ~isempty(Id2Dn)
                [~, Id2] = intersect(mzAxis(:,1), NoiseVector(Id2Dn, 1));
                cNoise(Id2) = NoiseVector(Id2Dn, 4);
            end
            S2NMat(:, MtU{2}+ii)  = fMS./cNoise;

        end
        vector = any(S2NMat > MtU{4}, 2);
        vector = smooth(vector,  2*MtU{3}+1);

        CurSeq = zeros(numel(mzAxis(:, 1)), 2);
        if ~options.ParallelMe, h = waitbar(0,'processing scans'); end
        for ii = 1:length(Profiles.Scan_Index)
            if ~options.ParallelMe, waitbar(ii/length(Profiles.Scan_Index)); end

            % 2. Check in line and column
            ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
            fidRead = fopen(ScanName, 'rb');
            myMS = fread(fidRead, ...
                [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                "double");
            fclose(fidRead);

            if isempty(myMS)
                MS2Cor = mzAxis(:,1);
                MS2Cor(:, 2) = 0;

            else
                [~, Id1] = intersect(mzAxis(:,1), myMS(:,1));
                if numel(Id1) < numel(myMS(:,1)), error(""); end
                MS2Cor = mzAxis(:,1);
                MS2Cor(Id1, 2) = myMS(:,2);

            end

            MS2Cor(vector == 0, 2) = 0;

            % Filter spikes if needed
            if ~isempty(MS2Cor) && options.RemSpks
                spkSz = options.SpkSz;
                MS2Cor    = spikesRemoval(MS2Cor, spkSz );
            end

            %Record MZSpectra
            IdNnz = MS2Cor(:, 2) > 0;
            newSpectra(:,2) = newSpectra(:,2) + double(IdNnz);
            IdStart = CurSeq(:, 2) == 0 & IdNnz;

            CurSeq(IdStart, 1) = ii;
            CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

            % Check Length
            IdMatch = newSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
            newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
            newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
            CurSeq(~IdNnz, 2) = 0;
            newSpectra(:,5) = newSpectra(:,5) + MS2Cor(:, 2);
            newSpectra(:,6) = max([newSpectra(:,6), MS2Cor(:, 2)], [], 2);

            % reduced trailing zero in excess
            if ~isempty(MS2Cor)
                provMat      = [MS2Cor(2:end, 2); 0];
                provMat(:,2) = MS2Cor(:, 2);
                provMat(:,3) = [0; MS2Cor(1:end-1, 2)];
                MS2Cor         = MS2Cor(sum(provMat, 2) > 0, :);
            end

            if min(MS2Cor(:,1)) < mzInterval(1), mzInterval(1) = min(MS2Cor(:,1)); end
            if max(MS2Cor(:,1)) > mzInterval(2), mzInterval(2) = max(MS2Cor(:,1)); end

            if isempty(MS2Cor) | size(MS2Cor, 1) == 1 | sum(MS2Cor(:, 2)) == 0
                newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
                    Profiles.Scan_Time(ii), ...
                    Profiles.Injection_Time(ii), 0, 0, 0, ...
                    nnz(MS2Cor(:,2)), size(MS2Cor), nan];

            else

                M = ChrMoment(MS2Cor);
                newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
                    Profiles.Scan_Time(ii), ...
                    Profiles.Injection_Time(ii), max(MS2Cor(:,2)), sum(MS2Cor(:,2)), ...
                    min(MS2Cor(MS2Cor(:, 2) ~= 0, 2)), nnz(MS2Cor(:,2)),...
                    size(MS2Cor), M(2)];
            end

            fileName = fullfile(obj.Path2Fin, newDataset, 'Scans', ['Scan#', ...
                num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
            [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
            fwrite(fidWriteDat, MS2Cor(:), "double");
            fclose(fidWriteDat);

            S2A  = mzAxis(:,1); S2A(:, 2) = 0;
            cNoise = mzAxis(:,2);
            % Load the  next scan
            if ii+MtU{2}+1 > length(Profiles.Scan_Index)
                S2A(:, 2) = 0;

            else
                ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii+MtU{2}+1, 1)), '.dat']);
                fidRead = fopen(ScanName, 'rb');
                myMS = fread(fidRead, ...
                    [Profiles.Length_MS_Scan_Profile(ii+MtU{2}+1) Profiles.Column_MS_Scan_Profile(ii+MtU{2}+1)],...
                    "double");
                fclose(fidRead);

                if isempty(myMS)
                    S2A(:, 2) = 0;

                else
                    [~, Id1] = intersect(S2A(:,1), myMS(:,1));
                    S2A(Id1, 2) =  myMS(:,2);

                end

                Id2Dn = find(NoiseVector(:, 2) <= Profiles.Scan_Time(ii+MtU{2}+1)...
                    & NoiseVector(:, 3) >= Profiles.Scan_Time(ii+MtU{2}+1));

                if ~isempty(Id2Dn)
                    [~, Id2] = intersect(mzAxis(:,1), NoiseVector(Id2Dn, 1));
                    cNoise(Id2) = NoiseVector(Id2Dn, 4);
                end
            end

            S2NMat = [S2NMat(:, 2:end), S2A(:, 2)./cNoise];
            vector = any(S2NMat > MtU{4}, 2);
            vector = smooth(vector,  2*MtU{3}+1);

        end

    otherwise
        error('%s is not a recognised method', MtU{1})
end
try close(h), catch, end
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
infoDataset.AdditionalInformation.nonEndingErrors{1} = '';
infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
infoDataset.Label = obj.Datasets.Labels{end};
infoDataset.mzInterval =  mzInterval;
infoDataset.Options = options;
infoDataset.Profiles.size = size(newProfiles);
infoDataset.Spectra.size = size(newSpectra);

obj.Datasets = [obj.Datasets; ...
    {newDataset, '', true, true, true, true, true, ...
    ['Dataset', num2str(dts)], obj.Datasets.Labels{Id2Dataset}, ...
    options,  'Noise_removal', {}}];

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
    function [options, MtU] = checkVarargin(method, dts, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.

        %
        % 2- Default parameters
        options.RemSpks = true;
        options.SpkSz   = 3;
        options.parentDataset = dts;
        options.wdw4corr        = 2;
        options.ParallelMe = false;
        options.Spk4noise       = 4;
        options.Sig2Nois        = 3;

        % check parameter for method
        MtU = strsplit(method, ':');
        MtU{1} = 'remNoise';
        if length(MtU) < 2
            MtU{2} = 10;
        else
            MtU{2} = str2double(MtU{2});
        end

        if length(MtU) < 3
            MtU{3} = 5;
        else
            MtU{3} = str2double(MtU{3});
        end

        if length(MtU) < 4
            MtU{4} = 10;
        else
            MtU{4} = str2double(MtU{4});
        end
        options.method = ['RemoveNoise:', num2str(MtU{2}),':',...
            num2str(MtU{3}),':', num2str(MtU{4})];

        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        tgtIx = input('spikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemSpks = false;
            else
                options.RemSpks = true;
                options.SpkSz   =  spks;
            end
        end
    end
end

