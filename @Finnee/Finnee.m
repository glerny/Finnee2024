%% DESCRIPTION
% CFINNEE is the entry object of the cFinnee toolbox.

%% LIST OF THE CLASS'S PROPERTIES

%% LIST OF THE CLASS'S METHODS

%% REFERENCES

%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2022-2023 G. Erny
% (guillaume.erny@outlook.com)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Finnee

    properties
        FileID         % The generic name
        FileIn         % Path to the original file (Single only)
        DateOfDataAcq  % The date of data acquisition
        DateOfCreation % The date of creation of the cFinnee object
        Datasets       % The array of datasets linked to this object
        Options        % The options that were used by the constructor
        % method
        Path2Fin       % The paths to the .fin folder
        Type           % "Single", "Replicate" or "Multiple"
        FinneesIn       % Paths to the Finnee files used for the "Replicate"
        % or "Multiple"
        Version = '6.01';

    end

    methods
        function obj = Finnee(Type, varargin)
            %% DESCRIPTION
            % FINNEE is the constructor method. The following options can
            % be used when creating a cFinnee object:
            %
            % *fileIn    For Single: Followed by a string that is the full
            %                        path of the mzML file.
            % *folderOut Followed by a tring that is the path to the
            %            destination folder for this FINNEE object.
            % *fileID    Followed by a string that is the generic name for
            %            this object. the Finnee objects and all associated
            %            data files will be recorded in the folder
            %            'folderOut\fileID.fin'.
            % *overwrite delete the folder 'folderOut\fileID.fin' and ALL
            %            files in this folder if it already exits.
            % *tLim      Followed by a 2x1 array of numbers
            %            (default [0 inf]). Only records scans between
            %            tLim(1) and tLim(2)
            % *mzLim     Followed by a 2x1 array of numbers
            %            (default [0 inf]). For each scan only kept mz
            %            values between mzLim(1) and mzLim(2).
            % *spikes    Followed by an integer between 0 and 10 (default
            %            3). Remove spikes in every MS scans If used, where
            %            spikes are any peaks in MS scans, delimited by
            %            double zeros, of length equal or lower that the
            %            integer. 'spikes' followed by 0 allows to to turn
            %            off spikes removal.
            % *MasterMz  followed by a structure that contains as principal
            %            field .Axis, the master mz axis as well as
            %            specific oprions (see... for more info)
            % *MasterTm  followed by a structure that contains as principal
            %            field .Axis, the master time axis as well as
            %            specific options (see... for more info)
            % *Fast      Assume MS1 data with only one type of recording.
            %
            %% EXAMPLES
            % * myFinnee = cFinnee;
            %   Will create a Finnee object without any options. The target
            %   mzML file, destination folder and generic name will be
            %   asked when running the script.
            % *	myFinnee = Finnee('fileIn'   , 'K:\Data\test.mzML', ...
            %                     'folderOut', 'K:Finnee',          ...
            %                     'Overwrite',                      ...
            %                     'fileID'   , 'test1',             ...
            %                     'Spikes'   , 0);
            %   Here all parameters are defined by the optional parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% CORE OF THE FUNCTION

            % 1- Initialisation and options

            if nargin == 0, Type = 'Single'; end
            if isempty(Type), Type = 'Single'; end

            switch Type

                %% CASE SINGLE %%
                case 'Single'

                    % Check the options and create the Finnee object
                    options = checkVarargin_single(varargin{:});
                    % Initialization of Finnee
                    obj.DateOfCreation = datetime;
                    obj.FileID         = options.FileID;
                    obj.Datasets       = {};
                    obj.Options        = options;
                    obj.Path2Fin       = '';
                    obj.Type           = Type;
                    recogniseFormats   = {'.mzML'};

                    % Check Options, ask for mzML file/folder and name if the
                    % corresponding options are empty
                    if isempty(options.FileIn)
                        ext = 'pwd/*.mzML';
                        txtStg = 'Select the file to load';
                        [fileName, pathName] = uigetfile(ext, txtStg);
                        if ~ischar(fileName) && ~ischar(pathName)
                            error('myApp:argChk', 'User cancel file selection');

                        end
                        obj.FileIn = fullfile(pathName, fileName);
                    else

                        obj.FileIn = options.FileIn;
                    end

                    if isempty(obj.FileID)
                        [~, dfltName, ~] = fileparts(obj.FileIn);
                        answer = inputdlg('Enter a name for the destination folder', ...
                            'fileID', 1, {dfltName});
                        if isempty(answer)
                            error('myApp:argChk', 'Cancel by user')

                        elseif isempty(answer{1})
                            error('myApp:argChk', 'Cancel by user')

                        end
                        obj.FileID = answer{1};

                    end

                    if isempty(options.FolderOut)
                        options.FolderOut = uigetdir(pwd, 'Select the folder of destination');
                        if ~ischar(options.FolderOut)

                            error('myApp:argChk', 'Cancel by user');
                        end
                    end

                    if exist(fullfile(options.FolderOut, [obj.FileID '.fin']), 'dir') == 7
                        if options.Overwrite
                            rmdir(fullfile(options.FolderOut, [obj.FileID '.fin']), 's')
                        else
                            error('The directory %s already exist. \nDelete it, change the name or use the option ''Overwrite'' \ntype ''help cFinnee'' or ''help cFinnee\\cFinnee'' for more information', ...
                                fullfile(options.FolderOut, [obj.FileID '.fin']));
                        end
                    end

                    % Create the .fin folder
                    obj.Path2Fin = fullfile(options.FolderOut, [obj.FileID '.fin']);
                    try
                        mkdir(obj.Path2Fin);

                    catch ME
                        rethrow(ME)

                    end

                    [~, ~, ext] = fileparts(obj.FileIn);
                    switch ext
                        % Script for each format are in subfunction
                        % TODO: Other format

                        case '.mzML'
                            % Scripts for other format could be implemented here if needed
                            obj = mk_mzML2cFinnee(obj, options);

                        otherwise
                            error('Error. \n''%s'' is not a recognised extension.', ext)

                    end



                    % 2- Check and Run docFinnee
                    % Scripts for other format could be implemented here if needed


                    % MSn are the MSn scans that were present in the mzML file,
                    % they are not used yet; but are accessible via the following
                    % variable that will be created in the workspace.

                case 'Multiple'
                    % Check the options and create the Finnee object
                    options = checkVarargin_multiple(varargin{:});
                    % Initialization of Finnee
                    obj.DateOfCreation = datetime;
                    obj.FileID         = options.FileID;
                    obj.Datasets       = {};
                    obj.Options        = options;
                    obj.Path2Fin       = '';
                    obj.Type           = Type;
                    obj.FinneesIn      = options.FinneesIn;
                    targetDataset      = options.TgtDatasets;

                    tStart = tic;
                    if isempty(obj.FileID)
                        answer = inputdlg('Enter a name for the destination folder', ...
                            'fileID', 1, {'Merged_'});
                        if isempty(answer)
                            error('myApp:argChk', 'Cancel by user')

                        elseif isempty(answer{1})
                            error('myApp:argChk', 'Cancel by user')

                        end
                        obj.FileID = answer{1};

                    end

                    if isempty(options.FolderOut)
                        options.FolderOut = uigetdir(pwd, 'Select the folder of destination');
                        if ~ischar(options.FolderOut)

                            error('myApp:argChk', 'Cancel by user');
                        end
                    end

                    if isempty(options.FinneesIn)
                        obj.FinneesIn = uigetdirs();
                    end

                    if exist(fullfile(options.FolderOut, [obj.FileID '.fin']), 'dir') == 7
                        if options.Overwrite
                            rmdir(fullfile(options.FolderOut, [obj.FileID '.fin']), 's')

                        else
                            error('The directory %s already exist. \nDelete it, change the name or use the option ''Overwrite'' \ntype ''help cFinnee'' or ''help cFinnee\\cFinnee'' for more information', ...
                                fullfile(options.FolderOut, [obj.FileID '.fin']));
                        end
                    end

                    % Create the .fin folder
                    obj.Path2Fin = fullfile(options.FolderOut, [obj.FileID '.fin']);
                    try
                        mkdir(obj.Path2Fin);

                    catch ME
                        rethrow(ME)

                    end

                    % 1. Do the first object
                    load(fullfile(obj.FinneesIn{1}, 'myFinnee.mat'), 'myFinnee');
                    fObj = myFinnee;
                    Id2Dataset = find(strcmp(['Dataset', num2str(targetDataset)], fObj.Datasets.Name));
                    infoDataset = load(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
                    infoDataset = infoDataset.infoDataset;

                    fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
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

                    if ~isempty(options.Alignment)
                        Profiles.Scan_Time = Profiles.Scan_Time +...
                            polyval(options.Alignment{1}{1}, Profiles.Scan_Time, options.Alignment{1}{2}, options.Alignment{1}{3});

                    end
                    MasterTm = Profiles.Scan_Time;

                    fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Spectra.dat'), 'rb');
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

                    fileName_axis = fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'MasterMZAxis.dat');
                    fidRead = fopen(fileName_axis, 'rb');
                    MasterMz = fread(fidRead, inf, "double");
                    fclose(fidRead);


                    myDatasets.Name{1, 1} = 'Dataset1';
                    myDatasets.CreationCode{1, 1} = '';
                    myDatasets.IsMS1{1, 1} = true;
                    myDatasets.isProfileScan{1, 1} = true;
                    myDatasets.HasMasterMZ(1, 1) = true;
                    myDatasets.isBaseDriftCorr{1, 1} = false;
                    myDatasets.HasNoiseRem{1, 1} = false;
                    myDatasets.IsDaugherOff{1, 1} = '';
                    myDatasets.Labels{1, 1}.AxisX.Label = 'Time';
                    myDatasets.Labels{1, 1}.AxisX.Units = 'min'; % TODO Check Uniys
                    myDatasets.Labels{1, 1}.AxisX.fo = '%0.2f';% TODO Check precision
                    myDatasets.Labels{1, 1}.AxisY.Label = 'm/z';
                    myDatasets.Labels{1, 1}.AxisY.Units = '';
                    myDatasets.Labels{1, 1}.AxisY.fo = '%0.4f';
                    myDatasets.Labels{1, 1}.AxisZ.Label = 'Intensity';
                    myDatasets.Labels{1, 1}.AxisZ.Units = 'Arb. units';
                    myDatasets.Labels{1, 1}.AxisZ.fo = '%0.0f';
                    myDatasets.Options4Creations{1, 1} = options;
                    myDatasets.PrimaryActions{1, 1} = 'Creation';
                    myDatasets.SecondaryActions{1, 1} = {};

                    mkdir(fullfile(obj.Path2Fin, myDatasets.Name{1, 1}));
                    mkdir(fullfile(obj.Path2Fin, myDatasets.Name{1, 1}, 'Scans'));
                    mzInterval = [inf 0];
                    newProfiles = [];

                    ScanNumber = 0;
                    vector = [0 0];
                    vq = [];

                    for ii = 1:numel(MasterTm)
                        if any(Profiles.Scan_Time == MasterTm(ii))
                            IdM = find(Profiles.Scan_Time == MasterTm(ii));
                            ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdM, 1)), '.dat']);
                            fidRead = fopen(ScanName, 'rb');
                            cMS = fread(fidRead, ...
                                [Profiles.Length_MS_Scan_Profile(IdM) Profiles.Column_MS_Scan_Profile(IdM)],...
                                "double");
                            fclose(fidRead);
                            vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            vq(~isfinite(vq)) = 0;
                            myMS = [MasterMz, vq];
                            vector = IdM;

                        else
                            vq_old = vq;

                            IdS = find(Profiles.Scan_Time < MasterTm(ii), 1, 'last');
                            IdE = find(Profiles.Scan_Time > MasterTm(ii), 1, 'first');
                            if isempty(IdS), IdS = IdE; IdE = IdE+1; end
                            if isempty(IdE), IdE = IdS; IdS = IdS-1; end

                            if any(vector == IdS)
                                try
                                    vq = vq_old(:, vector == IdS);
                                catch
                                    disp("pp")
                                end
                            else
                                ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS, 1)), '.dat']);
                                fidRead = fopen(ScanName, 'rb');
                                cMS = fread(fidRead, ...
                                    [Profiles.Length_MS_Scan_Profile(IdS) Profiles.Column_MS_Scan_Profile(IdS)],...
                                    "double");
                                fclose(fidRead);
                                vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            end

                            if  any(vector == IdE)
                                vq(:,2) = vq_old(:, vector == IdE);
                            else
                                ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdE, 1)), '.dat']);
                                fidRead = fopen(ScanName, 'rb');
                                cMS = fread(fidRead, ...
                                    [Profiles.Length_MS_Scan_Profile(IdE) Profiles.Column_MS_Scan_Profile(IdE)],...
                                    "double");
                                fclose(fidRead);
                                vq(:,2) =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            end

                            vq(~isfinite(vq)) = 0;
                            vq2D = interp1(Profiles.Scan_Time(IdS:IdE), vq', MasterTm(ii));
                            vq2D(~isfinite(vq2D)) = 0;
                            myMS = [MasterMz, vq2D'];
                            vector = [IdS IdE];
                        end

                        % Filter spikes if needed
                        if options.RemSpks
                            spkSz = options.SpkSz;
                            myMS    = spikesRemoval(myMS, spkSz );
                        end

                        % reduced trailing zero in excess
                        if ~isempty(myMS)
                            provMat      = [myMS(2:end, 2); 0];
                            provMat(:,2) = myMS(:, 2);
                            provMat(:,3) = [0; myMS(1:end-1, 2)];
                            myMS         = myMS(sum(provMat, 2) > 0, :);
                        end

                        if isempty(myMS)
                            newProfiles(ScanNumber+1, :) = [ScanNumber, ...
                                MasterTm(ii), ...
                                nan, 0, sum(myMS(:,2)), 0, ...
                                nnz(myMS(:,2)), size(myMS), nan];

                        else

                            M = ChrMoment(myMS, 2);
                            newProfiles(ScanNumber+1, :) = [ScanNumber, ...
                                MasterTm(ii), ...
                                nan, max(myMS(:,2)), sum(myMS(:,2)), ...
                                min(myMS(myMS(:, 2) ~= 0, 2)), nnz(myMS(:,2)), ...
                                size(myMS), M(2)];

                        end

                        fileName = fullfile(obj.Path2Fin, 'Dataset1', 'Scans', ['Scan#', ...
                            num2str(ScanNumber), '.dat']);
                        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                        fwrite(fidWriteDat, myMS(:), "double");
                        fclose(fidWriteDat);
                        ScanNumber = ScanNumber + 1;
                    end

                    % Middle files
                    for jj = 2:numel(obj.FinneesIn)-1

                        load(fullfile(obj.FinneesIn{jj}, 'myFinnee.mat'), 'myFinnee');
                        fObj = myFinnee;

                        Id2Dataset = find(strcmp(['Dataset', num2str(targetDataset)], fObj.Datasets.Name));
                        infoDataset = load(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
                        infoDataset = infoDataset.infoDataset;

                        fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
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
                        if ~isempty(options.Alignment)
                            Profiles.Scan_Time = Profiles.Scan_Time +...
                                polyval(options.Alignment{jj}{1}, Profiles.Scan_Time, options.Alignment{jj}{2}, options.Alignment{jj}{3});

                        end

                        fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Spectra.dat'), 'rb');
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

                        vector = [0 0];
                        vq = [];
                        for ii = 1:numel(MasterTm)
                            vq_old = vq;
                            % Open old MS
                            ScanName = fullfile(obj.Path2Fin, 'Dataset1',...
                                'Scans', ['Scan#', num2str(newProfiles(ii, 1)), '.dat']);
                            fidRead = fopen(ScanName, 'rb');
                            cMS = fread(fidRead, ...
                                [newProfiles(ii, 8) newProfiles(ii, 9)],...
                                "double");
                            fclose(fidRead);

                            if ~isempty(cMS)
                                vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                                vq(~isfinite(vq)) = 0;
                                oldMS = [MasterMz, vq];

                            else
                                oldMS = MasterMz; oldMS(:, 2) = 0;

                            end

                            if any(Profiles.Scan_Time == MasterTm(ii))
                                IdM = find(Profiles.Scan_Time == MasterTm(ii));
                                ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdM, 1)), '.dat']);
                                fidRead = fopen(ScanName, 'rb');
                                cMS = fread(fidRead, ...
                                    [Profiles.Length_MS_Scan_Profile(IdM) Profiles.Column_MS_Scan_Profile(IdM)],...
                                    "double");
                                fclose(fidRead);
                                vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                                vq(~isfinite(vq)) = 0;
                                myMS = [MasterMz, vq];
                                vector = IdM;

                            else

                                IdS = find(Profiles.Scan_Time < MasterTm(ii), 1, 'last');
                                IdE = find(Profiles.Scan_Time > MasterTm(ii), 1, 'first');
                                if isempty(IdS), IdS = IdE; IdE = IdE+1; end
                                if isempty(IdE), IdE = IdS; IdS = IdS-1; end

                                if any(vector == IdS)
                                    vq = vq_old(:, vector == IdS);
                                else
                                    ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                        'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS, 1)), '.dat']);
                                    fidRead = fopen(ScanName, 'rb');
                                    cMS = fread(fidRead, ...
                                        [Profiles.Length_MS_Scan_Profile(IdS) Profiles.Column_MS_Scan_Profile(IdS)],...
                                        "double");
                                    fclose(fidRead);
                                    vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                                end

                                if  any(vector == IdE)
                                    vq(:,2) = vq_old(:, vector == IdE);

                                else
                                    ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                        'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdE, 1)), '.dat']);
                                    fidRead = fopen(ScanName, 'rb');
                                    cMS = fread(fidRead, ...
                                        [Profiles.Length_MS_Scan_Profile(IdE) Profiles.Column_MS_Scan_Profile(IdE)],...
                                        "double");
                                    fclose(fidRead);
                                    vq(:,2) =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                                end

                                vq(~isfinite(vq)) = 0;

                                vq2D = interp1(Profiles.Scan_Time(IdS:IdE), vq', MasterTm(ii));
                                vq2D(~isfinite(vq2D)) = 0;

                                myMS = [MasterMz, vq2D'];

                                vector = [IdS IdE];
                            end
                            myMS(:, 2) = myMS(:, 2) +  oldMS(:, 2);

                            % Filter spikes if needed
                            if options.RemSpks
                                spkSz = options.SpkSz;
                                myMS    = spikesRemoval(myMS, spkSz );
                            end

                            % reduced trailing zero in excess
                            if ~isempty(myMS)
                                provMat      = [myMS(2:end, 2); 0];
                                provMat(:,2) = myMS(:, 2);
                                provMat(:,3) = [0; myMS(1:end-1, 2)];
                                myMS         = myMS(sum(provMat, 2) > 0, :);
                            end

                            if isempty(myMS)
                                newProfiles(ii, :) = [ii-1, MasterTm(ii), ...
                                    nan, 0, sum(myMS(:,2)), 0, ...
                                    nnz(myMS(:,2)), size(myMS), nan];

                            else

                                M = ChrMoment(myMS, 2);
                                newProfiles(ii, :) = [ii-1, MasterTm(ii), ...
                                    nan, max(myMS(:,2)), sum(myMS(:,2)), ...
                                    min(myMS(myMS(:, 2) ~= 0, 2)), ...
                                    nnz(myMS(:,2)), size(myMS), M(2)];

                            end

                            fileName = fullfile(obj.Path2Fin, 'Dataset1', 'Scans', ['Scan#', ...
                                num2str(ii-1), '.dat']);
                            [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                            fwrite(fidWriteDat, myMS(:), "double");
                            fclose(fidWriteDat);
                        end

                    end

                    % last files

                    load(fullfile(obj.FinneesIn{end}, 'myFinnee.mat'), 'myFinnee');
                    fObj = myFinnee;

                    Id2Dataset = find(strcmp(['Dataset', num2str(targetDataset)], fObj.Datasets.Name));
                    infoDataset = load(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
                    infoDataset = infoDataset.infoDataset;

                    fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
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
                    if ~isempty(options.Alignment)
                        Profiles.Scan_Time = Profiles.Scan_Time +...
                            polyval(options.Alignment{end}{1}, Profiles.Scan_Time, options.Alignment{end}{2}, options.Alignment{end}{3});

                    end

                    fidRead = fopen(fullfile(fObj.Path2Fin, fObj.Datasets.Name{Id2Dataset}, 'Spectra.dat'), 'rb');
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

                    CurSeq = zeros(numel(MasterMz), 2);
                    NewSpectra = MasterMz; NewSpectra(:, 7) = 0;
                    vector = [0 0];
                    vq = [];
                    for ii = 1:numel(MasterTm)
                        vq_old = vq;
                        % Open old MS
                        ScanName = fullfile(obj.Path2Fin, 'Dataset1',...
                            'Scans', ['Scan#', num2str(newProfiles(ii, 1)), '.dat']);
                        fidRead = fopen(ScanName, 'rb');
                        cMS = fread(fidRead, ...
                            [newProfiles(ii, 8) newProfiles(ii, 9)],...
                            "double");
                        fclose(fidRead);

                        if ~isempty(cMS)
                            vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            vq(~isfinite(vq)) = 0;
                            oldMS = [MasterMz, vq];

                        else
                            oldMS = MasterMz; oldMS(:, 2) = 0;

                        end

                        if any(Profiles.Scan_Time == MasterTm(ii))
                            IdM = find(Profiles.Scan_Time == MasterTm(ii));
                            ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdM, 1)), '.dat']);
                            fidRead = fopen(ScanName, 'rb');
                            cMS = fread(fidRead, ...
                                [Profiles.Length_MS_Scan_Profile(IdM) Profiles.Column_MS_Scan_Profile(IdM)],...
                                "double");
                            fclose(fidRead);
                            vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            vq(~isfinite(vq)) = 0;
                            myMS = [MasterMz, vq];

                            vector = IdM;

                        else

                            IdS = find(Profiles.Scan_Time < MasterTm(ii), 1, 'last');
                            IdE = find(Profiles.Scan_Time > MasterTm(ii), 1, 'first');
                            if isempty(IdS), IdS = IdE; IdE = IdE+1; end
                            if isempty(IdE), IdE = IdS; IdS = IdS-1; end

                            if any(vector == IdS) & size(vector, 2) == size(vq_old, 1)
                                vq = vq_old(:, vector == IdS);

                            else
                                ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdS, 1)), '.dat']);
                                fidRead = fopen(ScanName, 'rb');
                                cMS = fread(fidRead, ...
                                    [Profiles.Length_MS_Scan_Profile(IdS) Profiles.Column_MS_Scan_Profile(IdS)],...
                                    "double");
                                fclose(fidRead);
                                vq =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            end

                            if  any(vector == IdE)
                                vq(:,2) = vq_old(:, vector == IdE);
                            else
                                ScanName = fullfile(fObj.Path2Fin, ['Dataset', num2str(Id2Dataset)],...
                                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(IdE, 1)), '.dat']);
                                fidRead = fopen(ScanName, 'rb');
                                cMS = fread(fidRead, ...
                                    [Profiles.Length_MS_Scan_Profile(IdE) Profiles.Column_MS_Scan_Profile(IdE)],...
                                    "double");
                                fclose(fidRead);
                                vq(:,2) =  interp1(cMS(:,1), cMS(:,2), MasterMz);
                            end

                            vq(~isfinite(vq)) = 0;
                            vq2D = interp1(Profiles.Scan_Time(IdS:IdE), vq', MasterTm(ii));
                            vq2D(~isfinite(vq2D)) = 0;
                            myMS = [MasterMz, vq2D'];
                        end
                        myMS(:, 2) = (myMS(:, 2) +  oldMS(:, 2))/numel(obj.FinneesIn);

                        % Filter spikes if needed
                        if options.RemSpks
                            spkSz = options.SpkSz;
                            myMS    = spikesRemoval(myMS, spkSz );
                        end

                        %Record MZSpectra
                        IdNnz = myMS(:, 2) > 0;
                        NewSpectra(:,2) = NewSpectra(:,2) + double(IdNnz);
                        IdStart = CurSeq(:, 2) == 0 & IdNnz;
                        CurSeq(IdStart, 1) = ScanNumber+1;
                        CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

                        % Check Length
                        % newSpectra(~IdNnz,3) = max(newSpectra(~IdNnz,3), CurSeq(~IdNnz, 2));
                        IdMatch = NewSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
                        NewSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
                        NewSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
                        CurSeq(~IdNnz, 2) = 0;
                        NewSpectra(:,5) = NewSpectra(:,5) + myMS(:, 2);
                        NewSpectra(:,6) = max([NewSpectra(:,6), myMS(:, 2)], [], 2);
                        %

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
                            newProfiles(ii, 1:10) = [ii-1,...
                                MasterTm(ii), nan, sum(myMS(:,2)), 0, ...
                                nnz(myMS(:,2)), size(myMS), NaN];
                        else

                            M = ChrMoment(myMS, 2);
                            newProfiles(ii, 1:10) = [ii-1, ...
                                MasterTm(ii), ...
                                nan,  max(myMS(:, 2)), sum(myMS(:, 2)),...
                                min(myMS(myMS(:, 2) ~= 0, 2)), ...
                                nnz(myMS(:,2)), size(myMS), M(2)];

                        end

                        fileName = fullfile(obj.Path2Fin, 'Dataset1', 'Scans', ['Scan#', ...
                            num2str(ii-1), '.dat']);
                        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                        fwrite(fidWriteDat, myMS(:), "double");
                        fclose(fidWriteDat);

                        vector = [IdS IdE];

                    end

                    IdMatch = NewSpectra(:,3) < CurSeq(:, 2);
                    NewSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
                    NewSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
                    NewSpectra(:,7) = NewSpectra(:,5)./NewSpectra(:,2);

                    provMat = [NewSpectra(2:end, 4); 0];
                    provMat(:,2) = NewSpectra(:, 4);
                    provMat(:,3) = [0; NewSpectra(1:end-1, 2)];
                    NewSpectra = NewSpectra(sum(provMat, 2) > 0, :);
                    short_mzAxis =  NewSpectra(:,1);

                    % Find mz lines with only spikes
                    IdSpikes = ~isnan(NewSpectra(:, 6)) ...
                        & NewSpectra(:, 3) <= options.Spk4noise;
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
                        ScanName = fullfile(obj.Path2Fin, 'Dataset1',...
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

                    obj.Options = options;
                    % 3- Save and done

                    infoDataset.dateofcreation = datetime("now");
                    infoDataset.AdditionalInformation.nonEndingErrors{1} = '';
                    infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
                    infoDataset.Label = myDatasets.Labels{1};
                    infoDataset.mzInterval =  mzInterval;
                    infoDataset.minNoise =  nan;
                    infoDataset.Options = obj.Options;
                    infoDataset.Profiles.size = size(newProfiles);
                    infoDataset.Profiles.column{1}.name = 'Scan_Index';
                    infoDataset.Profiles.column{1}.labelUnits = '#';
                    infoDataset.Profiles.column{1}.units = '';
                    infoDataset.Profiles.column{1}.formatUnits = '%0.0f';
                    infoDataset.Profiles.column{1}.ShowMe = false;
                    infoDataset.Profiles.column{2}.name = 'Scan_Time';
                    infoDataset.Profiles.column{2}.labelUnits = myDatasets.Labels{1}.AxisX.Label;
                    infoDataset.Profiles.column{2}.units = myDatasets.Labels{1}.AxisX.Units;
                    infoDataset.Profiles.column{2}.formatUnits = myDatasets.Labels{1}.AxisX.fo;
                    infoDataset.Profiles.column{2}.ShowMe = false;
                    infoDataset.Profiles.column{3}.name = 'Injection_Time';
                    infoDataset.Profiles.column{3}.labelUnits = 'Time';
                    infoDataset.Profiles.column{3}.units = '';
                    infoDataset.Profiles.column{3}.formatUnits = myDatasets.Labels{1}.AxisX.fo;
                    infoDataset.Profiles.column{3}.ShowMe = true;
                    infoDataset.Profiles.column{4}.name = 'Base_Peak_Profile';
                    infoDataset.Profiles.column{4}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Profiles.column{4}.units = myDatasets.Labels{1}.AxisZ.Units;
                    infoDataset.Profiles.column{4}.formatUnits = myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Profiles.column{4}.ShowMe = true;
                    infoDataset.Profiles.column{5}.name = 'Total_Ions_Profile';
                    infoDataset.Profiles.column{5}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Profiles.column{5}.units = myDatasets.Labels{1}.AxisZ.Units;
                    infoDataset.Profiles.column{5}.formatUnits = myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Profiles.column{5}.ShowMe = true;
                    infoDataset.Profiles.column{6}.name = 'Lowest_Intensity_Profile';
                    infoDataset.Profiles.column{6}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Profiles.column{6}.units = myDatasets.Labels{1}.AxisZ.Units;
                    infoDataset.Profiles.column{6}.formatUnits = myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Profiles.column{6}.ShowMe = false;
                    infoDataset.Profiles.column{7}.name = 'Non_Null_Elements_Profile';
                    infoDataset.Profiles.column{7}.labelUnits = '';
                    infoDataset.Profiles.column{7}.units = '';
                    infoDataset.Profiles.column{7}.formatUnits = '%0.0f';
                    infoDataset.Profiles.column{7}.ShowMe = false;
                    infoDataset.Profiles.column{8}.name = 'Length_MS_Scan_Profile';
                    infoDataset.Profiles.column{8}.labelUnits = '';
                    infoDataset.Profiles.column{8}.units = '';
                    infoDataset.Profiles.column{8}.formatUnits = '%0.0f';
                    infoDataset.Profiles.column{8}.ShowMe = false;
                    infoDataset.Profiles.column{9}.name = 'Column_MS_Scan_Profile';
                    infoDataset.Profiles.column{9}.labelUnits = '';
                    infoDataset.Profiles.column{9}.units = '';
                    infoDataset.Profiles.column{9}.formatUnits = '%0.0f';
                    infoDataset.Profiles.column{9}.ShowMe = false;
                    infoDataset.Profiles.column{10}.name = 'MS_centroid_profile';
                    infoDataset.Profiles.column{10}.labelUnits = myDatasets.Labels{1}.AxisY.fo;
                    infoDataset.Profiles.column{10}.units = myDatasets.Labels{1}.AxisY.Units;
                    infoDataset.Profiles.column{10}.formatUnits = myDatasets.Labels{1}.AxisY.fo;
                    infoDataset.Profiles.column{10}.ShowMe = true
                    
                    infoDataset.Spectra.size = size(NewSpectra);
                    infoDataset.Spectra.column{1}.name = 'm/z axis';
                    infoDataset.Spectra.column{1}.labelUnits = myDatasets.Labels{1}.AxisY.Label;
                    infoDataset.Spectra.column{1}.units = '';
                    infoDataset.Spectra.column{1}.formatUnits = myDatasets.Labels{1}.AxisY.fo;
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
                    infoDataset.Spectra.column{5}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Spectra.column{5}.units =  myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Spectra.column{5}.formatUnits =  myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Spectra.column{5}.ShowMe = true;
                    infoDataset.Spectra.column{6}.name = 'Base_Intensity_Spectrum';
                    infoDataset.Spectra.column{6}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Spectra.column{6}.units = myDatasets.Labels{1}.AxisZ.Units;
                    infoDataset.Spectra.column{6}.formatUnits = myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Spectra.column{6}.ShowMe = true;
                    infoDataset.Spectra.column{7}.name = 'Average_Intensity';
                    infoDataset.Spectra.column{7}.labelUnits = myDatasets.Labels{1}.AxisZ.Label;
                    infoDataset.Spectra.column{7}.units = myDatasets.Labels{1}.AxisZ.Units;
                    infoDataset.Spectra.column{7}.formatUnits = myDatasets.Labels{1}.AxisZ.fo;
                    infoDataset.Spectra.column{7}.ShowMe = true;


                    save(fullfile(obj.Path2Fin, 'Dataset1', 'infoDataset.mat'), 'infoDataset');
                    %Record Spectra
                    fileName = fullfile(obj.Path2Fin, 'Dataset1', 'Profiles.dat');
                    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                    fwrite(fidWriteDat, newProfiles(:), "double");
                    fclose(fidWriteDat);

                    fileName = fullfile(obj.Path2Fin, 'Dataset1', 'Spectra.dat');
                    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                    fwrite(fidWriteDat, NewSpectra(:), "double");
                    fclose(fidWriteDat);

                    fileName = fullfile(obj.Path2Fin, 'Dataset1', 'MasterMZAxis.dat');
                    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                    fwrite(fidWriteDat, MasterMz, "double");
                    fclose(fidWriteDat);
                    myDatasets = struct2table(myDatasets);
                    obj.Datasets = myDatasets;

                case 'Replicate'
                    %% CASE REPLICATE %%
            end

            obj.Options = options;
            % 3- Save and done
            myFinnee = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')

            %% SUB FUNCTIONS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% CHECKVARARGIN
            function options = checkVarargin_single(varargin)
                % CHECKVARARGIN is used to check the optional paramters
                % and create the options parameter.

                % 1- Defaults optional parameters
                options.FileIn               = '';
                options.FolderOut            = '';
                options.FileID               = '';
                options.Path2Fin             = '';
                options.Overwrite            = false;
                options.FileFormat           = 'mzML';
                options.XLim                 = [0 inf];
                options.YLim                 = [0 inf];
                options.RemSpks              = false;
                options.SpkSz                = 0;
                options.Precisions.Time      = '%0.2f';
                options.Precisions.Intensity = '%0.0f';
                options.Precisions.MZ        = '%0.4f';
                options.Precisions.Other     = '%0.0f';
                options.ParallelMe           = false;

                if nargin == 1 &  isstruct(varargin{1})
                    OPT = varargin{1};
                    if isfield(OPT, 'FileIn'), options.FileIn = OPT.FileIn; end
                    if isfield(OPT, 'FolderOut'), options.FolderOut = OPT.FolderOut; end
                    if isfield(OPT, 'FileID'), options.FileID = OPT.FileID; end
                    if isfield(OPT, 'Overwrite'), options.Overwrite = OPT.Overwrite; end
                    if isfield(OPT, 'XLim'), options.XLim = OPT.XLim; end
                    if isfield(OPT, 'YLim'), options.YLim = OPT.YLim; end
                    if isfield(OPT, 'RemSpks'), options.RemSpks = OPT.RemSpks; end
                    if isfield(OPT, 'SpkSz'), options.SpkSz = OPT.SpkSz; end
                    if isfield(OPT, 'UseParallel'), options.ParallelMe = true; end
                    options.ParallelMe = true;

                else

                    % 2- Decipher varargin and update options when relevamt
                    input = @(x) find(strcmpi(varargin,x),1);

                    tgtIx = input('fileIn');
                    if ~isempty(tgtIx)
                        options.FileIn = varargin{tgtIx +1};
                    end

                    tgtIx = input('folderOut');
                    if ~isempty(tgtIx)
                        options.FolderOut = varargin{tgtIx +1};
                    end

                    tgtIx = input('fileID');
                    if ~isempty(tgtIx)
                        options.FileID = varargin{tgtIx +1};
                    end

                    tgtIx = input('overwrite');
                    if ~isempty(tgtIx)
                        options.Overwrite = true;
                    end

                    tgtIx = input('tLim');
                    if ~isempty(tgtIx)
                        tLim         = varargin{tgtIx +1};
                        options.XLim = [min(tLim) max(tLim)];
                    end

                    tgtIx = input('mzLim');
                    if ~isempty(tgtIx)
                        mzLim        = varargin{tgtIx +1};
                        options.YLim = [min(mzLim) max(mzLim)];
                    end

                    tgtIx = input('spikes');
                    if ~isempty(tgtIx)
                        spks = varargin{tgtIx +1};
                        if spks == 0
                            options.RemSpks = false;
                            options.SpkSz   =  spks;
                        else
                            options.RemSpks = true;
                            options.SpkSz   =  spks;
                        end
                    end
                end
            end

            %% CHECKVARARGIN
            function options = checkVarargin_multiple(varargin)
                % CHECKVARARGIN is used to check the optional paramters
                % and create the options parameter.

                % 1- Defaults optional parameters
                options.FileIn       = '';
                options.FolderOut    = '';
                options.FileID       = '';
                options.Path2Fin     = '';
                options.Overwrite    = false;
                options.XLim         = [0 inf];
                options.YLim         = [0 inf];
                options.RemSpks      = false;
                options.SpkSz        = 3;
                options.FinneesIn     = {};
                options.Alignment    = {};
                options.MasterFinnee = {};
                options.TgtDatasets  = 2;
                options.wdw4noise       = 2;
                options.p4noise         = 2;
                options.minLen4Noise    = 100;
                options.wdw4corr        = 2;
                options.ParallelMe      = false;
                options.Spk4noise       = 3;
                options.Sig2Nois        = 3;
                options.MinMinNoise     = [100, 10];
                options.getDerivatives  = true;

                options.Precisions.Time = '%0.2f';
                options.Precisions.Intensity = '%0.0f';
                options.Precisions.MZ = '%0.4f';
                options.Precisions.Other = '%0.0f';

                % 2- Decipher varargin and update options when relevamt
                input = @(x) find(strcmpi(varargin,x),1);

                tgtIx = input('TgtDatasets');
                if ~isempty(tgtIx)
                    options.TgtDatasets = varargin{tgtIx +1};
                end

                tgtIx = input('MasterFinnee');
                if ~isempty(tgtIx)
                    options.MasterFinnee = varargin{tgtIx +1};
                end

                tgtIx = input('Alignment');
                if ~isempty(tgtIx)
                    options.Alignment = varargin{tgtIx +1};
                end

                tgtIx = input('FinneesIn');
                if ~isempty(tgtIx)
                    options.FinneesIn = varargin{tgtIx +1};
                end

                tgtIx = input('fileIn');
                if ~isempty(tgtIx)
                    options.FileIn = varargin{tgtIx +1};
                end

                tgtIx = input('folderOut');
                if ~isempty(tgtIx)
                    options.FolderOut = varargin{tgtIx +1};
                end

                tgtIx = input('fileID');
                if ~isempty(tgtIx)
                    options.FileID = varargin{tgtIx +1};
                end

                tgtIx = input('overwrite');
                if ~isempty(tgtIx)
                    options.Overwrite = true;
                end

                tgtIx = input('tLim');
                if ~isempty(tgtIx)
                    tLim         = varargin{tgtIx +1};
                    options.XLim = [min(tLim) max(tLim)];
                end

                tgtIx = input('mzLim');
                if ~isempty(tgtIx)
                    mzLim        = varargin{tgtIx +1};
                    options.YLim = [min(mzLim) max(mzLim)];
                end

                tgtIx = input('spikes');
                if ~isempty(tgtIx)
                    spks = varargin{tgtIx +1};
                    if spks == 0
                        options.RemSpks = false;
                        options.SpkSz   =  spks;
                    else
                        options.RemSpks = true;
                        options.SpkSz   =  spks;
                    end
                end

            end


            %% mk_mzML2cFinnee
            function obj = mk_mzML2cFinnee(obj, options)


                % 1- Initialisations
                % Open and check the mzML file
                tStart = tic;
                nonEndingErrors = {};

                fidRead = fopen(obj.FileIn, 'r'); % original mzML file
                if fidRead == -1
                    error('myApp:argChk', ...
                        'Error while opening the files')
                end

                metadata = getNodes(fidRead, 'mzML', 'spectrum');
                allfield = unpackAttributes(metadata, {});

                %Check if MS1 spectra are present
                if any(sum(strcmp(allfield, 'cvParam/accession') ...
                        | strcmp(allfield, 'MS:1000579'), 2)== 2)

                    %Check if MSn spectra are als0 present
                    if any(sum(strcmp(allfield, 'cvParam/accession') ...
                            | strcmp(allfield, 'MS:1000580'), 2)== 2)
                        %TODO: Warning for MSn

                    end
                    save(fullfile(obj.Path2Fin, 'AquisitionData.mat'), 'metadata')
                else

                    error("Finnee has been designed for MS1 spectra")
                end

                % Find the date of data acquision
                if any(strcmp(allfield(:, 1), 'run/startTimeStamp'))
                    dt = allfield{strcmp(allfield(:, 1), 'run/startTimeStamp'), 2};
                    dt = strrep(strrep(dt, 'T', ' '), 'Z', '');
                    I2C = strfind(dt, '+');
                    if ~isempty(I2C), dt = dt(1:I2C(1)-1); end
                    %TODO: include time zone?
                    obj.DateOfDataAcq = datetime(dt);
                end

                Profiles = {};
                myDatasets = {};

                ScanNbr = str2double(allfield{strcmp(allfield(:,1), 'spectrumList/count'), 2});
                if ~options.ParallelMe, h = waitbar(0,'processing scans'); end

                while 1

                    % GET AND DECIPHER EACH SPECTRA
                    mySpectra = getNodes(fidRead, 'spectrum', 'spectrum');
                    if isempty(fieldnames(mySpectra))
                        break
                    end

                    [myMS, myLabels, code] = mkSpectra(mySpectra);
                    if ~options.ParallelMe, waitbar(myLabels.index/ScanNbr); end

                    if  myLabels.ScanTime.value < obj.Options.XLim(1) | myLabels.ScanTime.value > obj.Options.XLim(2)
                        continue;
                    end

                    % Filter spikes if asked
                    if ~isempty(myMS) && obj.Options.RemSpks
                        spkSz = obj.Options.SpkSz;
                        myMS    = spikesRemoval(myMS, spkSz );
                    end

                    % Set to zeros the mz outside the limits
                    myMS(myMS(:,1) < options.YLim(1) | myMS(:,1) > options.YLim(2), 2) = 0;

                    % reduced trailing zero in excess
                    if ~isempty(myMS)
                        provMat      = [myMS(2:end, 2); 0];
                        provMat(:,2) = myMS(:, 2);
                        provMat(:,3) = [0; myMS(1:end-1, 2)];
                        myMS         = myMS(sum(provMat, 2) > 0, :);
                    end

                    % INFO: code will be used when different type of
                    % spectra are recorded within the original mzML file.

                    if code(1) == '1' & code(7) == '1'

                        M = ChrMoment(myMS, 2);
                        % AddVector: Scan Index | scan time | injection
                        % time | Base peak ion | sum Ion | min Intensity |
                        % non zeros values | Scan centroid | size data
                        % | centroid MS scan
                        addVector = [myLabels.index, myLabels.ScanTime.value, ...
                            myLabels.InjectionTime.value, ...
                            max(myMS(:,2)), sum(myMS(:,2)), ...
                            min(myMS(myMS(:, 2) ~= 0, 2)), ...
                            nnz(myMS(:,2)), size(myMS), M(2)];

                    elseif code(2) == '1'
                        % TODO: MSn Spectra
                        %error("MSn spectra not recognised yet")
                        %                         if ~options.ParallelMe
                        %                             sprintf('Scan %i discarded (MSn Spectra not processed by Finnee)', myLabels.index)
                        %                         end
                        nonEndingErrors{end+1} = 'MSn spectra';
                        continue;

                    elseif code(6) == '1'
                        % TODO: MSn Spectra
                        if ~options.ParallelMe
                            sprintf('Scan %i discarded (centroid spectrum detected - Finnee works with profile spectrum)', myLabels.index)
                        end
                        nonEndingErrors{end+1} = 'Centroid spectra';
                        continue;

                    else
                        error("error in the mzML file");

                    end

                    if isempty(myDatasets)
                        % IF FIRST SPECTRA IN DATASET

                        mzInterval{1} = [inf 0];
                        if min(myMS(:,1)) < mzInterval{1}(1), mzInterval{1}(1) = min(myMS(:,1)); end
                        if max(myMS(:,1)) > mzInterval{1}(2), mzInterval{1}(2) = max(myMS(:,1)); end

                        cD = 1;
                        myDatasets.Name{1, 1} = 'Dataset1';
                        myDatasets.CreationCode{1, 1} = code;
                        myDatasets.IsMS1{1, 1} = myLabels.IsMS1;
                        myDatasets.isProfileScan{1, 1} = myLabels.isProfileScan;
                        myDatasets.HasMasterMZ(1, 1) = false;
                        myDatasets.isBaseDriftCorr{1, 1} = false;
                        myDatasets.HasNoiseRem{1, 1} = false;
                        myDatasets.IsDaugherOff{1, 1} = '';
                        myDatasets.Labels{1, 1}.AxisX.Label = 'Time';
                        myDatasets.Labels{1, 1}.AxisX.Units = myLabels.ScanTime.unit;
                        myDatasets.Labels{1, 1}.AxisX.fo = options.Precisions.Time;
                        myDatasets.Labels{1, 1}.AxisY.Label = myLabels.Axis{1, 1}.name;
                        myDatasets.Labels{1, 1}.AxisY.Units = myLabels.Axis{1, 1}.name;
                        myDatasets.Labels{1, 1}.AxisY.fo = options.Precisions.MZ;
                        myDatasets.Labels{1, 1}.AxisZ.Label = myLabels.Axis{1, 2}.name;
                        myDatasets.Labels{1, 1}.AxisZ.Units = myLabels.Axis{1, 2}.unit;
                        myDatasets.Labels{1, 1}.AxisZ.fo = options.Precisions.Intensity;
                        myDatasets.Options4Creations{1, 1} = options;
                        myDatasets.PrimaryActions{1, 1} = 'Creation';
                        myDatasets.SecondaryActions{1, 1} = {};


                        mkdir(fullfile(obj.Path2Fin, myDatasets.Name{1, 1}));
                        mkdir(fullfile(obj.Path2Fin, myDatasets.Name{1, 1}, 'Scans'));

                        %TIP, BPP and more
                        Profiles{cD} = addVector;

                        %Record Spectra
                        fileName = fullfile(obj.Path2Fin, myDatasets.Name{1, 1}, 'Scans', ['Scan#', num2str(myLabels.index), '.dat']);
                        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                        fwrite(fidWriteDat, myMS(:), "double");
                        fclose(fidWriteDat);

                    else

                        if any(strcmp(code, myDatasets.CreationCode(:)))
                            % NOTE: Chek if the code correspond to the ones
                            % of an existing dataset

                            IcD = find(strcmp(code, myDatasets.CreationCode));
                            Profiles{IcD} = [Profiles{IcD}; addVector];
                            if min(myMS(:,1)) < mzInterval{IcD}(1), mzInterval{IcD}(1) = min(myMS(:,1)); end
                            if max(myMS(:,1)) > mzInterval{IcD}(2), mzInterval{IcD}(2) = max(myMS(:,1)); end

                            %Record Spectra
                            fileName = fullfile(obj.Path2Fin, myDatasets.Name{IcD, 1}, 'Scans', ['Scan#', num2str(myLabels.index), '.dat']);
                            [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                            fwrite(fidWriteDat, myMS(:), "double");
                            fclose(fidWriteDat);

                        else
                            % TODO: different type of spectra in an uniue
                            % dataset

                            error("Only file with single type of spectra is recognised at the moment")
                        end
                    end
                end
                try close(h); catch, end
                myDatasets = struct2table(myDatasets);
                obj.Datasets = myDatasets;

                for ii = 1:height(myDatasets)
                    if myDatasets.CreationCode{ii}(1) == '1'

                        clear infoDataset
                        infoDataset.dateofcreation = datetime("now");
                        infoDataset.AdditionalInformation.nonEndingErrors = nonEndingErrors;
                        infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
                        infoDataset.Label = myDatasets.Labels{ii};
                        infoDataset.mzInterval =  mzInterval{ii};
                        infoDataset.minNoise =  nan;
                        infoDataset.Options = obj.Options;
                        infoDataset.Profiles.location = size(Profiles{ii});
                        infoDataset.Profiles.size = size(Profiles{ii});
                        infoDataset.Profiles.column{1}.name = 'Scan_Index';
                        infoDataset.Profiles.column{1}.labelUnits = '#';
                        infoDataset.Profiles.column{1}.units = '';
                        infoDataset.Profiles.column{1}.formatUnits = '%0.0f';
                        infoDataset.Profiles.column{1}.ShowMe = false;
                        infoDataset.Profiles.column{2}.name = 'Scan_Time';
                        infoDataset.Profiles.column{2}.labelUnits = myDatasets.Labels{ii}.AxisX.Label;
                        infoDataset.Profiles.column{2}.units = myDatasets.Labels{ii}.AxisX.Units;
                        infoDataset.Profiles.column{2}.formatUnits = myDatasets.Labels{ii}.AxisX.fo;
                        infoDataset.Profiles.column{2}.ShowMe = false;
                        infoDataset.Profiles.column{3}.name = 'Injection_Time';
                        infoDataset.Profiles.column{3}.labelUnits = 'Time';
                        infoDataset.Profiles.column{3}.units = myLabels.InjectionTime.unit;
                        infoDataset.Profiles.column{3}.formatUnits = myDatasets.Labels{ii}.AxisX.fo;
                        infoDataset.Profiles.column{3}.ShowMe = true;
                        infoDataset.Profiles.column{4}.name = 'Base_Peak_Profile';
                        infoDataset.Profiles.column{4}.labelUnits = myDatasets.Labels{ii}.AxisZ.Label;
                        infoDataset.Profiles.column{4}.units = myDatasets.Labels{ii}.AxisZ.Units;
                        infoDataset.Profiles.column{4}.formatUnits = myDatasets.Labels{ii}.AxisZ.fo;
                        infoDataset.Profiles.column{4}.ShowMe = true;
                        infoDataset.Profiles.column{5}.name = 'Total_Ions_Profile';
                        infoDataset.Profiles.column{5}.labelUnits = myDatasets.Labels{ii}.AxisZ.Label;
                        infoDataset.Profiles.column{5}.units = myDatasets.Labels{ii}.AxisZ.Units;
                        infoDataset.Profiles.column{5}.formatUnits = myDatasets.Labels{ii}.AxisZ.fo;
                        infoDataset.Profiles.column{5}.ShowMe = true;
                        infoDataset.Profiles.column{6}.name = 'Lowest_Intensity_Profile';
                        infoDataset.Profiles.column{6}.labelUnits = myDatasets.Labels{ii}.AxisZ.Label;
                        infoDataset.Profiles.column{6}.units = myDatasets.Labels{ii}.AxisZ.Units;
                        infoDataset.Profiles.column{6}.formatUnits = myDatasets.Labels{ii}.AxisZ.fo;
                        infoDataset.Profiles.column{6}.ShowMe = false;
                        infoDataset.Profiles.column{7}.name = 'Non_Null_Elements_Profile';
                        infoDataset.Profiles.column{7}.labelUnits = '';
                        infoDataset.Profiles.column{7}.units = '';
                        infoDataset.Profiles.column{7}.formatUnits = '%0.0f';
                        infoDataset.Profiles.column{7}.ShowMe = false;
                        infoDataset.Profiles.column{8}.name = 'Length_MS_Scan_Profile';
                        infoDataset.Profiles.column{8}.labelUnits = '';
                        infoDataset.Profiles.column{8}.units = '';
                        infoDataset.Profiles.column{8}.formatUnits = '%0.0f';
                        infoDataset.Profiles.column{8}.ShowMe = false;
                        infoDataset.Profiles.column{9}.name = 'Column_MS_Scan_Profile';
                        infoDataset.Profiles.column{9}.labelUnits = '';
                        infoDataset.Profiles.column{9}.units = '';
                        infoDataset.Profiles.column{9}.formatUnits = '%0.0f';
                        infoDataset.Profiles.column{9}.ShowMe = false;
                        infoDataset.Profiles.column{10}.name = 'MS_centroid_profile';
                        infoDataset.Profiles.column{10}.labelUnits = myDatasets.Labels{ii}.AxisY.fo;
                        infoDataset.Profiles.column{10}.units = myDatasets.Labels{ii}.AxisY.Units;
                        infoDataset.Profiles.column{10}.formatUnits = myDatasets.Labels{ii}.AxisY.fo;
                        infoDataset.Profiles.column{10}.ShowMe = true;
                        infoDataset.Spectra = {};
                        save(fullfile(obj.Path2Fin, myDatasets.Name{ii}, 'infoDataset.mat'), 'infoDataset');

                        %Record Spectra
                        fileName = fullfile(obj.Path2Fin, myDatasets.Name{ii}, 'Profiles.dat');
                        [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
                        fwrite(fidWriteDat, Profiles{ii}(:), "double");
                        fclose(fidWriteDat);

                    end
                end
                fclose(fidRead);

            end

            %% SUBFUNCTION TO mk_mzML2cFinnee
            function [Spectra, Label, code] = mkSpectra(SpectraData)
                % Decifer first level

                Label.index = str2double(SpectraData.spectrum.Attributes.index);
                Label.InjectionTime.value = nan;
                Label.InjectionTime.unit = '';

                % NOTE: code is used to recognised the type of spectra.
                % If MS1 spectrum      code = 10XXXXXXXXXX
                % If MSn spectrum      code = 01nXXXXXXXXX (n MS level)
                % If negative polarity code = XXX10XXXXXXX
                % If positive polarity code = XXX01XXXXXXX
                % If unknow polarity   code = XXX00XXXXXXX
                % If centroid scan     code = XXXXX10XXXXX
                % If profile scan      code = XXXXX01XXXXX
                % If unknown mode      code = XXXXX00XXXXX
                code = '000000000000';
                Spectra = [];

                for ii = 1:numel(SpectraData.spectrum.subElements)
                    if isfield(SpectraData.spectrum.subElements{ii}, 'cvParam')
                        cvParam = SpectraData.spectrum.subElements{ii}.cvParam.Attributes;

                        switch cvParam.accession

                            case 'MS:1000579'
                                % MS1 spectrum
                                code(1) = '1';
                                Label.IsMS1 = true;

                            case 'MS:1000580'
                                % MSn spectrum
                                code(2) = '1';
                                Label.IsMS1 = false;

                            case 'MS:1000511'
                                % MS level (1 for MS1)
                                level = cvParam.value(1);
                                if ~isfinite(level), error(""); end
                                code(3) = num2str(level);
                                Label.level = level;

                            case 'MS:1000129'
                                % negative scan
                                code(4) = '1';
                                Label.isPositivePolarity = false;

                            case 'MS:1000130'
                                % Positive scan
                                code(5) = '1';
                                Label.isPositivePolarity = true;

                            case 'MS:1000127'
                                % centroid spectrum
                                code(6) = '1';
                                Label.isProfileScan = false;

                            case 'MS:1000128'
                                % profile spectrum
                                code(7) = '1';
                                Label.isProfileScan = true;

                        end

                    elseif isfield(SpectraData.spectrum.subElements{ii}, 'scanList')

                        for jj = 1:numel(SpectraData.spectrum.subElements{ii}.scanList.subElements)

                            if isfield(SpectraData.spectrum.subElements{ii}.scanList.subElements{jj}, 'scan')
                                scan = SpectraData.spectrum.subElements{ii}.scanList.subElements{jj}.scan;

                                for kk = 1:numel(scan.subElements)

                                    if isfield(scan.subElements{kk}, 'cvParam')

                                        switch scan.subElements{kk}.cvParam.Attributes.accession

                                            case 'MS:1000016'
                                                Label.ScanTime.value = str2double(scan.subElements{kk}.cvParam.Attributes.value);

                                                switch scan.subElements{kk}.cvParam.Attributes.unitAccession
                                                    case 'UO:0000031'
                                                        Label.ScanTime.unit = 'min';

                                                    case 'UO:0000010'
                                                        Label.ScanTime.unit = 's';

                                                    otherwise
                                                        Label.ScanTime.unit = '';

                                                end

                                            case 'MS:1000927'
                                                Label.InjectionTime.value = str2double(scan.subElements{kk}.cvParam.Attributes.value);

                                                switch scan.subElements{kk}.cvParam.Attributes.unitAccession
                                                    case 'UO:0000028'
                                                        Label.InjectionTime.unit = 'ms';

                                                    case 'UO:0000010'
                                                        Label.InjectionTime.unit = 's';

                                                    otherwise
                                                        Label.InjectionTime.unit = '';

                                                end
                                        end
                                    end
                                end
                            end
                        end

                    elseif isfield(SpectraData.spectrum.subElements{ii}, 'binaryDataArrayList')

                        for jj = 1:numel(SpectraData.spectrum.subElements{ii}.binaryDataArrayList.subElements)
                            if isfield(SpectraData.spectrum.subElements{ii}.binaryDataArrayList.subElements{jj}, 'binaryDataArray')
                                bda = SpectraData.spectrum.subElements{ii}.binaryDataArrayList.subElements{jj}.binaryDataArray;

                                for kk = 1:numel(bda.subElements)
                                    if isfield(bda.subElements{kk}, 'cvParam')

                                        switch bda.subElements{kk}.cvParam.Attributes.accession

                                            case 'MS:1000523'
                                                dataFormat= 'MS:1000523';

                                            case 'MS:1000574'
                                                compression = 'MS:1000574';

                                            case 'MS:1000514'
                                                Label.Axis{jj}.name = 'm/z';

                                                switch bda.subElements{kk}.cvParam.Attributes.unitAccession
                                                    case 'MS:1000040'
                                                        Label.Axis{jj}.unit = '';

                                                end

                                            case 'MS:1000515'
                                                Label.Axis{jj}.name = 'Intensity';

                                                switch bda.subElements{kk}.cvParam.Attributes.unitAccession
                                                    case 'MS:1000131'
                                                        Label.Axis{jj}.unit = 'arb. unit';

                                                end
                                        end
                                    elseif isfield(bda.subElements{kk}, 'binary')

                                        input = bda.subElements{kk}.binary.Content;
                                        if ~isempty(input)

                                            switch dataFormat
                                                case 'MS:1000523'
                                                    output = base64decode(input);

                                                otherwise
                                                    error('precision not recognized')

                                            end

                                            switch compression
                                                case 'MS:1000574'
                                                    output = zlibdecode(output);

                                                otherwise
                                                    error('compression not recognized')

                                            end

                                            Spectra(:, end+1) = typecast(uint8(output),'double');
                                        end
                                    end
                                end
                            end
                        end
                    else

                        error("")
                    end

                end

                %CHECK CODE FOR INCOHERENCE
                if code(1) == '1' && code(2) == '1', error('incoherence on the mzML file'); end
                if code(1) == '0' && code(2) == '0', error('incoherence on the mzML file'); end
                if code(4) == '1' && code(5) == '1', error('incoherence on the mzML file'); end
                if code(6) == '1' && code(7) == '1', error('incoherence on the mzML file'); end
                if code(6) == '0' && code(7) == '0', error('incoherence on the mzML file'); end
            end
        end

        function myProfile = getProfile(obj, dts, Target)

            infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'infoDataset.mat'));
            infoDataset = infoDataset.infoDataset;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Profiles.dat'), 'rb');
            TimeSpectra = fread(fidRead, [infoDataset.Profiles.size], "double");
            fclose(fidRead);
            TimeSpectra = array2table(TimeSpectra);
            ColumnNames = {}; ColumnUnits = {}; LabelUnits = {}; FormatUnits = {}; IdShow = [];
            for ii = 1:numel(infoDataset.Profiles.column)
                ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
                LabelUnits{ii} = infoDataset.Profiles.column{ii}.labelUnits;
                FormatUnits{ii} = infoDataset.Profiles.column{ii}.formatUnits;
                IdShow(ii) = infoDataset.Profiles.column{ii}.ShowMe;
            end
            TimeSpectra.Properties.VariableNames = ColumnNames;
            TimeSpectra.Properties.VariableUnits = ColumnUnits;

            if Target == "?"
                fprintf('%s \n', ColumnNames{IdShow == 1});
                return;
            end

            myFinnee = obj;
            AdiPrm = struct();
            AdiPrm.Finnee.object = myFinnee;
            AdiPrm.Finnee.dataset = dts;
            AdiPrm.Finnee.Target = Target;

            IdTarget = find(strcmp(ColumnNames, Target));
            IdST = find(strcmp(ColumnNames, 'Scan_Time'));
            data2write = [TimeSpectra.Scan_Time, TimeSpectra.(Target)];
            infoTrc.Title = Target;
            infoTrc.FT = '';
            infoTrc.TT = 'EMP';
            infoTrc.AxisX.Label = LabelUnits{IdST};
            infoTrc.AxisX.Unit = ColumnUnits{IdST};
            infoTrc.AxisX.fo = FormatUnits{IdST};
            infoTrc.AxisY.Label = LabelUnits{IdTarget};
            infoTrc.AxisY.Unit = ColumnUnits{IdTarget};
            infoTrc.AxisY.fo = FormatUnits{IdTarget};
            infoTrc.Data = [];
            infoTrc.AdiPrm  = AdiPrm;
            myProfile = Profile(infoTrc, data2write);
        end

        function mySpectrum = getSpectrum(obj, dts, Target)

            if dts > size(obj.Datasets, 1)
                error('The target dataset does not exist yet')
            end

            if ~obj.Datasets.HasMasterMZ(dts)
                error('Align all MS scan to a master MZ axis before using this function')

            end

            infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'infoDataset.mat'));
            infoDataset = infoDataset.infoDataset;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Spectra.dat'), 'rb');
            Spectra = fread(fidRead, [infoDataset.Spectra.size], "double");
            fclose(fidRead);
            Spectra = array2table(Spectra);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Spectra.column)
                ColumnNames{ii} = infoDataset.Spectra.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Spectra.column{ii}.units;
                LabelUnits{ii} = infoDataset.Spectra.column{ii}.labelUnits;
                FormatUnits{ii} = infoDataset.Spectra.column{ii}.formatUnits;
                IdShow(ii) = infoDataset.Spectra.column{ii}.ShowMe;
            end
            Spectra.Properties.VariableNames = ColumnNames;
            Spectra.Properties.VariableUnits = ColumnUnits;
            if Target == "?"
                fprintf('%s \n', ColumnNames{IdShow == 1});
                return;
            end
            if Target == "?all"
                fprintf('%s \n', ColumnNames{:});
                return;
            end

            IdTarget = find(strcmp(ColumnNames, Target));
            IdST = find(strcmp(ColumnNames, 'm/z axis'));
            data2write = [Spectra.('m/z axis'), Spectra.(Target)];

            myFinnee = obj;
            AdiPrm = struct();
            AdiPrm.Finnee.object = myFinnee;
            AdiPrm.Finnee.dataset = dts;
            AdiPrm.Finnee.Target = Target;

            infoTrc.Title = Target;
            infoTrc.FT = '';
            infoTrc.TT = 'EMP';
            infoTrc.AxisX.Label = LabelUnits{IdST};
            infoTrc.AxisX.Unit = ColumnUnits{IdST};
            infoTrc.AxisX.fo = FormatUnits{IdST};
            infoTrc.AxisY.Label = LabelUnits{IdTarget};
            infoTrc.AxisY.Unit = ColumnUnits{IdTarget};
            infoTrc.AxisY.fo = FormatUnits{IdTarget};
            infoTrc.Data = [];
            infoTrc.AdiPrm  = AdiPrm;
            mySpectrum = Spectrum(infoTrc, data2write);
        end

        function mySpectrum = mkSpectrum(obj, dts, TmLim, mzInt)

            narginchk(3, 10)

            if dts > size(obj.Datasets, 1)
                error('The target dataset does not exist yet')
            end

            if ~obj.Datasets.HasMasterMZ(dts)
                error('Align all MS scan to a master MZ axis before using this function')

            end

            if nargin == 3
                mzInt = [ 0 inf];

            end

            if numel (TmLim) == 1
                TmLim = [TmLim TmLim];

            elseif numel (TmLim) > 2
                error('');
            end

            if numel (mzInt) > 2
                error('');
            end

            infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'infoDataset.mat'));
            infoDataset = infoDataset.infoDataset;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Profiles.dat'), 'rb');
            Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
            fclose(fidRead);
            Profiles = array2table(Profiles);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Profiles.column)
                ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
                LabelUnits{ii} = infoDataset.Profiles.column{ii}.labelUnits;
                FormatUnits_1{ii} = infoDataset.Profiles.column{ii}.formatUnits;
            end
            Profiles.Properties.VariableNames = ColumnNames;
            Profiles.Properties.VariableUnits = ColumnUnits;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Spectra.dat'), 'rb');
            Spectra = fread(fidRead, [infoDataset.Spectra.size], "double");
            fclose(fidRead);
            Spectra = array2table(Spectra);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Spectra.column)
                ColumnNames{ii} = infoDataset.Spectra.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Spectra.column{ii}.units;
                LabelUnits{ii} = infoDataset.Spectra.column{ii}.labelUnits;
                FormatUnits_2{ii} = infoDataset.Spectra.column{ii}.formatUnits;
                IdShow(ii) = infoDataset.Spectra.column{ii}.ShowMe;
            end
            Spectra.Properties.VariableNames = ColumnNames;
            Spectra.Properties.VariableUnits = ColumnUnits;

            IdS = find(Profiles.Scan_Time <= TmLim(1), 1, 'last');
            if isempty(IdS), IdS = 1; end

            IdE = find(Profiles.Scan_Time >= TmLim(2), 1, 'first');
            if isempty(IdS), IdS = numel(Profiles.ScanTime); end

            fileName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)], 'MasterMZAxis.dat');
            [fidReadDat, errmsg]  = fopen(fileName, 'rb');
            XY = fread(fidReadDat, inf, "double");
            fclose(fidReadDat);
            XY(:,2) = 0;

            for ii = IdS:IdE
                ScanName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)],...
                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
                fidRead = fopen(ScanName, 'rb');
                myMS = fread(fidRead, ...
                    [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                    "double");
                fclose(fidRead);
                [~, IdX] = intersect(XY(:,1), myMS(:,1));
                XY(IdX, 2) = XY(IdX, 2) + myMS(:,2);

            end

            XY(XY(:,1) < mzInt(1) | XY(:,1) > mzInt(2), :) = [];

            % reduced trailing zero in excess
            if ~isempty(XY)
                provMat      = [XY(2:end, 2); 0];
                provMat(:,2) = XY(:, 2);
                provMat(:,3) = [0; XY(1:end-1, 2)];
                XY         = XY(sum(provMat, 2) > 0, :);
            end

            myFinnee = obj;
            AdiPrm = struct();
            AdiPrm.Finnee.object = myFinnee;
            AdiPrm.Finnee.dataset = dts;
            AdiPrm.Finnee.Target = TmLim;

            Text = {sprintf('%s, dataset %i', obj.FileID, dts)...
                sprintf(['Averaged Spectrum between ', FormatUnits_1{2},' and ', FormatUnits_1{2},' ', ColumnUnits{2}], Profiles.Scan_Time([IdS, IdE]))};

            infoTrc.Title = Text;
            infoTrc.FT = '';
            infoTrc.TT = 'EMP';
            infoTrc.AxisX.Label = infoDataset.Label.AxisY.Label;
            infoTrc.AxisX.Unit = '';
            infoTrc.AxisX.fo = infoDataset.Label.AxisY.fo;
            infoTrc.AxisY.Label = infoDataset.Label.AxisZ.Label;
            infoTrc.AxisY.Unit = infoDataset.Label.AxisZ.Units;
            infoTrc.AxisY.fo = infoDataset.Label.AxisZ.fo;
            infoTrc.AdiPrm  = AdiPrm;
            mySpectrum = Spectrum(infoTrc, XY);


        end

        function myProfile = mkProfile(obj, dts, MzLim, TmLim)

            narginchk(3, 10)

            if dts > size(obj.Datasets, 1)
                error('The target dataset does not exist yet')
            end

            if ~obj.Datasets.HasMasterMZ(dts)
                error('Align all MS scan to a master MZ axis before using this function')

            end

            if nargin == 3
                TmLim = [ 0 inf];

            end

            if numel (MzLim) == 1
                MzLim = [MzLim MzLim];

            elseif numel (TmLim) > 2
                error('');
            end

            if numel (MzLim) > 2
                error('');
            end

            infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'infoDataset.mat'));
            infoDataset = infoDataset.infoDataset;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Profiles.dat'), 'rb');
            Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
            fclose(fidRead);
            Profiles = array2table(Profiles);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Profiles.column)
                ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
                LabelUnits{ii} = infoDataset.Profiles.column{ii}.labelUnits;
                FormatUnits{ii} = infoDataset.Profiles.column{ii}.formatUnits;
            end
            Profiles.Properties.VariableNames = ColumnNames;
            Profiles.Properties.VariableUnits = ColumnUnits;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Spectra.dat'), 'rb');
            Spectra = fread(fidRead, [infoDataset.Spectra.size], "double");
            fclose(fidRead);
            Spectra = array2table(Spectra);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Spectra.column)
                ColumnNames{ii} = infoDataset.Spectra.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Spectra.column{ii}.units;
                LabelUnits{ii} = infoDataset.Spectra.column{ii}.labelUnits;
                FormatUnits{ii} = infoDataset.Spectra.column{ii}.formatUnits;
                IdShow(ii) = infoDataset.Spectra.column{ii}.ShowMe;
            end
            Spectra.Properties.VariableNames = ColumnNames;
            Spectra.Properties.VariableUnits = ColumnUnits;

            ItS = find(Profiles.Scan_Time <= TmLim(1), 1, 'last');
            if isempty(ItS), ItS = 1; end
            ItE = find(Profiles.Scan_Time >= TmLim(2), 1, 'first');
            if isempty(ItE), ItE = numel(Profiles.Scan_Time); end

            XY2 = zeros(ItE-ItS+1, 3);

            ImS = findCloser(MzLim(1), Spectra.("m/z axis"));
            if isempty(ImS), ImS = 1; end
            ImE = findCloser(MzLim(2), Spectra.("m/z axis"));
            if isempty(ImE), ImE = numel(Spectra.("m/z axis")); end
            for ii = ItS:ItE
                ScanName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)],...
                    'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
                fidRead = fopen(ScanName, 'rb');
                myMS = fread(fidRead, ...
                    [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                    "double");
                fclose(fidRead);

                [~, Id2Keep] = intersect(myMS(:, 1), Spectra.("m/z axis")(ImS:ImE));

                if ~isempty(Id2Keep)

                    XY2(ii-ItS+1, 1) = Profiles.Scan_Time(ii);
                    XY2(ii-ItS+1, 2) = sum(myMS(Id2Keep,2));
                    XY2(ii-ItS+1, 3) = max(myMS(Id2Keep,2));

                else
                    XY2(ii-ItS+1, 1) = Profiles.Scan_Time(ii);

                end
            end

            myFinnee = obj;
            AdiPrm = struct();
            AdiPrm.Finnee.object = myFinnee;
            AdiPrm.Finnee.dataset = dts;
            AdiPrm.Finnee.Target = MzLim;

            Text = {sprintf('%s, dataset %i', obj.FileID, dts)...
                sprintf(['Extracted Ions profile m/z: [', obj.Options.Precisions.MZ,':', obj.Options.Precisions.MZ,'] '], MzLim)};

            infoTrc.Title = Text;
            infoTrc.FT = '';
            infoTrc.TT = 'EMP';
            infoTrc.AxisX.Label = infoDataset.Label.AxisX.Label;
            infoTrc.AxisX.Unit = infoDataset.Label.AxisX.Units;
            infoTrc.AxisX.fo = infoDataset.Label.AxisX.fo;
            infoTrc.AxisY.Label = infoDataset.Label.AxisZ.Label;
            infoTrc.AxisY.Unit = infoDataset.Label.AxisZ.Units;
            infoTrc.AxisY.fo = infoDataset.Label.AxisZ.fo;
            infoTrc.Data = [];
            infoTrc.AdiPrm  = AdiPrm;
            myProfile = Profile(infoTrc, XY2);
        end

        function myROI = mkROI(obj, dts, TmInt, mzInt)

            narginchk(4, 10)

            if dts > size(obj.Datasets, 1)
                error('The target dataset does not exist yet')
            end

            if ~obj.Datasets.HasMasterMZ(dts)
                error('Align all MS scan to a master MZ axis before using this function')

            end

            if numel (TmInt) ~= 2

                error('error')

            end

            if numel (mzInt) ~= 2

                error('error');
            end

            infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'infoDataset.mat'));
            infoDataset = infoDataset.infoDataset;

            fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{dts}, 'Profiles.dat'), 'rb');
            TimeSpectra = fread(fidRead, [infoDataset.Profiles.size], "double");
            fclose(fidRead);
            TimeSpectra = array2table(TimeSpectra);
            ColumnNames = {}; ColumnUnits = {};
            for ii = 1:numel(infoDataset.Profiles.column)
                ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
                ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
                LabelUnits{ii} = infoDataset.Profiles.column{ii}.labelUnits;
                FormatUnits{ii} = infoDataset.Profiles.column{ii}.formatUnits;
            end
            TimeSpectra.Properties.VariableNames = ColumnNames;
            TimeSpectra.Properties.VariableUnits = ColumnUnits;

            IdS = find(TimeSpectra.Scan_Time <= TmInt(1), 1, 'last');
            if isempty(IdS), IdS = 1; end

            IdE = find(TimeSpectra.Scan_Time >= TmInt(2), 1, 'first');
            if isempty(IdS), IdS = numel(TimeSpectra.ScanTime); end

            fileName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)], 'MasterMZAxis.dat');
            [fidReadDat, errmsg]  = fopen(fileName, 'rb');

            AxisMZ = fread(fidReadDat, inf, "double");
            fclose(fidReadDat);
            AxisMZ(AxisMZ < mzInt(1) | AxisMZ > mzInt(2)) = [];

            XYZ = zeros(numel(AxisMZ), IdE-IdS+1);

            for ii = IdS:IdE
                ScanName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)],...
                    'Scans', ['Scan#', num2str(TimeSpectra.Scan_Index(ii, 1)), '.dat']);
                fidRead = fopen(ScanName, 'rb');
                myMS = fread(fidRead, ...
                    [TimeSpectra.Length_MS_Scan_Profile(ii) TimeSpectra.Column_MS_Scan_Profile(ii)],...
                    "double");
                fclose(fidRead);
                [~, IdX1] = intersect(myMS(:,1), AxisMZ);
                [~, IdX2] = intersect(AxisMZ, myMS(IdX1,1));
                XYZ(IdX2, ii-IdS+1) = myMS(IdX1,2);

            end

            myROI = XYZ;
        end
    end
end
