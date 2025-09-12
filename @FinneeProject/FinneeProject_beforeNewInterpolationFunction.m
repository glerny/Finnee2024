%% DESCRIPTION
% PROFILE is the class that deals with two-dimensional representations.
% Those can either be electropherograms, chromatograms, MS spectra or
% others.
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef FinneeProject

    properties
        Title        %
        Description  %
        ListOfFiles  %
        OutputFolder % I
        Path2Project % I
        MasterObject %
        Summary
        FeaturesLists
        QuantificationLists

    end


    methods

        function obj = FinneeProject(Finnee_Master, type, interval)

            narginchk(1, 3)
            if nargin <= 1
                type = 'original';
                interval = [1, size(Finnee_Master.Datasets, 1)];

            elseif nargin <= 2
                interval = [1, size(Finnee_Master.Datasets, 1)];
            end

            switch type
                case 'original'

                    [file, path] = uigetfile('*.mzML', 'MultiSelect', 'on');
                    folderOut = uigetdir();
                    Folder1 = fullfile(folderOut, 'Finnee Files');
                    mkdir(Folder1)
                    obj.Title = '';
                    obj.Description = '';
                    obj.ListOfFiles = file;
                    obj.OutputFolder = Folder1;
                    obj.Path2Project = folderOut;
                    obj.MasterObject = Finnee_Master;
                    obj.Summary = table();
                    options.Dataset4Quali = 4;
                    options.Dataset4Quant = 2;


                    Summary = table();
                    FileID = {}; StartTime = {}; EndTime = {};
                    FolderID = {}; FinneeType = {};

                    for ii = 1:length(file)
                    % parfor ii = 1:length(file)
                        StartTime{ii, 1} = datetime;
                        FinneeType{ii, 1} = 'original';
                        [~, fileID, ~] = fileparts(file{ii});
                        disp([fileID, ' ... STARTING ...']),
                        FileID{ii, 1} = fileID;

%                         try

                            for jj =  interval(1):interval(2)

                                options = Finnee_Master.Datasets.Options4Creations{jj};
                                switch  Finnee_Master.Datasets.PrimaryActions{jj}

                                    case 'Creation'
                                        options.FileIn = fullfile(path, file{ii});
                                        options.FolderOut = Folder1;
                                        options.FileID = fileID;
                                        myFinnee = Finnee('Single', options);

                                    case 'Master_mz_axis'
                                        myFinnee = myFinnee.Interpolate2D(options);

                                    case 'Baseline_correction'
                                        myFinnee = myFinnee.BaselineCorrection(options);

                                    case 'Noise_removal'
                                        myFinnee = myFinnee.filterDataset(options);

                                    otherwise
                                        error()

                                end

                                if ~isempty(Finnee_Master.Datasets.SecondaryActions{jj})

                                    for kk = 1:size(Finnee_Master.Datasets.SecondaryActions{jj}, 1)
                                        options4SA = Finnee_Master.Datasets.SecondaryActions{jj}.Options4Creations{kk};

                                        switch Finnee_Master.Datasets.SecondaryActions{jj}.TypeOfAction{kk}
                                            case 'getROIs'
                                                if isfield(options4SA.filterROIs, 'model')
                                                    myModel = options4SA.filterROIs.model;

                                                else
                                                    myModel = [];
                                                end
                                                [myROIs, myFinnee] = myFinnee.getROIs(4, options4SA);

                                            case 'mkPeakTable'
                                              [myPeakList, myFinnee] = mkPeakList(myROIs, myFinnee, 'analysis_1', options4SA);


                                            otherwise
                                                error()

                                        end
                                    end
                                end
                            end
                            disp([fileID, ' ... DONE!']),

                            FolderID{ii, 1} = Folder1;
                            EndTime{ii, 1} = datetime;

%                         catch
%                             disp([fileID, ' ... ERROR!']),
% 
%                             FolderID{ii, 1} = Folder1;
%                             EndTime{ii, 1} = NaT;
% 
%                         end
                    end

            end

            obj.Summary = table(FileID, FolderID, FinneeType, StartTime, EndTime);
            obj.Summary.Dataset4Quali = ones(size(obj.Summary, 1), 1)*options.Dataset4Quali;
            obj.Summary.Dataset4Quant = ones(size(obj.Summary, 1), 1)*options.Dataset4Quant;
            myProject = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')

        end

    end

end

