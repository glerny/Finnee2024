function obj = addFiles(obj)

%% INTRODUCTION

[file, path] = uigetfile('*.mzML', 'MultiSelect', 'on');
Folder1 = obj.OutputFolder;
Summary = obj.Summary;
FileID = {}; StartTime = {}; EndTime = {};
FolderID = {}; FinneeType = {};
Finnee_Master = obj.MasterObject;


% parfor ii = 1:length(file)
for ii = 1:length(file)
    StartTime{ii, 1} = datetime;
    FinneeType{ii, 1} = 'original';
    [~, fileID, ~] = fileparts(file{ii});
    disp([fileID, ' ... STARTING ...']),
    FileID{ii, 1} = fileID;
    interval = [1, size(Finnee_Master.Datasets, 1)];

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

   
end

newTable = table(FileID, FolderID, FinneeType, StartTime, EndTime);
newTable.Dataset4Quali = ones(numel(EndTime), 1)*Summary.Dataset4Quali(1);
newTable.Dataset4Quant = ones(numel(EndTime), 1)*Summary.Dataset4Quant(1);
obj.Summary = [Summary; newTable];
myProject = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')



