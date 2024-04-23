function [obj, myMrgPeakList] = mkMergedPeakList(obj, Tag2Fin, options)

%% INTRODUCTION
if isempty(options)
    options.remPeakTable        = false;

    options.TimeInterval        = [0 inf];
    options.LimitResolution     = [0.8 1.3 2];
    options.alfaClustering      = 2;
    options.PearsonCluster      = 0.8;
    options.KolmogorovCluster   = 0.05;
    options.Sig2PkLim           = 2;
    options.IntThres            = 0.01;
    options.maxRs               = 10;
    options.maxPPM              = 3;
    options.maxTm               = 0.5;

    options.Merge.do            = false;
    options.Filter.do           = false;

    options.Analysis.minMZ      = 3;
    options.Analysis.minTime    = 4;
    options.Analysis.WindowMZ   = 2;
    options.Analysis.WindowTime = 3;
    options.Analysis.thr4bsl    = 0.95;
    options.Analysis.maxSkim    = 0.50;
    options.Analysis.Sig2Noise  = 3;
    options.Analysis.doBslCorr  = false;

    options.cleanMe.Thresh      = 0.01;

    options.Merge.SigMZ         = 3;
    options.Merge.SigTm         = 3;
    options.Merge.minRepl       = 3;
    options.Merge.S2N           = [3, 10, 25];
    options.Merge.minOcc4Stat   = 3;

    
end
optionsGMM = statset('Display', 'off', 'MaxIter', 1000);
MrgPeakList = table();
ListID      = [];
alfa        = options.alfaClustering;
tStart      = tic;
cFile       = 0;
IdFi        = [];

%% MERGE LISTS
for ii = 1:size(obj.Summary, 1)
    if strcmp(obj.Summary.FinneeType{ii}, Tag2Fin)
        cFile = cFile + 1;
        FilesUsed{cFile} = obj.Summary.FileID{ii};
        FilesID(cFile)   = ii;
        Dataset4Peaklist = ['Dataset', num2str(obj.Summary.Dataset4Quali(ii))];
        IdfPT = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], Dataset4Peaklist, 'myPeakList.mat');
        if exist(IdfPT, 'file') == 2 && ~ options.remPeakTable
            continue

        else
            myFile = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], Dataset4Peaklist, 'myROIs.mat');
            LMF = load(myFile);
            myROIs = LMF.myROIs;

            myPeakList = mkPeakList(myROIs, 'analysis_1', options);
            save(IdfPT, 'myPeakList')
        end
    end
    
end

if  options.Merge.do
    disp("pp")
end

if  options.Filter.do
    disp("pp")
end

for ii = 1:numel(FilesUsed)
    IdF = FilesID(ii);
    Dataset4Peaklist = ['Dataset', num2str(obj.Summary.Dataset4Quali(IdF))];
    IdfPT = fullfile(obj.Summary.FolderID{IdF}, [obj.Summary.FileID{IdF}, '.fin'], Dataset4Peaklist, 'myPeakList.mat');
    load(IdfPT, "myPeakList");
    myPeakList = sortrows(myPeakList,"chrom_Time_IA","ascend");

    myPeakList.Id2Dts = ii*ones(size(myPeakList, 1), 1);
    MrgPeakList = [MrgPeakList; myPeakList];
end


%% CREATE FEATURES FROM SUPERPOSING 3D PEAKS
Nodes = {};
countMe = 1;
while ~isempty(MrgPeakList)
    MrgPeakList = sortrows(MrgPeakList,"chrom_intensity_apex","descend");
    RS = (MrgPeakList.chrom_centroid(1) - MrgPeakList.chrom_centroid)./...
        (2*(sqrt(MrgPeakList.chrom_variance(1)) + sqrt(MrgPeakList.chrom_variance)));
    RS(:, 2) = (MrgPeakList.ms_Accu_mass_2(1) - MrgPeakList.ms_Accu_mass_2)./...
        (2*(sqrt(MrgPeakList.ms_variance(1)) + sqrt(MrgPeakList.ms_variance)));
    RS(:, 3) = sqrt(RS(:, 2).^2 + RS(:, 1).^2);
    [~, idX] = sort(RS(:, 3));
    MrgPeakList = MrgPeakList(idX, :);
    RS = RS(idX, :);

    mList = 1;
    LDst = MrgPeakList.Id2Dts(1);
    for ii = 2:numel(idX)

        if RS(ii, 3) > options.maxRs
            break

        end

        if abs((mean(MrgPeakList.ms_Accu_mass_2(mList)) - MrgPeakList.ms_Accu_mass_2(ii)))/mean(MrgPeakList.ms_Accu_mass_2(mList))*1000000 < options.maxPPM ...
                & abs(mean(MrgPeakList.chrom_centroid(mList)) - MrgPeakList.chrom_centroid(ii)) < options.maxTm
            if any(LDst == MrgPeakList.Id2Dts(ii))
                break

            else
                mList(end+1) = ii;
                LDst(end + 1) = MrgPeakList.Id2Dts(ii);

            end
        end

    end

    while 1
        stopMe = true;

        if size(mList, 2) > 1
            dList = [];
            for ii = 1:numel(mList)

                RS = (MrgPeakList.chrom_centroid(mList(ii)) - MrgPeakList.chrom_centroid)./...
                    (2*(sqrt(MrgPeakList.chrom_variance(mList(ii))) + sqrt(MrgPeakList.chrom_variance)));
                RS(:, 2) = (MrgPeakList.ms_Accu_mass_2(mList(ii)) - MrgPeakList.ms_Accu_mass_2)./...
                    (2*(sqrt(MrgPeakList.ms_variance(mList(ii))) + sqrt(MrgPeakList.ms_variance)));
                RS(:, 3) = sqrt(RS(:, 2).^2 + RS(:, 1).^2);


                IdA = find(RS(:, 3) == min(RS(RS(:, 3) ~= 0, 3)));
                if ~any(mList == IdA)
                    dList(ii) = 1;
                    stopMe = false;

                else
                    dList(ii) = 0;

                end
            end
            mList = mList(dList == 0);
        end

        if stopMe, break; end
    end

    if numel(mList) > options.Merge.minRepl
        ia = isoutlier(MrgPeakList.ms_Accu_mass_2(mList), "ThresholdFactor", 5);
        mList(ia) = [];

    end
   
    if isempty(mList), mList = 1; end

    if numel(mList) > options.Merge.minRepl
        fPeaks = MrgPeakList(mList, :);
        Nodes{countMe} = fPeaks;

        SNodes(countMe, :) = [countMe, size(fPeaks, 1), mean(fPeaks.chrom_area), ...
            mean(fPeaks.chrom_Time_IA), std(fPeaks.chrom_Time_IA),...
            mean(fPeaks.chrom_centroid), std(fPeaks.chrom_centroid),...
            mean(fPeaks.chrom_variance), std(fPeaks.chrom_variance), ...
            mean(fPeaks.ms_Accu_mass_1), std(fPeaks.ms_Accu_mass_1), ...
            mean(fPeaks.ms_Accu_mass_2), std(fPeaks.ms_Accu_mass_2), ...
            mean(fPeaks.ms_variance), std(fPeaks.ms_variance), ...
            std(fPeaks.chrom_area)/mean(fPeaks.chrom_area)];
        countMe = countMe+1;


    end

    MrgPeakList(mList, :) = [];
end

SummaryNodes = array2table(SNodes);
SummaryNodes.Properties.VariableNames = {'idMF'; 'sizeMF'; 'Mean_Area';  ...
    'Mean_Chrom_PeakMaxima'; 'Std_Chrom_PeakMaxima';...
    'Mean_Chrom_PeakCenter'; 'Std_Chrom_PeakCenter'; ...
    'Mean_Chrom_PeakVariance'; 'Std_Chrom_PeakVariance'; ...
    'Mean_MS_AccuMass1'; 'Std_MS_AccuMass1'; ...
    'Mean_MS_AccuMass2'; 'Std_MS_AccuMass2'; ...
    'Mean_MS_PeakVariance'; 'Std_MS_PeakVariance'; ...
    'CV_Area'};

if ~exist(fullfile(obj.Path2Project, 'myProjectPeakLists'), 'dir') == 7
    mkdir(fullfile(obj.Path2Project, 'myProjectPeakLists'))
end

file = table;
file.ID = FilesID';
file.Name = FilesUsed';

listName = tempname(fullfile(obj.Path2Project, 'myProjectPeakLists'));
if exist(fullfile(obj.Path2Project, 'myProjectPeakLists'), 'dir') ~= 7
    mkdir(fullfile(obj.Path2Project, 'myProjectPeakLists'));
end

myMrgPeakList = {};
myMrgPeakList.Name = listName;
myMrgPeakList.AllFeatures = MrgPeakList;
myMrgPeakList.MergedFeatures = Nodes;
myMrgPeakList.Summary = SummaryNodes;
myMrgPeakList.AdditionalInformation.options = options;
myMrgPeakList.AdditionalInformation.ListFinnee = file;
myMrgPeakList.ComputingTime = toc(tStart);
myMrgPeakList.AdditionalInformation.Data = datetime;

obj.FeaturesLists{end+1}.Original.Name = myMrgPeakList.Name;
save(myMrgPeakList.Name, 'myMrgPeakList')
myProject = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')


