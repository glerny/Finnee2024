function [obj, myMrgPeakList] = mkMergedPeakList(obj, Tag2Fin, align, options)

%% INTRODUCTION
if align
    AlignMoi = true;

else
    AlignMoi = false;

end

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
    options.maxPPM              = 4;
    options.maxTm               = 0.01;

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
    options.Merge.MahalThr      = 5;
    options.Merge.SigTm         = 3;
    options.Merge.minRepl       = 10;
    options.Merge.S2N           = [3, 10, 25];
    options.Merge.method        = "method2";
    options.Merge.threshold     = 0.05;
    options.Merge.minOcc4Stat   = 2;
    options.Merge.MahalThr      = 5;
    options.Merge.LocMaxThre    = max(0.15*options.Merge.minRepl(1)*0.25, 0.5);


    options.method1.Dmz = 10; %ppm
    options.method1.Tm = 0.02; %min
    options.method1.outliers = 5; %min


    options.method2.Dmz = 10; %ppm
    options.method2.Tm = 0.02; %min


end
optionsGMM = statset('Display', 'off', 'MaxIter', 1000);
inputList = table();
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
            error("")
            %             myFile = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], Dataset4Peaklist, 'myROIs.mat');
            %             LMF = load(myFile);
            %             myROIs = LMF.myROIs;
            %
            %             myPeakList = mkPeakList(myROIs, 'analysis_1', options);
            %             save(IdfPT, 'myPeakList')
        end
    end

end

for ii = 1:numel(FilesUsed)
    IdF = FilesID(ii);
    Dataset4Peaklist = ['Dataset', num2str(obj.Summary.Dataset4Quali(IdF))];
    IdfPT = fullfile(obj.Summary.FolderID{IdF}, [obj.Summary.FileID{IdF}, '.fin'], Dataset4Peaklist, 'myPeakList.mat');
    load(IdfPT, "myPeakList");
    myPeakList = sortrows(myPeakList,"chrom_Time_IA","ascend");

    if AlignMoi
        correctMe = obj.Summary.AlignMe{IdF};
        myPeakList.chrom_Time_IA = myPeakList.chrom_Time_IA + ...
            polyval(correctMe{1}, myPeakList.chrom_Time_IA, correctMe{2}, correctMe{3});
        myPeakList = sortrows(myPeakList,"chrom_Time_IA","ascend");

    end

    myPeakList.Id2Dts = ii*ones(size(myPeakList, 1), 1);
    inputList = [inputList; myPeakList];
end

%% CREATE FEATURES FROM SUPERPOSING 3D PEAKS
switch  options.Merge.method
    case "method1"
        Nodes = {};
        countMe = 1;
        Impur = {};
        mPeaks = {};

        % STep 1: Quick merge
        OutputList = MyRecursive(inputList, options.method2.Dmz, options.method2.Tm, options.Merge.minRepl, {});
        for ii = 1:numel(OutputList)
            fPeaks = OutputList{ii};

            if numel(unique(fPeaks.Id2Dts)) == numel(fPeaks.Id2Dts)
                fPeaks(isoutlier(fPeaks.chrom_Time_IA, 'ThresholdFactor', options.method1.outliers) |...
                    isoutlier(fPeaks.ms_Accu_mass_1, 'ThresholdFactor', options.method1.outliers), :) = [];

                if height(fPeaks) >= options.Merge.minRepl
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

            else
                mPeaks{end + 1} = fPeaks;

            end
        end

        for ii = 1:numel(mPeaks)
            fPeaks = mPeaks{ii};


            loop = 1;

            while loop < 10

                fPeaks(isoutlier(fPeaks.chrom_Time_IA, 'ThresholdFactor', 3) |...
                    isoutlier(fPeaks.ms_Accu_mass_1, 'ThresholdFactor', 3), :) = [];

                if height(fPeaks) < options.Merge.minRepl; break; end
                if numel(unique(fPeaks.Id2Dts)) == numel(fPeaks.Id2Dts); break; end

                loop = loop + 1;
            end

            if height(fPeaks) >= options.Merge.minRepl & numel(unique(fPeaks.Id2Dts)) == numel(fPeaks.Id2Dts)
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


    case 'method2'
        Nodes = {};
        countMe = 1;
        Impur = {};

        OutputList = MyRecursive(inputList, options.method2.Dmz, options.method2.Tm, options.Merge.minRepl, {});

        for ii = 1:numel(OutputList)
            fPeaks = OutputList{ii};
            Xedges = min(fPeaks.chrom_Time_IA)-std(fPeaks.chrom_Time_IA):std(fPeaks.chrom_Time_IA)/2:max(fPeaks.chrom_Time_IA)+std(fPeaks.chrom_Time_IA);
            Yedges = min(fPeaks.ms_Accu_mass_1)-std(fPeaks.ms_Accu_mass_1):std(fPeaks.ms_Accu_mass_1)/2:max(fPeaks.ms_Accu_mass_1)+std(fPeaks.ms_Accu_mass_1);

            N = histcounts2(fPeaks.chrom_Time_IA, fPeaks.ms_Accu_mass_1, Xedges, Yedges);
            h = sgsdf_2d(-2:2, -2:2, 2, 2, 0);
            Data = filter2(h, N, 'same');
            Data(Data < 0) = 0;

            LocMax = [];

            for jj = 2:size(Data, 1)-1
                for kk = 2:size(Data, 2)-1
                    Int = [max(1, jj-2) min(jj+2, size(Data, 1)) ...
                        max(1, kk-2) min(kk+2, size(Data, 2))];
                    if Data(jj, kk) == max(Data(Int(1):Int(2), Int(3):Int(4)), [], 'all', 'omitnan')
                        LocMax(end+1, :) = [Data(jj, kk), jj, kk];

                    end
                end
            end
            LocMax(LocMax(:,1) < options.Merge.LocMaxThre, :) = [];

            XY = [fPeaks.chrom_Time_IA, fPeaks.ms_Accu_mass_1];
            if height(LocMax) <= 1 | numel(fPeaks.Id2Dts) - numel(unique(fPeaks.Id2Dts)) <= round(0.01*numel(unique(fPeaks.Id2Dts)))
                gmd = gmdistribution(mean(XY), cov(XY));
                % fPeaks(mahal(gmd, XY) > options.Merge.MahalThr, :) = [];
                XY = [fPeaks.chrom_Time_IA, fPeaks.ms_Accu_mass_1];
                gmd = gmdistribution(mean(XY), cov(XY));

                ListDoublons = find(histcounts(fPeaks.Id2Dts, 1:numel(FilesUsed)+1) > 1);
                Id2Del = [];
                for jj = 1:numel(ListDoublons)
                    IdD = find(fPeaks.Id2Dts == ListDoublons(jj));
                    mahald2 = mahal(gmd, [fPeaks.chrom_Time_IA(IdD), fPeaks.ms_Accu_mass_1(IdD)]);
                    [~, IdM] = min(mahald2);
                    IdD(IdM) = [];
                    Id2Del = [Id2Del; IdD];

                end

                fPeaks(Id2Del, :) = [];
                if height(fPeaks) > options.Merge.minRepl
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

            else
                IdGmd = kmeans(XY, height(LocMax), 'Start', [Xedges(LocMax(:, 2))' Yedges(LocMax(:, 3))']);
                % gmd = fitgmdist(XY, height(LocMax), 'CovarianceType','diagonal','SharedCovariance', true, 'Start', gmd);
                % IdGmd = gmd.cluster(XY);

                for jj = 1:max(IdGmd)
                    gPeaks = fPeaks(IdGmd==jj, :);
                    if height(gPeaks) <= options.Merge.minRepl, continue; end

                    XYg = [gPeaks.chrom_Time_IA, gPeaks.ms_Accu_mass_1];
                    gmd = gmdistribution(mean(XYg), cov(XYg));
                    % gPeaks(mahal(gmd, XYg) > options.Merge.MahalThr, :) = [];
                    XYg = [gPeaks.chrom_Time_IA, gPeaks.ms_Accu_mass_1];
                    gmd = gmdistribution(mean(XYg), cov(XYg));

                    ListDoublons = find(histcounts(gPeaks.Id2Dts, 1:numel(FilesUsed)+1) > 1);
                    Id2Del = [];
                    for jj = 1:numel(ListDoublons)
                        IdD = find(gPeaks.Id2Dts == ListDoublons(jj));
                        mahald2 = mahal(gmd, [fPeaks.chrom_Time_IA(IdD), gPeaks.ms_Accu_mass_1(IdD)]);
                        [~, IdM] = min(mahald2);
                        IdD(IdM) = [];
                        Id2Del = [Id2Del; IdD];

                    end

                    gPeaks(Id2Del, :) = [];
                    if height(gPeaks) > options.Merge.minRepl
                        Nodes{countMe} = gPeaks;

                        SNodes(countMe, :) = [countMe, size(gPeaks, 1), mean(gPeaks.chrom_area), ...
                            mean(gPeaks.chrom_Time_IA), std(gPeaks.chrom_Time_IA),...
                            mean(gPeaks.chrom_centroid), std(gPeaks.chrom_centroid),...
                            mean(gPeaks.chrom_variance), std(gPeaks.chrom_variance), ...
                            mean(gPeaks.ms_Accu_mass_1), std(gPeaks.ms_Accu_mass_1), ...
                            mean(gPeaks.ms_Accu_mass_2), std(gPeaks.ms_Accu_mass_2), ...
                            mean(gPeaks.ms_variance), std(gPeaks.ms_variance), ...
                            std(gPeaks.chrom_area)/mean(gPeaks.chrom_area)];
                        countMe = countMe+1;

                    end
                end

            end
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

end

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
myMrgPeakList.AllFeatures = inputList;
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

    function OutputList = MyRecursive(inputList, Dmz, tm, minsize, OutputList)
        ThisTag = true;

        inputList = sortrows(inputList, "ms_Accu_mass_1");
        MySplits = ((inputList.ms_Accu_mass_1(2:end)-inputList.ms_Accu_mass_1(1: end-1))./...
            ((inputList.ms_Accu_mass_1(1: end-1)+ inputList.ms_Accu_mass_1(2:end))/2))*1000000;
        MySplits(:, 2) = MySplits(:, 1) > Dmz;
        if any(MySplits(:, 2)), ThisTag = false; end

        list1 = unique([0; find(MySplits(:, 2)); height(inputList)]);
        for iMR = 1:numel(list1)-1
            SplitList1 = inputList(list1(iMR)+1:list1(iMR+1), :);

            if height(SplitList1) >= minsize
                SplitList1 = sortrows(SplitList1, "chrom_Time_IA");
                My2nSplits = (SplitList1.chrom_Time_IA(2:end) - SplitList1.chrom_Time_IA(1:end-1));
                My2nSplits(:, 2) = My2nSplits(:, 1) > tm;
                if any(My2nSplits(:, 2)), ThisTag = false; end

                list2 = unique([0; find(My2nSplits(:, 2)); height(SplitList1)]);

                for jMR = 1:numel(list2)-1
                    SplitList2 = SplitList1(list2(jMR)+1:list2(jMR+1), :);

                    if height(SplitList2) >= minsize
                        if ~ ThisTag
                            OutputList = MyRecursive(SplitList2, Dmz, tm, minsize, OutputList);

                        else
                            OutputList{end+1} = inputList;

                        end
                    end
                end
            end
        end
    end

end


