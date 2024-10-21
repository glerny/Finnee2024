function [obj, myMrgPeakList] = doQuantitative(obj, Tag2Fin, Id2PeakList, replace, timeAlign, options)

%% INTRODUCTION
if isempty(options)

    options.ForROIs.SigMZ         = 4;
    options.ForROIs.SigTm         = 20;
    options.ForROIs.TgtDataset    = 'Dataset2';
    options.ForROIs.IdDataset     = 2;
    options.ForROIs.Folder        = 'ROI4Quant';
    options.ForROIs.alfa          = 0.95;

    options.Analysis.sigMZ        = 1;
    options.Analysis.minTime      = 5;
    options.Analysis.WindowTime   = 2;
    options.Analysis.WindowMZ     = 2;
    options.Analysis.sgfd         = 'average';
    options.Analysis.critRes      = 0.6;
    options.Analysis.DynRange     = 0.05;
    
    options.Merge.maxPPM          = 10;
    options.Merge.maxTm           = 0.02;
    options.Merge.maxRs           = 10;
    options.Merge.S2N             = [3, 3];
    options.Merge.n2n             = [0.9, 0.75];
    options.Merge.minRepl         = [10, 8];
    options.Merge.MahalThr        = 3;
    options.Merge.LocMaxThre      = max(0.15*options.Merge.minRepl(1)*0.25, 0.5);
    options.Merge.MahalThres      = 7;

    options.Doublons.CritRes      = [1.0, 0.5];

end

%% INITIALISATION %%

load(obj.FeaturesLists{Id2PeakList}.Original.Name);
mySummary = myMrgPeakList.Summary;

%% END INITIALISATION %%


target = table();
target.ID = mySummary.idMF;
target.TargetTime = mySummary.Mean_Chrom_PeakCenter;
target.StdDev_time = mySummary.Std_Chrom_PeakMaxima;
target.PeakStdDev_time = sqrt(mySummary.Mean_Chrom_PeakVariance);
target.TimeMin = mySummary.Mean_Chrom_PeakCenter ... % center of the peak
    -3*mySummary.Std_Chrom_PeakCenter...             % variation of the center of the peak
    -options.ForROIs.SigTm*sqrt(mySummary.Mean_Chrom_PeakVariance);
target.TimeMax = mySummary.Mean_Chrom_PeakCenter ...
    +3*mySummary.Std_Chrom_PeakCenter...
    +options.ForROIs.SigTm*sqrt(mySummary.Mean_Chrom_PeakVariance);
target.Targetmz = mySummary.Mean_MS_AccuMass2;
target.StdDev_MS = mySummary.Std_MS_AccuMass1;
target.PeakStdDev_MS = sqrt(mySummary.Mean_MS_PeakVariance);
target.mzMin = mySummary.Mean_MS_AccuMass2 ...          % center of the peak
    -3*mySummary.Std_MS_AccuMass2...                    % variation of the center of the peak
    -options.ForROIs.SigMZ*sqrt(mySummary.Mean_MS_PeakVariance);
target.mzMax = mySummary.Mean_MS_AccuMass2 ...
    +3*mySummary.Std_MS_AccuMass2...
    +options.ForROIs.SigMZ*sqrt(mySummary.Mean_MS_PeakVariance);
%% CHECK EXISTSING ROI
% ATTENTION: Check if ROI exist and if ~replace

myMapFiles = {}; myMapFiles.Id2ii = []; myMapFiles.Id = [];
myMapFiles.OriFile = []; myMapFiles.FinneeFile = [];  myMapFiles.toDO = []; myMapFiles.CorrectMe = {};
myMapFiles.TgtFolder = [];  

for ii = 1:size(obj.Summary, 1)
    if strcmp(obj.Summary.FinneeType{ii}, Tag2Fin)
        myTgtFolder = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], ...
            options.ForROIs.TgtDataset);
        myMapFiles.Id2ii(end+1, 1) = ii;
        myMapFiles.Id(end+1, 1) = numel(myMapFiles.Id2ii);
        myMapFiles.OriFile{end+1, 1} = obj.Summary.FileID{ii};
        myMapFiles.FinneeFile{end+1, 1} = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], 'myFinnee.mat');

        if timeAlign
            myMapFiles.CorrectMe{end+1, 1} = obj.Summary.AlignMe{ii};

        else
            myMapFiles.CorrectMe{end+1, 1} = {};

        end


        if ~(exist(fullfile(myTgtFolder, options.ForROIs.Folder), "dir") == 7)
            myMapFiles.toDO(end+1, 1) = true;

        elseif replace
            myMapFiles.toDO(end+1, 1) = true;

        else
            myMapFiles.toDO(end+1, 1) = false;

        end
        myMapFiles.TgtFolder{end+1, 1} = myTgtFolder;

    end
end
myMapFiles = struct2table(myMapFiles);

%% MAKE THE ROIs

IdFile = find(myMapFiles.toDO);
if ~isempty(IdFile)
    parfor ii = 1:numel(IdFile)
        myFile =  myMapFiles.FinneeFile{ii};
        CFN = load(myFile);
        CFN = CFN.myFinnee;

        if timeAlign
            corTarget = target;

        else
            corTarget = target;
            for jj = 1:size(corTarget, 1)
                corTarget.TargetTime(jj) = reversepolyval(corTarget.TargetTime(jj), obj.Summary.AlignMe{ii});
                corTarget.TimeMin(jj) = reversepolyval(corTarget.TimeMin(jj), obj.Summary.AlignMe{ii});
                corTarget.TimeMax(jj) = reversepolyval(corTarget.TimeMax(jj), obj.Summary.AlignMe{ii});

            end

        end

        myTgtFolder = fullfile(myMapFiles.TgtFolder{ii}, options.ForROIs.Folder);

        CFN.mkMnROI(options.ForROIs.IdDataset, corTarget, [0 0 0], myTgtFolder);
        disp([obj.Summary.FileID{ii}, ' done!'])

    end
end

for ii = 1:numel(myMapFiles.Id2ii)
    myTgtFolder = fullfile(myMapFiles.TgtFolder{ii}, 'myROIs.mat');
    MR = load(myTgtFolder);
    myROIs{ii} = MR.myROIs;

end

variable2_names_types = [["chrom_apex", "double", "%.4f", "m/z"]; ...
    ["chrom_intensity_apex", "double", "%.0f", ""]; ...
    ["chrom_area", "double", "%.2f", "Arb. units"]; ...
    ["chrom_centroid", "double", "%.4f", ""]; ...
    ["chrom_variance", "double", "%.2e", ""];...
    ["chrom_assym", "double", "%.2e", ""];...
    ["Peak_start", "double", "%.2e", ""];...
    ["Peak_end", "double", "%.2e", ""];...
    ["Accu_mass", "double", "%.4f", ""]; ...
    ["ms_variance", "double", "%.4f", ""]; ...
    ["ms_assym", "double", "%.4f", ""]; ...
    ["ms_Peak_start", "double", "%.2e", ""];...
    ["ms_Peak_end", "double", "%.2e", ""];...
    ["Signal", "double", "%.4f", ""]; ...
    ["Noise", "double", "%.4f", ""]; ...
    ["n2n", "double", "%.4f", ""]; ...
    ["Id2Tgt", "double", "%.0f", ""]];
IsEmpty = table('Size', [0, size(variable2_names_types,1)],...
    'VariableNames', variable2_names_types(:,1),...
    'VariableTypes', variable2_names_types(:,2));
IsEmpty.Properties.VariableUnits = variable2_names_types(:,4);
IsEmpty.Properties.VariableDescriptions = variable2_names_types(:,3);
IsEmpty(1, :) = {NaN, NaN, NaN, NaN,  NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};

Area = {}; Area.ID = target.ID; Area.Values = [];
Features = {}; Features.ID = target.ID; Features.Values = {};
Summary = {};  Summary.ID = [];
Summary.mean_IApex = []; Summary.Std_IApex = [];
Summary.mean_Apex = []; Summary.Std_Apex = [];
Summary.mean_ctr = []; Summary.Std_ctr = [];
Summary.mean_ChrV = []; Summary.Std_ChrV = [];
Summary.mean_Accu = []; Summary.Std_Accu = [];
Summary.mean_MZV = []; Summary.Std_MZV = [];
Summary.chr_Assym = []; Summary.MS_Assym = [];
Summary.n = [];

for ii = 1:numel(target.ID)

    ii
    AN = table();
    for jj = 1:numel(myROIs)

        myROI = ROI(myROIs{jj}, target.ID(ii), target(ii, :));

        myROIAnalysis = myROI.analysis_2(options.Analysis);
        if ~isempty(myROIAnalysis)
            myROIAnalysis.Id2Tgt = jj*ones(size(myROIAnalysis.chrom_apex));

        else
            continue

        end

        if timeAlign
            myROIAnalysis.chrom_apex = myROIAnalysis.chrom_apex + polyval(myMapFiles.CorrectMe{jj}{1}, myROIAnalysis.chrom_apex, myMapFiles.CorrectMe{jj}{2}, myMapFiles.CorrectMe{jj}{3});
            myROIAnalysis.chrom_centroid = myROIAnalysis.chrom_centroid + polyval(myMapFiles.CorrectMe{jj}{1}, myROIAnalysis.chrom_centroid, myMapFiles.CorrectMe{jj}{2}, myMapFiles.CorrectMe{jj}{3});

        end
        AN = [AN; myROIAnalysis];
    end

    % FILTER DAT
    AN(find(isnan(AN.chrom_apex)), :) = [];
    if size(AN, 1) <  options.Merge.minRepl(1)
        continue

    end

    %% CLUSTER DATA
    Nodes = {}; 
    StdMz = std(AN.Accu_mass);
    StdTm = std(AN.chrom_apex);
    nLoop = 1;

    while 1
        stepIn = false;
        Xedges = min(AN.chrom_apex)-StdTm:StdTm/2:max(AN.chrom_apex)+StdTm;
        Yedges = min(AN.Accu_mass)-StdMz:StdMz/2:max(AN.Accu_mass)+StdMz;

        N = histcounts2(AN.chrom_apex, AN.Accu_mass, Xedges, Yedges);
        h = sgsdf_2d(-2:2, -2:2, 2, 2, 0);
        Data = filter2(h, N, 'same');
        Data(Data < 0) = 0;
        Nodes = {};

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
        LocMax(LocMax(:,1) <  options.Merge.LocMaxThre, :) = [];

        if height(LocMax) == 0
            Nodes{end + 1} = AN(IdGmd == jj, :);
            break
        end

        XY = [AN.chrom_apex, AN.Accu_mass];
        IdGmd = kmeans(XY, height(LocMax), 'Start', [Xedges(LocMax(:, 2))' Yedges(LocMax(:, 3))']);

        for jj = 1:max(IdGmd)
            if numel(AN.Id2Tgt(IdGmd == jj)) - numel(unique(AN.Id2Tgt(IdGmd == jj))) > max(round(0.05*numel(unique(AN.Id2Tgt(IdGmd == jj)))), 1) & nLoop <= 5
                stepIn = true;
                nLoop = nLoop + 1;
                break

            else
                Nodes{end + 1} = AN(IdGmd == jj, :);
            end
        end

        if stepIn
            StdMz = StdMz/2;
            StdTm = StdTm/2;
        
        else
            break

        end

    end

    NNodes = {};
    SNodes = [];
    countMe = 1;

    for jj = 1:numel(Nodes)
        fPeaks = Nodes{jj};

        if height(fPeaks) < options.Merge.minRepl(1), continue; end

        XY = [fPeaks.chrom_apex, fPeaks.Accu_mass];
        gmd = gmdistribution(mean(XY), cov(XY));
        fPeaks(mahal(gmd, XY) > options.Merge.MahalThres, :) = [];
        XY = [fPeaks.chrom_apex, fPeaks.Accu_mass];
        gmd = gmdistribution(mean(XY), cov(XY));

        ListDoublons = find(histcounts(fPeaks.Id2Tgt, 1:max(fPeaks.Id2Tgt)+1) > 1);
        Id2Del = [];
        for kk = 1:numel(ListDoublons)
            IdD = find(fPeaks.Id2Tgt == ListDoublons(kk));
            mahald2 = mahal(gmd, [fPeaks.chrom_apex(IdD), fPeaks.Accu_mass(IdD)]);
            [~, IdM] = min(mahald2);
            IdD(IdM) = [];
            Id2Del = [Id2Del; IdD];

        end
        fPeaks(Id2Del, :) = [];

        if height(fPeaks) > options.Merge.minRepl(1)

            NNodes{countMe} = fPeaks;
            SNodes(countMe, :) = [countMe, size(fPeaks, 1), mean(fPeaks.chrom_area), ...
                mean(fPeaks.chrom_apex), std(fPeaks.chrom_apex),...
                mean(fPeaks.chrom_centroid), std(fPeaks.chrom_centroid),...
                mean(fPeaks.chrom_variance), std(fPeaks.chrom_variance), ...
                mean(fPeaks.Accu_mass), std(fPeaks.Accu_mass), ...
                mean(fPeaks.ms_variance), std(fPeaks.ms_variance)];
            countMe = countMe+1;
        end
    end

    % END CLUSTER

    if ~isempty(SNodes)
        [~, idX] = min(...
            ((SNodes(:, 10)-target.Targetmz(ii))./SNodes(:, 11)).^2+...
            ((SNodes(:, 4)-target.TargetTime(ii))./SNodes(:, 5)).^2);
        cAN = NNodes{idX};

        if numel(cAN.Id2Tgt) ~= numel(unique(cAN.Id2Tgt))
            disp("pp")
        end

        if height(cAN) > options.Merge.minRepl(1)
            V = (1:numel(obj.ListOfFiles))';
            [~, ia, ib] = intersect(V, cAN.Id2Tgt);
            V(ia, 2:4) = [cAN.chrom_area(ib), cAN.chrom_centroid(ib), cAN.Accu_mass(ib)];
            IdROI = target.ID(ii);
            Area.Values(Area.ID == IdROI, :) = V(:, 2)';
            Features.Values{Features.ID == IdROI} = cAN;
            Summary.ID(end+1, 1) = IdROI;
            Summary.mean_IApex(end+1, 1) = mean(cAN.chrom_intensity_apex, 'omitnan');
            Summary.Std_IApex(end+1, 1)  = std(cAN.chrom_intensity_apex, [], 'omitnan');
            Summary.mean_Apex(end+1, 1) = mean(cAN.chrom_apex, 'omitnan');
            Summary.Std_Apex(end+1, 1)  = std(cAN.chrom_apex, [], 'omitnan');
            Summary.mean_ctr(end+1, 1) = mean(cAN.chrom_centroid, 'omitnan');
            Summary.Std_ctr(end+1, 1) = std(cAN.chrom_centroid, [], 'omitnan');
            Summary.mean_ChrV(end+1, 1) = mean(cAN.chrom_variance, 'omitnan');
            Summary.Std_ChrV(end+1, 1) = std(cAN.chrom_variance, [], 'omitnan');
            Summary.chr_Assym(end+1, 1) = mean(cAN.chrom_assym, 'omitnan');
            Summary.mean_Accu(end+1, 1) =  mean(cAN.Accu_mass, 'omitnan');
            Summary.Std_Accu(end+1, 1) =  std(cAN.Accu_mass, [], 'omitnan');
            Summary.mean_MZV(end+1, 1) =  mean(cAN.ms_variance, 'omitnan');
            Summary.Std_MZV(end+1, 1) = std(cAN.ms_variance, [], 'omitnan');
            Summary.MS_Assym(end+1, 1) = mean(cAN.ms_assym, 'omitnan');
            Summary.n(end+1, 1) = height(cAN);
        end
    end
end
Summary = struct2table(Summary);

%% PROVISORY
assignin('base', 'target', target)
assignin('base', 'options', options)
assignin('base', 'Features', Features)
assignin('base', 'Area', Area)
assignin('base', 'Summary', Summary)

%% CHECK FOR DOUBLONS
while 1
    doStop = true;
    StdROIs = [Summary.mean_Apex Summary.Std_Apex  Summary.mean_Accu Summary.Std_Accu];

    for ii = 1:height(StdROIs)
        RS = abs((StdROIs(:, 1) - StdROIs(ii, 1))./(2*((StdROIs(:, 2) + StdROIs(ii, 2)))));
        RS(:, 2) = abs((StdROIs(:, 3) - StdROIs(ii, 3))./(2*((StdROIs(:, 4) + StdROIs(ii, 4)))));
        RS(:, 3) = sqrt(RS(:, 1).^2 + RS(:, 2).^2);

        IdDoublons = find(RS(:, 3) < options.Doublons.CritRes(1));
        if numel(IdDoublons) > 1

            TAN = []; NAN = [];
            [val, IdFeature] = intersect(Features.ID, Summary.ID(IdDoublons));

            if numel(val) ~= numel(IdDoublons)
                error("line 339");
            end

            for jj = 1:numel(IdFeature)
                TAN = [TAN; Features.Values{IdFeature(jj)}];
            end
            myDts = unique(TAN.Id2Tgt);

            for kk = 1:numel(myDts)
                sAN = TAN(TAN.Id2Tgt == myDts(kk), :);

                if height(sAN) > 1
                    IdStart = 1;
                    while 1
                        RR = abs((sAN.chrom_apex(IdStart) - sAN.chrom_apex)./(2*(sqrt(sAN.chrom_variance(IdStart)) + sqrt(sAN.chrom_variance))));
                        RR(:, 2) = abs((sAN.Accu_mass(IdStart) - sAN.Accu_mass)./(2*(sqrt(sAN.ms_variance(IdStart)) + sqrt(sAN.ms_variance))));
                        RR(:, 3) = sqrt(RR(:, 1).^2 + RR(:, 2).^2);

                        Id2Merge = find(RR(:, 3) < options.Doublons.CritRes(2));
                        if numel(Id2Merge) == 1
                            IdStart = IdStart + 1;

                        else
                            [~, IM] = max(sAN.chrom_area(Id2Merge));
                            sAN(Id2Merge~=Id2Merge(IM), :) = [];

                        end
                        if IdStart+1 > height(sAN), break; end
                    end
                end
                NAN = [NAN; sAN];
            end

            StdMz = std(NAN.Accu_mass);
            StdTm =  std(NAN.chrom_apex);
            nLoop = 1;

            while 1
                stepIn = false;
                Xedges = min(NAN.chrom_apex)-StdTm:StdTm/2:max(NAN.chrom_apex)+StdTm;
                Yedges = min(NAN.Accu_mass)-StdMz:StdMz/2:max(NAN.Accu_mass)+StdMz;

                N = histcounts2(NAN.chrom_apex, NAN.Accu_mass, Xedges, Yedges);
                h = sgsdf_2d(-2:2, -2:2, 2, 2, 0);
                Data = filter2(h, N, 'same');
                Data(Data < 0) = 0;
                NNodes = {};

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
                LocMax(LocMax(:,1) <  options.Merge.LocMaxThre, :) = [];

                if height(LocMax) == 0
                    NNodes{end + 1} = NAN(IdGmd == jj, :);
                    break
                end

                XY = [NAN.chrom_apex, NAN.Accu_mass];
                IdGmd = kmeans(XY, height(LocMax), 'Start', [Xedges(LocMax(:, 2))' Yedges(LocMax(:, 3))']);

                for jj = 1:max(IdGmd)
                    if numel(NAN.Id2Tgt(IdGmd == jj)) - numel(unique(NAN.Id2Tgt(IdGmd == jj))) > max(round(0.05*numel(unique(NAN.Id2Tgt(IdGmd == jj)))), 1) & nLoop <= 5
                        stepIn = true;
                        nLoop = nLoop + 1;
                        break

                    else
                        NNodes{end + 1} = NAN(IdGmd == jj, :);
                    end
                end

                if stepIn
                    StdMz = StdMz/2;
                    StdTm = StdTm/2;

                else
                    break

                end

            end

            ToKtoD = IdDoublons;
            sTarget = target(Summary.ID(IdDoublons), :);
            ToKtoD(:, 2) = NaN;
            ToKtoD(:, 3) = Summary.ID(IdDoublons);
            for jj = 1:numel(NNodes)
                pAN = NNodes{jj};

                if height(pAN) <= options.Merge.minRepl(1), continue; end

                XY = [pAN.chrom_apex, pAN.Accu_mass];
                gmd = gmdistribution(mean(XY), cov(XY));
                pAN(mahal(gmd, XY) > options.Merge.MahalThres, :) = [];
                XY = [pAN.chrom_apex, pAN.Accu_mass];
                gmd = gmdistribution(mean(XY), cov(XY));

                LD = find(histcounts(pAN.Id2Tgt, 1:max(pAN.Id2Tgt)+1) > 1);
                ID = [];
                for kk = 1:numel(LD)
                    IdD = find(pAN.Id2Tgt == LD(kk));
                    mahald2 = mahal(gmd, [pAN.chrom_apex(IdD), pAN.Accu_mass(IdD)]);
                    [~, IdM] = min(mahald2);
                    IdD(IdM) = [];
                    ID = [ID; IdD];

                end
                pAN(ID, :) = [];

                if height(pAN) > options.Merge.minRepl(1)
                    [~, IdT] = min(((mean(pAN.chrom_apex)-sTarget.TargetTime)/std(pAN.chrom_apex)).^2+((mean(pAN.Accu_mass)-sTarget.Targetmz)/std(pAN.Accu_mass)).^2);

                    if isnan(ToKtoD(IdT, 2))
                        IdTgt = ToKtoD(IdT, 1);
                        V = (1:sum(strcmp(obj.Summary.FinneeType, Tag2Fin)))';
                        [~, ia, ib] = intersect(V, pAN.Id2Tgt);
                        V(ia, 2:4) = [pAN.chrom_area(ib), pAN.chrom_centroid(ib), pAN.Accu_mass(ib)];
                        Area.Values(Area.ID ==  Summary.ID(IdTgt), :) = V(:, 2)';
                        Features.Values{Features.ID ==  Summary.ID(IdTgt)} = pAN;
                        Summary.mean_IApex(IdTgt) = mean(pAN.chrom_intensity_apex, 'omitnan');
                        Summary.Std_IApex(IdTgt)  = std(pAN.chrom_intensity_apex, [], 'omitnan');
                        Summary.mean_Apex(IdTgt) = mean(pAN.chrom_apex, 'omitnan');
                        Summary.Std_Apex(IdTgt)  = std(pAN.chrom_apex, [], 'omitnan');
                        Summary.mean_ctr(IdTgt) = mean(pAN.chrom_centroid, 'omitnan');
                        Summary.Std_ctr(IdTgt) = std(pAN.chrom_centroid, [], 'omitnan');
                        Summary.mean_ChrV(IdTgt) = mean(pAN.chrom_variance, 'omitnan');
                        Summary.Std_ChrV(IdTgt) = std(pAN.chrom_variance, [], 'omitnan');
                        Summary.chr_Assym(IdTgt) = mean(pAN.chrom_assym, 'omitnan');
                        Summary.mean_Accu(IdTgt) =  mean(pAN.Accu_mass, 'omitnan');
                        Summary.Std_Accu(IdTgt) =  std(pAN.Accu_mass, [], 'omitnan');
                        Summary.mean_MZV(IdTgt) =  mean(pAN.ms_variance, 'omitnan');
                        Summary.Std_MZV(IdTgt) = std(pAN.ms_variance, [], 'omitnan');
                        Summary.MS_Assym(IdTgt) = mean(pAN.ms_assym, 'omitnan');
                        Summary.n(IdTgt) = height(pAN);
                        ToKtoD(IdT, 2) = 1;

                    else
                        IdTgt = ToKtoD(IdT, 1);
                        if Summary.n(IdTgt) < height(pAN)
                           V = (1:sum(strcmp(obj.Summary.FinneeType, Tag2Fin)))';
                            [~, ia, ib] = intersect(V, pAN.Id2Tgt);
                            V(ia, 2:4) = [pAN.chrom_area(ib), pAN.chrom_centroid(ib), pAN.Accu_mass(ib)];
                            Area.Values(Area.ID ==  Summary.ID(IdTgt), :) = V(:, 2)';
                            Features.Values{Features.ID ==  Summary.ID(IdTgt)} = pAN;
                            Summary.mean_IApex(IdTgt) = mean(pAN.chrom_intensity_apex, 'omitnan');
                            Summary.Std_IApex(IdTgt)  = std(pAN.chrom_intensity_apex, [], 'omitnan');
                            Summary.mean_Apex(IdTgt) = mean(pAN.chrom_apex, 'omitnan');
                            Summary.Std_Apex(IdTgt)  = std(pAN.chrom_apex, [], 'omitnan');
                            Summary.mean_ctr(IdTgt) = mean(pAN.chrom_centroid, 'omitnan');
                            Summary.Std_ctr(IdTgt) = std(pAN.chrom_centroid, [], 'omitnan');
                            Summary.mean_ChrV(IdTgt) = mean(pAN.chrom_variance, 'omitnan');
                            Summary.Std_ChrV(IdTgt) = std(pAN.chrom_variance, [], 'omitnan');
                            Summary.chr_Assym(IdTgt) = mean(pAN.chrom_assym, 'omitnan');
                            Summary.mean_Accu(IdTgt) =  mean(pAN.Accu_mass, 'omitnan');
                            Summary.Std_Accu(IdTgt) =  std(pAN.Accu_mass, [], 'omitnan');
                            Summary.mean_MZV(IdTgt) =  mean(pAN.ms_variance, 'omitnan');
                            Summary.Std_MZV(IdTgt) = std(pAN.ms_variance, [], 'omitnan');
                            Summary.MS_Assym(IdTgt) = mean(pAN.ms_assym, 'omitnan');
                            Summary.n(IdTgt) = height(pAN);

                        end
                    end
                end
            end

            Id2Del = find(isnan(ToKtoD(:, 2)));
            Summary(ToKtoD(Id2Del, 1), :) = [];
            for il =1:numel(Id2Del)

                Area.Values(Area.ID == ToKtoD(Id2Del(il), 3), :) = NaN;
                Features.Values{Features.ID ==  ToKtoD(Id2Del(il), 3)} = table();
            end

            doStop = false;
            break;
        end
    end

    if doStop, break; end
end

if ~isfield(myMrgPeakList, 'TargetedAnalysis')
    myMrgPeakList.TargetedAnalysis = {};
end
myMrgPeakList.TargetedAnalysis{end+1}.Date = datetime('now');
myMrgPeakList.TargetedAnalysis{end}.options = options;
myMrgPeakList.TargetedAnalysis{end}.Target = target;
myMrgPeakList.TargetedAnalysis{end}.Summary =  Summary;
myMrgPeakList.TargetedAnalysis{end}.Area = Area;
myMrgPeakList.TargetedAnalysis{end}.Features = Features;
myMrgPeakList.TargetedAnalysis{end}.mapFiles = myMapFiles;
save(obj.FeaturesLists{Id2PeakList}.Original.Name, 'myMrgPeakList')


