function [obj, myMrgPeakList] = doQuantitative_small_single(obj, Tag2Fin, PeakList, replace, timeAlign, options)

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
    options.Merge.MahalThres      = 3;
    options.Doublons.CritRes      = [1, 0.5];

end



%% END INITIALISATION %%


target = table();
target.ID = PeakList.ID;
target.TargetTime = PeakList.chrom_Time_IA;
target.PeakStdDev_time = sqrt(PeakList.chrom_variance);
target.TimeMin = PeakList.chrom_centroid ... % center of the peak
    -options.ForROIs.SigTm*sqrt(PeakList.chrom_variance);
target.TimeMax = PeakList.chrom_centroid ... % center of the peak
    +options.ForROIs.SigTm*sqrt(PeakList.chrom_variance);
target.Targetmz = PeakList.ms_Accu_mass_2;
target.PeakStdDev_MS = sqrt(PeakList.ms_variance);
target.mzMin = PeakList.ms_Accu_mass_2 ...         
    -options.ForROIs.SigMZ*sqrt(PeakList.ms_variance);
target.mzMax = PeakList.ms_Accu_mass_2 ...         
    +options.ForROIs.SigMZ*sqrt(PeakList.ms_variance);
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
    for ii = 1:numel(IdFile)
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

        try

        CFN.mkMnROI(options.ForROIs.IdDataset, corTarget, [0 0 0], myTgtFolder);
        catch
            dip("pp")

        end
        disp([myMapFiles.FinneeFile{ii}, ' done!'])

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
    cAN = table();

    OutputList = MyRecursive(AN, options.Merge.maxPPM, options.Merge.maxTm, options.Merge.minRepl(1), {});
    if isempty(OutputList), continue; end

    TestMe = (1:numel(OutputList))';

    for jj = 1:numel(OutputList)
        TestMe(jj, 2) = (mean(OutputList{jj}.Accu_mass)- target.Targetmz(ii)).^2 +...
            (mean(OutputList{jj}.chrom_apex - target.TargetTime(ii))).^2;

    end
    [~, Idmin] = min(TestMe(:, 2));
    cAN = OutputList{Idmin};

    while numel(cAN.Id2Tgt) ~= numel(unique(cAN.Id2Tgt))
        [~, id] = max(((cAN.Accu_mass - target.Targetmz(ii))/std(cAN.Accu_mass)).^2 +...
            ((cAN.chrom_apex - target.TargetTime(ii))/std(cAN.chrom_apex)).^2);
        cAN(id, :) = [];

    end

    if height(cAN) >= options.Merge.minRepl(1)
        V = (1:sum(strcmp(obj.Summary.FinneeType, Tag2Fin)))';
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
Summary = struct2table(Summary);

%% PROVISORY
assignin('base', 'target', target)
assignin('base', 'options', options)
assignin('base', 'Features', Features)
assignin('base', 'Area', Area)
assignin('base', 'Summary', Summary)

Summary.NotDone = true(height(Summary), 1);
while 1

    myOutliers = isoutlier(Summary.Std_Accu./Summary.mean_Accu);
    [~, IdFeature] = intersect(Features.ID, Summary.ID(find(myOutliers & Summary.NotDone)));
    if isempty(IdFeature), break, end

    for ii = 1:numel(IdFeature)
        testedFeat = Features.Values{IdFeature(ii)};
        if height(testedFeat) == 2

            Id2Del = find(Summary.ID == IdFeature(ii));
            Summary(Id2Del, :) = [];
            Area.Values(IdFeature(ii), :) = NaN;
            Features.Values{IdFeature(ii)} = table();

        else
            XY = [testedFeat.Accu_mass, testedFeat.chrom_intensity_apex];
            gmd = gmdistribution(mean(XY), cov(XY));
            testedFeat(mahal(gmd, XY) > options.Merge.MahalThres, :) = [];

            if height(testedFeat) == 2

                Id2Del = find(Summary.ID == IdFeature(ii));
                Summary(Id2Del, :) = [];
                Area.Values(IdFeature(ii), :) = NaN;
                Features.Values{IdFeature(ii)} = table();

            else

                Id2Rpl = find(Summary.ID == IdFeature(ii));
                V = (1:sum(strcmp(obj.Summary.FinneeType, Tag2Fin)))';
                [~, ia, ib] = intersect(V, testedFeat.Id2Tgt);
                V(ia, 2:4) = [testedFeat.chrom_area(ib), testedFeat.chrom_centroid(ib), testedFeat.Accu_mass(ib)];
                Area.Values(IdFeature(ii), :) = V(:, 2)';
                Features.Values{IdFeature(ii)} = testedFeat;
                Summary.mean_IApex(Id2Rpl) = mean(testedFeat.chrom_intensity_apex, 'omitnan');
                Summary.Std_IApex(Id2Rpl)  = std(testedFeat.chrom_intensity_apex, [], 'omitnan');
                Summary.mean_Apex(Id2Rpl) = mean(testedFeat.chrom_apex, 'omitnan');
                Summary.Std_Apex(Id2Rpl)  = std(testedFeat.chrom_apex, [], 'omitnan');
                Summary.mean_ctr(Id2Rpl) = mean(testedFeat.chrom_centroid, 'omitnan');
                Summary.Std_ctr(Id2Rpl) = std(testedFeat.chrom_centroid, [], 'omitnan');
                Summary.mean_ChrV(Id2Rpl) = mean(testedFeat.chrom_variance, 'omitnan');
                Summary.Std_ChrV(Id2Rpl) = std(testedFeat.chrom_variance, [], 'omitnan');
                Summary.chr_Assym(Id2Rpl) = mean(testedFeat.chrom_assym, 'omitnan');
                Summary.mean_Accu(Id2Rpl) =  mean(testedFeat.Accu_mass, 'omitnan');
                Summary.Std_Accu(Id2Rpl) =  std(testedFeat.Accu_mass, [], 'omitnan');
                Summary.mean_MZV(Id2Rpl) =  mean(testedFeat.ms_variance, 'omitnan');
                Summary.Std_MZV(Id2Rpl) = std(testedFeat.ms_variance, [], 'omitnan');
                Summary.MS_Assym(Id2Rpl) = mean(testedFeat.ms_assym, 'omitnan');
                Summary.n(Id2Rpl) = height(testedFeat);
                Summary.NotDone(Id2Rpl) = false;
            end
        end
    end
end

Summary.NotDone = true(height(Summary), 1);
while 1

    myOutliers = isoutlier(Summary.Std_Apex);
    [~, IdFeature] = intersect(Features.ID, Summary.ID(find(myOutliers & Summary.NotDone)));
    if isempty(IdFeature), break, end

    for ii = 1:numel(IdFeature)
        testedFeat = Features.Values{IdFeature(ii)};
        if height(testedFeat) == 2

            Id2Del = find(Summary.ID == IdFeature(ii));
            Summary(Id2Del, :) = [];
            Area.Values(IdFeature(ii), :) = NaN;
            Features.Values{IdFeature(ii)} = table();

        else
            XY = [testedFeat.Accu_mass, testedFeat.chrom_intensity_apex];
            gmd = gmdistribution(mean(XY), cov(XY));
            testedFeat(mahal(gmd, XY) > options.Merge.MahalThres, :) = [];

            if height(testedFeat) == 2

                Id2Del = find(Summary.ID == IdFeature(ii));
                Summary(Id2Del, :) = [];
                Area.Values(IdFeature(ii), :) = NaN;
                Features.Values{IdFeature(ii)} = table();

            else

                Id2Rpl = find(Summary.ID == IdFeature(ii));
                V = (1:sum(strcmp(obj.Summary.FinneeType, Tag2Fin)))';
                [~, ia, ib] = intersect(V, testedFeat.Id2Tgt);
                V(ia, 2:4) = [testedFeat.chrom_area(ib), testedFeat.chrom_centroid(ib), testedFeat.Accu_mass(ib)];
                Area.Values(IdFeature(ii), :) = V(:, 2)';
                Features.Values{IdFeature(ii)} = testedFeat;
                Summary.mean_IApex(Id2Rpl) = mean(testedFeat.chrom_intensity_apex, 'omitnan');
                Summary.Std_IApex(Id2Rpl)  = std(testedFeat.chrom_intensity_apex, [], 'omitnan');
                Summary.mean_Apex(Id2Rpl) = mean(testedFeat.chrom_apex, 'omitnan');
                Summary.Std_Apex(Id2Rpl)  = std(testedFeat.chrom_apex, [], 'omitnan');
                Summary.mean_ctr(Id2Rpl) = mean(testedFeat.chrom_centroid, 'omitnan');
                Summary.Std_ctr(Id2Rpl) = std(testedFeat.chrom_centroid, [], 'omitnan');
                Summary.mean_ChrV(Id2Rpl) = mean(testedFeat.chrom_variance, 'omitnan');
                Summary.Std_ChrV(Id2Rpl) = std(testedFeat.chrom_variance, [], 'omitnan');
                Summary.chr_Assym(Id2Rpl) = mean(testedFeat.chrom_assym, 'omitnan');
                Summary.mean_Accu(Id2Rpl) =  mean(testedFeat.Accu_mass, 'omitnan');
                Summary.Std_Accu(Id2Rpl) =  std(testedFeat.Accu_mass, [], 'omitnan');
                Summary.mean_MZV(Id2Rpl) =  mean(testedFeat.ms_variance, 'omitnan');
                Summary.Std_MZV(Id2Rpl) = std(testedFeat.ms_variance, [], 'omitnan');
                Summary.MS_Assym(Id2Rpl) = mean(testedFeat.ms_assym, 'omitnan');
                Summary.n(Id2Rpl) = height(testedFeat);
                Summary.NotDone(Id2Rpl) = false;
            end
        end
    end
end

Summary.NotDone = [];
%% CHECK FOR DOUBLONS
Summary.toTest = true(height(Summary), 1);

while 1
    doStop = true;
    StdROIs = [Summary.mean_Apex Summary.Std_Apex  Summary.mean_Accu Summary.Std_Accu Summary.toTest];

    for ii = 1:height(StdROIs)
        RS = abs((StdROIs(:, 1) - StdROIs(ii, 1))./(2*((StdROIs(:, 2) + StdROIs(ii, 2)))));
        RS(:, 2) = abs((StdROIs(:, 3) - StdROIs(ii, 3))./(2*((StdROIs(:, 4) + StdROIs(ii, 4)))));
        RS(:, 3) = sqrt(RS(:, 1).^2 + RS(:, 2).^2);

        IdDoublons = find(RS(:, 3) < options.Doublons.CritRes(1));
        if numel(IdDoublons) > 1 & all(Summary.toTest(IdDoublons))

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

            ToKtoD = IdDoublons;
            sTarget = target(Summary.ID(IdDoublons), :);
            ToKtoD(:, 2) = NaN;
            ToKtoD(:, 3) = Summary.ID(IdDoublons);
            for jj = 1:10
                NNodes{jj} = table(); 

            end
            
            for jj = 1:height(NAN)

                [~, idX] = min(...
                    ((NAN.Accu_mass(jj) - sTarget.Targetmz)).^2+...
                    ((NAN.chrom_apex(jj) - sTarget.TargetTime)).^2);
                NNodes{idX} = [NNodes{idX}; NAN(jj, :)];
            end

            for jj = 1:numel(NNodes)
                pAN = NNodes{jj};

                try

                    if isempty(pAN), continue; end

                    while numel(pAN.Id2Tgt) ~= numel(unique(pAN.Id2Tgt))
                        [~, id] = max(((pAN.Accu_mass - target.Targetmz(ii))/std(pAN.Accu_mass)).^2 +...
                            ((pAN.chrom_apex - target.TargetTime(ii))/std(pAN.chrom_apex)).^2);
                        pAN(id, :) = [];

                    end
                catch
                    disp("pp")

                end

                if height(pAN) < options.Merge.minRepl(1), continue; end
                XY = [pAN.chrom_apex, pAN.Accu_mass];
               
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
                        Summary.toTest(IdTgt) = false;
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
                            Summary.toTest(IdTgt) = false;

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
Summary.toTest = [];

listName = tempname(fullfile(obj.Path2Project, 'myProjectPeakLists'));
myMrgPeakList = {};
myMrgPeakList.Name = listName;
myMrgPeakList.PeakList = PeakList;
myMrgPeakList.TargetedAnalysis{1}.Date = datetime('now');
myMrgPeakList.TargetedAnalysis{1}.options = options;
myMrgPeakList.TargetedAnalysis{1}.Target = target;
myMrgPeakList.TargetedAnalysis{1}.Summary =  Summary;
myMrgPeakList.TargetedAnalysis{1}.Area = Area;
myMrgPeakList.TargetedAnalysis{1}.Features = Features;
myMrgPeakList.TargetedAnalysis{1}.mapFiles = myMapFiles;

obj.FeaturesLists{end+1}.Original.Name = myMrgPeakList.Name;
save(myMrgPeakList.Name, 'myMrgPeakList')
myProject = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Project, 'myProject.mat'), 'myProject')

function OutputList = MyRecursive(inputList, Dmz, tm, minsize, OutputList)
        ThisTag = true;

        inputList = sortrows(inputList, "Accu_mass");
        MySplits = ((inputList.Accu_mass(2:end)-inputList.Accu_mass(1: end-1))./...
            ((inputList.Accu_mass(1: end-1)+ inputList.Accu_mass(2:end))/2))*1000000;
        MySplits(:, 2) = MySplits(:, 1) > Dmz;
        if any(MySplits(:, 2)), ThisTag = false; end

        list1 = unique([0; find(MySplits(:, 2)); height(inputList)]);
        for iMR = 1:numel(list1)-1
            SplitList1 = inputList(list1(iMR)+1:list1(iMR+1), :);

            if height(SplitList1) >= minsize
                SplitList1 = sortrows(SplitList1, "chrom_apex");
                My2nSplits = (SplitList1.chrom_apex(2:end) - SplitList1.chrom_apex(1:end-1));
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


