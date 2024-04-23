function myMrgPeakList = doQuantitative(obj, Tag2Fin, Id2PeakList, replace, timeAlign, options)

%% INTRODUCTION
if isempty(options)

    options.ForROIs.SigMZ       = 3;
    options.ForROIs.SigTm       = 12;
    options.ForROIs.TgtDataset  = 'Dataset2';
    options.ForROIs.IdDataset   = 2;
    options.ForROIs.Folder      = 'ROI4Quant';
    options.ForROIs.alfa        = 0.99;

    options.Analysis.sigMZ        = 1;
    options.Analysis.minTime      = 5;
    options.Analysis.WindowTime   = 4;
    options.Analysis.WindowMZ     = 2;
    options.Analysis.sgfd         = 'average';
    options.Analysis.critRes      = 0.8;
    options.Analysis.DynRange     = 0.05;
    options.Analysis.Orders       = 1;

    options.Merge.maxPPM          = 3;
    options.Merge.maxTm           = 0.5;
    options.Merge.maxRs           = 10;
    options.Merge.S2N             = [3, 10];
    options.Merge.n2n             = [0.9, 0.75];
    options.Merge.minRepl         = [4, 3];

    options.Doublons.mzPpm        = 0.5;
    options.Doublons.tmMin        = 0.1;

end

%% INITIALISATION %%

load(obj.FeaturesLists{Id2PeakList}.Original.Name);
mySummary = myMrgPeakList.Summary;

%% END INITIALISATION %%


target = table();
target.ID = mySummary.idMF; 
target.TargetTime = mySummary.Mean_Chrom_PeakCenter;
target.StdDev_time = tinv(.99, mySummary.sizeMF-1).*mySummary.Std_Chrom_PeakCenter;
target.PeakStdDev_time = sqrt(mySummary.Mean_Chrom_PeakVariance);
target.TimeMin = mySummary.Mean_Chrom_PeakCenter ... % center of the peak
    -3*mySummary.Std_Chrom_PeakCenter...             % variation of the center of the peak
    -options.ForROIs.SigTm*sqrt(mySummary.Mean_Chrom_PeakVariance);
target.TimeMax = mySummary.Mean_Chrom_PeakCenter ...
    +3*mySummary.Std_Chrom_PeakCenter...                           
    +options.ForROIs.SigTm*sqrt(mySummary.Mean_Chrom_PeakVariance);
target.Targetmz = mySummary.Mean_MS_AccuMass2;
target.StdDev_MS = tinv(.99, mySummary.sizeMF-1).*mySummary.Std_MS_AccuMass2;
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
myMapFiles.OriFile = []; myMapFiles.FinneeFile = [];  myMapFiles.toDO = [];
myMapFiles.TgtFolder = [];  myMapFiles.tempname = [];

for ii = 1:size(obj.Summary, 1)
    if strcmp(obj.Summary.FinneeType{ii}, Tag2Fin)
        myTgtFolder = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], ...
            options.ForROIs.TgtDataset);
        myMapFiles.Id2ii(end+1, 1) = ii;
        myMapFiles.Id(end+1, 1) = numel(myMapFiles.Id2ii);
        myMapFiles.OriFile{end+1, 1} = obj.Summary.FileID{ii};
        myMapFiles.FinneeFile{end+1, 1} = fullfile(obj.Summary.FolderID{ii}, [obj.Summary.FileID{ii}, '.fin'], 'myFinnee.mat');

        if ~(exist(fullfile(myTgtFolder, options.ForROIs.Folder), "dir") == 7) 
             myMapFiles.toDO(end+1, 1) = true;

        elseif replace
             myMapFiles.toDO(end+1, 1) = true;

        else
             myMapFiles.toDO(end+1, 1) = false;

        end
        myMapFiles.TgtFolder{end+1, 1} = myTgtFolder;
        myMapFiles.tempname{end+1, 1} = tempname;

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

        if isempty(timeAlign)
            corTarget = target;

        else
            corTarget = target;
            for jj = 1:size(corTarget, 1)
                corTarget.TargetTime(jj) = reversepolyval(corTarget.TargetTime(jj), timeAlign{ii});
                corTarget.TimeMin(jj) = reversepolyval(corTarget.TimeMin(jj), timeAlign{ii});
                corTarget.TimeMax(jj) = reversepolyval(corTarget.TimeMax(jj), timeAlign{ii});

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
Centers = {}; Centers.ID = target.ID; Centers.Values = [];
AccMass = {}; AccMass.ID = target.ID; AccMass.Values = [];
Summary = {};  Summary.ID = []; 
Summary.mean_IApex = []; Summary.Std_IApex = [];
Summary.mean_Apex = []; Summary.Std_Apex = [];
Summary.mean_ctr = []; Summary.Std_ctr = [];
Summary.mean_ChrV = []; Summary.Std_ChrV = [];
Summary.mean_Accu = []; Summary.Std_Accu = [];
Summary.mean_MZV = []; Summary.Std_MZV = []; 
Summary.chr_Assym = []; Summary.MS_Assym = [];
Summary.n2n05 = [];

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
        AN = [AN; myROIAnalysis];
    end

    % FILTER DATA
    if all(isnan(AN.Noise))
        AN(AN.n2n > options.Merge.n2n(1), :) = [];
    
    else
        Noise = mean(AN.Noise, 'omitnan');
        AN(AN.chrom_intensity_apex < options.Merge.S2N(1)*Noise, :) = [];
    
    end

    if size(AN, 1) <  options.Merge.minRepl(1)
        continue

    end

    %% CLUSTER DATA

    Nodes = {}; SNodes = [];
    countMe = 1;
    while ~isempty(AN)
        AN = sortrows(AN,"chrom_intensity_apex","descend");
        RS = (AN.chrom_centroid(1) - AN.chrom_centroid)./...
            (2*(sqrt(AN.chrom_variance(1)) + sqrt(AN.chrom_variance)));
        RS(:, 2) = (AN.Accu_mass(1) - AN.Accu_mass)./...
            (2*(sqrt(AN.ms_variance(1)) + sqrt(AN.ms_variance)));
        RS(:, 3) = sqrt(RS(:, 2).^2 + RS(:, 1).^2);
        [~, idX] = sort(RS(:, 3));
        AN = AN(idX, :);
        RS = RS(idX, :);

        mList = 1;
        LDst = AN.Id2Tgt(1);
        for jj = 2:numel(idX)

            if RS(jj, 3) > options.Merge.maxRs
                break

            end

            if abs((mean(AN.Accu_mass(mList)) - AN.Accu_mass(jj)))/mean(AN.Accu_mass(mList))*1000000 < options.Merge.maxPPM ...
                    & abs(mean(AN.chrom_centroid(mList)) - AN.chrom_centroid(jj)) < options.Merge.maxTm
                if any(LDst == AN.Id2Tgt(jj))
                    break

                else
                    mList(end+1) = jj;
                    LDst(end + 1) = AN.Id2Tgt(jj);

                end
            end

        end

        while 1
            stopMe = true;

            if size(mList, 2) > 1
                dList = [];
                for jj = 1:numel(mList)

                    RS = (AN.chrom_centroid(mList(jj)) - AN.chrom_centroid)./...
                        (2*(sqrt(AN.chrom_variance(mList(jj))) + sqrt(AN.chrom_variance)));
                    RS(:, 2) = (AN.Accu_mass(mList(jj)) - AN.Accu_mass)./...
                        (2*(sqrt(AN.ms_variance(mList(jj))) + sqrt(AN.ms_variance)));
                    RS(:, 3) = sqrt(RS(:, 2).^2 + RS(:, 1).^2);


                    IdA = find(RS(:, 3) == min(RS(RS(:, 3) ~= 0, 3)));
                    if ~any(mList == IdA)
                        dList(jj) = 1;
                        stopMe = false;

                    else
                        dList(jj) = 0;

                    end
                end
                mList = mList(dList == 0);
            end

            if stopMe, break; end
        end

        if numel(mList) > options.Merge.minRepl
            ia = isoutlier(AN.Accu_mass(mList), "ThresholdFactor", 5);
            mList(ia) = [];

        end

        if isempty(mList), mList = 1; end

        if numel(mList) > options.Merge.minRepl(1)
            fPeaks = AN(mList, :);
            Nodes{countMe} = fPeaks;

            SNodes(countMe, :) = [countMe, size(fPeaks, 1), mean(fPeaks.chrom_area), ...
                mean(fPeaks.chrom_apex), std(fPeaks.chrom_apex),...
                mean(fPeaks.chrom_centroid), std(fPeaks.chrom_centroid),...
                mean(fPeaks.chrom_variance), std(fPeaks.chrom_variance), ...
                mean(fPeaks.Accu_mass), std(fPeaks.Accu_mass), ...
                mean(fPeaks.ms_variance), std(fPeaks.ms_variance)];
            countMe = countMe+1;

        end
        AN(mList, :) = [];
    end
   % END CLUSTER

   if ~isempty(SNodes)
       [~, idX] = min(...
           ((SNodes(:, 10)-target.Targetmz(ii))./SNodes(:, 11)).^2+...
           ((SNodes(:, 4)-target.TargetTime(ii))./SNodes(:, 5)).^2);
       cAN = Nodes{idX};
       V = (1:numel(obj.ListOfFiles))';
       [~, ia, ib] = intersect(V, cAN.Id2Tgt);
       V(ia, 2:4) = [cAN.chrom_area(ib), cAN.chrom_centroid(ib), cAN.Accu_mass(ib)];

       if all(isnan(cAN.Noise))
           IdM = cAN.n2n < options.Merge.n2n(2);
           n2n05 = sum(IdM);
           if sum(IdM) >= options.Merge.minRepl(2)
               IdROI = target.ID(ii);
               Area.Values(IdROI, :) = V(:, 2)';
               Centers.Values(IdROI, :) = V(:, 3)';
               AccMass.Values(IdROI, :) = V(:, 4)';
               Summary.ID(end+1, 1) = IdROI; 
               Summary.mean_IApex(end+1, 1) = mean(cAN.chrom_intensity_apex, 'omitnan');
               Summary.Std_IApex(end+1, 1)  = std(cAN.chrom_intensity_apex, [], 'omitnan');
               Summary.mean_Apex(end+1, 1) = mean(cAN.chrom_apex(IdM), 'omitnan');
               Summary.Std_Apex(end+1, 1)  = std(cAN.chrom_apex(IdM), [], 'omitnan');
               Summary.mean_ctr(end+1, 1) = mean(cAN.chrom_centroid(IdM), 'omitnan');
               Summary.Std_ctr(end+1, 1) = std(cAN.chrom_centroid(IdM), [], 'omitnan');
               Summary.mean_ChrV(end+1, 1) = mean(cAN.chrom_variance(IdM), 'omitnan');
               Summary.Std_ChrV(end+1, 1) = std(cAN.chrom_variance(IdM), [], 'omitnan');
               Summary.chr_Assym(end+1, 1) = mean(cAN.chrom_assym(IdM), 'omitnan');
               Summary.mean_Accu(end+1, 1) =  mean(cAN.Accu_mass(IdM), 'omitnan');
               Summary.Std_Accu(end+1, 1) =  std(cAN.Accu_mass(IdM), [], 'omitnan');
               Summary.mean_MZV(end+1, 1) =  mean(cAN.ms_variance(IdM), 'omitnan');
               Summary.Std_MZV(end+1, 1) = std(cAN.ms_variance(IdM), [], 'omitnan');
               Summary.MS_Assym(end+1, 1) = mean(cAN.ms_assym(IdM), 'omitnan');
               Summary.n2n05(end+1, 1) = n2n05;
           end

       else
           Noise = mean(cAN.Noise, 'omitnan');
           IdM = cAN.chrom_intensity_apex >= options.Merge.S2N(2)*Noise;
           n2n05 = sum(IdM);
           if sum(IdM) >= options.Merge.minRepl(2)
               IdROI = target.ID(ii);
               Area.Values(IdROI, :) = V(:, 2)';
               Centers.Values(IdROI, :) = V(:, 3)';
               AccMass.Values(IdROI, :) = V(:, 4)';
               Summary.ID(end+1, 1) = IdROI; 
               Summary.mean_IApex(end+1, 1) = mean(cAN.chrom_intensity_apex, 'omitnan');
               Summary.Std_IApex(end+1, 1)  = std(cAN.chrom_intensity_apex, [], 'omitnan');
               Summary.mean_Apex(end+1, 1) = mean(cAN.chrom_apex(IdM), 'omitnan');
               Summary.Std_Apex(end+1, 1)  = std(cAN.chrom_apex(IdM), [], 'omitnan');
               Summary.mean_ctr(end+1, 1) = mean(cAN.chrom_centroid(IdM), 'omitnan');
               Summary.Std_ctr(end+1, 1) = std(cAN.chrom_centroid(IdM), [], 'omitnan');
               Summary.mean_ChrV(end+1, 1) = mean(cAN.chrom_variance(IdM), 'omitnan');
               Summary.Std_ChrV(end+1, 1) = std(cAN.chrom_variance(IdM), [], 'omitnan');
               Summary.chr_Assym(end+1, 1) = mean(cAN.chrom_assym(IdM), 'omitnan');
               Summary.mean_Accu(end+1, 1) =  mean(cAN.Accu_mass(IdM), 'omitnan');
               Summary.Std_Accu(end+1, 1) =  std(cAN.Accu_mass(IdM), [], 'omitnan');
               Summary.mean_MZV(end+1, 1) =  mean(cAN.ms_variance(IdM), 'omitnan');
               Summary.Std_MZV(end+1, 1) = std(cAN.ms_variance(IdM), [], 'omitnan');
               Summary.MS_Assym(end+1, 1) = mean(cAN.ms_assym(IdM), 'omitnan');
               Summary.n2n05(end+1, 1) = n2n05;
           end

       end

   end
end

%% CHECK FOR DOUBLON
In2Rem = [];
Summary =  struct2table(Summary);
for ii = 1:size(Summary, 1)
    I2T = find(abs(Summary.mean_ctr(ii) - Summary.mean_ctr) <=  options.Doublons.tmMin ... 
        & abs(Summary.mean_Accu(ii) - Summary.mean_Accu)./Summary.mean_Accu(ii)*1000000 <= options.Doublons.mzPpm);

    if numel(I2T) > 1
        I2I2T = find(Summary.n2n05(I2T) == max(Summary.n2n05(I2T)));
        new2Keep = I2T(I2I2T);
        [~, I2K] = max(Summary.mean_IApex(new2Keep));
        new2Keep = new2Keep(I2K);
        
        In2Rem = unique([In2Rem; I2T(I2T ~= new2Keep)]);
    end
end
Summary(In2Rem, :) = [];

if ~isfield(myMrgPeakList, 'TargetedAnalysis')
    myMrgPeakList.TargetedAnalysis = {};
end
myMrgPeakList.TargetedAnalysis{end+1}.Date = datetime('now');
myMrgPeakList.TargetedAnalysis{end}.options = options;
myMrgPeakList.TargetedAnalysis{end}.Target = target;
myMrgPeakList.TargetedAnalysis{end}.Summary =  Summary;
myMrgPeakList.TargetedAnalysis{end}.Area = Area;
myMrgPeakList.TargetedAnalysis{end}.Centers = Centers;
myMrgPeakList.TargetedAnalysis{end}.AccMass = AccMass;
myMrgPeakList.TargetedAnalysis{end}.mapFiles = myMapFiles;
save(obj.FeaturesLists{Id2PeakList}.Original.Name, 'myMrgPeakList')


