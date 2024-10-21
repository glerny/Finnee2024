function [myPeakList, myFinnee, options] = mkPeakList(myROIs, myFinnee, method, options)

if nargin == 2
    method = 'analysis_1';
    options = [];
elseif nargin == 3
    options = [];
end

if isempty(options)

    switch method
        case 'analysis_1'
            options.TimeInterval = [0 inf];
            options.doFilterRois = true;

            options.Analysis1.minMZ           = 3;
            options.Analysis1.minTime         = 5;
            options.Analysis1.WindowMZ        = 2;
            options.Analysis1.WindowTime      = 2;
            options.Analysis1.Min_nnz         = 10;
            options.Analysis1.sgfd            = 'average';
            options.Analysis1.StepSize        = 1.1;
            options.Analysis1.CritRes         = 0.8;

    end
end

AdditionalInformation = struct();
tStart = tic;
Id2Dataset = myROIs.Dataset;

%% Filter ROIs
if isfield(myROIs, 'Filter') & options.doFilterRois
    DesVar = myROIs.DescriptiveVars;
    Id2ROI = DesVar.IdROI;
    X = [DesVar.nnz, DesVar.CV, DesVar.Kurtosis, ...
        DesVar.ROIQual1, DesVar.ROIQual2, DesVar.LocalMaxDens, DesVar.LocalMax];
    Id2Rem = find(isnan(any(X, 2)));
    Id2ROI(Id2Rem) = [];
    X(Id2Rem, :) = [];
    FilteredData = [Id2ROI, str2double(myROIs.Filter.model.predict(X))];

    ROI2Keep = FilteredData(FilteredData(:, 2) == 1, 1);
    [~,ia] = intersect(myROIs.Targets.ID, ROI2Keep);
    myTargets =  myROIs.Targets(ia, :);
    ROIDescr  =  myROIs.DescriptiveVars(ia, :);

else
    myTargets =  myROIs.Targets;
    ROIDescr  =  myROIs.DescriptiveVars;

end

myPeakList = [];
for ii = 1:size(myTargets, 1)
    if myTargets.TimeMax(ii) < options.TimeInterval(1) ...
            | myTargets.TimeMin(ii) > options.TimeInterval(2)
        continue

    end

    if all(myTargets.sizeROI(ii, :) == [0 0]), continue; end

    myROI = ROI(myROIs, myTargets.ID(ii));
    myAnalysis = myROI.analysis_1(options.Analysis1);

    if ~isempty(myAnalysis)
        myAnalysis.ID2ROI = ones(size(myAnalysis.ID2ROI))*myTargets.ID(ii);
        myPeakList = [myPeakList; myAnalysis];
    end
end
myPeakList.ID = (1:size(myPeakList, 1))';


save(fullfile(myROIs.myPath, 'myPeakList.mat'), 'myPeakList')

AdditionalInformation.nonEndingErrors{1} = '';
AdditionalInformation.ComputingTime = toc(tStart);
AdditionalInformation.CreatedObject = myROIs;
newLineTable.TypeOfAction = 'mkPeakTable';
newLineTable.Options4Creations = {options};
newLineTable.DateOfCreation = datetime;
newLineTable.Output = fullfile(myROIs.myPath, 'myROIs.mat');
newLineTable.AdditionalInformation = AdditionalInformation;
newLineTable = struct2table(newLineTable, 'AsArray', true);

if isempty(myFinnee.Datasets.SecondaryActions{Id2Dataset})
    myFinnee.Datasets.SecondaryActions{Id2Dataset} = newLineTable;

else
    myFinnee.Datasets.SecondaryActions{Id2Dataset} =  [myFinnee.Datasets.SecondaryActions{Id2Dataset} ...
        ; newLineTable];
end
save(fullfile(myFinnee.Path2Fin, 'myFinnee.mat'), 'myFinnee');





