%% DESCRIPTION 
%% DESCRIPTION 
% The BASELINECORRECTION method is used to correct a full profile dataset 
% from baseline drift. BASELINECORRECTION used a graphic user interface 
% (GUI)to select the m/z values that are most relevant for correcting and 
% to select and optimize the baseline correction methods.
% BASELINECORRECTION will create a new dataset. BASELINECORRECTION can be
% used with the original profile dataset. However better results will be
% obtained if the dataset has already been corrected from background noise.
% For more information see 
% <https://github.com/glerny/Finnee2016/wiki/Basic-operation-with-Dataset>
%
%% INPUT PARAMETERS
% *Compulsory 
%   _obj_, the Finnee object 
%   _dts_, the indice to the dts to correct (i.e. myFinnee.Datasets{dts}
% * Optional 
%   _ spikes_ (see the method @Finnee\FilterDataset) for additional 
%       information. Remove spikes in every MS scan. If used, the options 
%       should be followed by 0,1, 2 ot 3 to define the size of spikes to
%       be considered for removal. The default value is 2; 0 allows to turn
%       off spikes removal.
%
%% EXAMPLES
% * myFinnee = myFinnee.BaselineCorrection(2); Will correct dataset 2 for Path2cFin
%   baseline drift. A new dataset will be created. 
%
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = BaselineCorrection(obj, dts, varargin)
    
%% CORE OF THE FUNCTION
% 1- Initialisation and options
% TODO: check database

if nargin == 2 & isstruct(dts)
    options = dts;
    dts = options.parentDataset;
    options.ParallelMe = true;

else
    options = checkVarargin(dts, varargin{:});
    
end

% 2- Selection of m/z values
Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name);
infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;

if isempty(infoDataset.minNoise) | ~isfinite(infoDataset.minNoise)
    infoDataset.minNoise = 1;
end

minNoise = infoDataset.minNoise;

fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
Profiles = fread(fidRead, [infoDataset.Profiles.size], "double");
fclose(fidRead);
Profiles = array2table(Profiles);
ColumnNames = {}; ColumnUnits = {};
for ii = 1:numel(infoDataset.Profiles.column)
        ColumnNames{ii} = infoDataset.Profiles.column{ii}.name;
        ColumnLabels{ii} = infoDataset.Profiles.column{ii}.labelUnits;
        ColumnUnits{ii} = infoDataset.Profiles.column{ii}.units;
        ColumnFO{ii}    = infoDataset.Profiles.column{ii}.formatUnits;
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

FIS = [Spectra.("m/z axis"), Spectra.Max_Continuous_Profile];
bool = false;
if isempty(options.size4BslCor)
    bool = true;
    f1 = figure;
    [N, edges] = histcounts(FIS(:,2));
    bar(edges, [0, N]);
    axis([0 inf 0 inf])
    hPlot = gcf;
    title('Frequency profile');
    h = msgbox('Select the minumum continuous profile length in a SmzP for it to be selected for baseline correction', 'Correct','modal');
    uiwait(h)
    figure(hPlot);
    [fmax, ~] = ginput(1);
    if isempty(fmax)
        return
    end
    fmax = inputdlg('Confirm the value', 'Correct', 1, {num2str(fmax)});
    options.size4BslCor = str2double(fmax{1});
    try close(f1), catch, end
end
IndBsl = FIS(:, 2) >= options.size4BslCor;

if isempty(options.size4Noise)
    fmax = inputdlg('Enter the maximum continuous profile length over noise in a SmzP for it to be selected as noise', 'OK', 1, "2");
    options.size4Noise = str2double(fmax{1});
    try close(f1), catch, end
end
IndNoise = FIS(:, 2) > options.size4Noise;
IndNoise = smoothdata(IndNoise, 'movmean', 10);
IndNoise = IndNoise == 0;

if ~options.ParallelMe, h = waitbar(0, 'Generating matrix of profiles'); end
profile = Profiles.Scan_Time;
matRes = zeros(sum(IndBsl), numel(Profiles.Scan_Time));
profile(:, 7) = 0;

% 3- Creating the matrix of profiles to be corrected
for ii = 1:length(profile(:,1))
    if ~options.ParallelMe, waitbar(ii/length(profile(:,1))); end
    ScanName = fullfile(obj.Path2Fin, ['Dataset', num2str(dts)],...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);
    
    if isempty(myMS)
        MS = Spectra.("m/z axis"); MS(:, 2) =0;

    else
        [~, Id1] = intersect(Spectra.("m/z axis"), myMS(:,1));
        if numel(Id1) < numel(myMS(:,1)), error(""); end
        MS = Spectra.("m/z axis"); MS(Id1, 2) = myMS(:,2);
        
    end
    
    matRes(:, ii) = MS(IndBsl, 2); %#ok<*AGROW>
    profile(ii,2) = sum( MS(IndBsl, 2));
    profile(ii,3) = max(MS(IndBsl, 2));
    profile(ii,4) = sum( MS(IndNoise, 2));

    try
        profile(ii,5) = max(MS(IndNoise, 2));
    catch
        profile(ii,5) = 0;
    end

    profile(ii,6) = sum( MS(~IndBsl & ~IndNoise, 2));
    profile(ii,7) = max( MS(~IndBsl & ~IndNoise, 2));
end

try close(h), catch, end

if bool
    f2 = figure('Name', 'Selected Profiles for Baseline Corrections');
    subplot(3, 1, 1)
    plot(profile(:,1), profile(:,3));
    title(['Selected single m/z profiles for correction (max plot): ', num2str(sum(IndBsl))])
    ylabel('Time /min'); %TODO: check units

    subplot(3, 1, 2)
    plot(profile(:,1), profile(:,5));
    title(['Single m/z profiles set for noise removal (mean plot): ', num2str(sum(IndNoise))])
    xlabel('Time /min');
    ylabel('Absorbance /Arb. units');

    subplot(3, 1, 3)
    plot(profile(:,1), profile(:,7));
    title(['Single m/z profiles that will not be changed (max plot): ', num2str(sum(~IndNoise & ~IndBsl))])
    xlabel('Time /min');
    ylabel('Absorbance /Arb. units');

    uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    button = questdlg('Are you happy with those conditions?');
    if ~strcmp(button, 'Yes')
        return
    end

    try close(f2), catch, end
end

% 5- Select the baseline method
if isempty(options.baseline)
     warning('off') %#ok<*WNOFF>
     [myBaselineFunction, ExitFlag] = gui4basCor(profile(:,1), MS(IndBsl, 1), matRes, options.XLim, options.HRMS, minNoise, options.size4BslCor);
     if ExitFlag == -1
         return
     end
     warning('on') %#ok<*WNON>
     
else
    Oo = strsplit(options.baseline, ":");
    myBaselineFunction.Algo = Oo{1};
    if numel(Oo) >= 2; myBaselineFunction.pm1 = str2double(Oo{2}); end
    if numel(Oo) >= 3; myBaselineFunction.pm2 = str2double(Oo{3}); end
    
end
AdditionalInformation = struct();
AdditionalInformation.BaselineCorrectionParameters.options = options;
AdditionalInformation.BaselineCorrectionParameters.XLim = options.XLim;
AdditionalInformation.BaselineCorrectionParameters.BaselineFunction = myBaselineFunction;
tStart = tic;

% 6- corecting the profiles
corProf = zeros(size(matRes));
if ~options.ParallelMe, h = waitbar(0,'Correcting profiles'); end

Noise = [];
myCount = 1;
IdN = find(IndBsl);

for ii = 1:size(corProf, 1)

   if ~options.ParallelMe, waitbar(ii/size(corProf, 1)); end
    
    prf      = profile(:,1);
    prf(:,2) = matRes(ii, :);
    isNan = isnan(prf(:,2));
    z = zeros(size(prf(:,1)));
    bslPts = false(size(prf(:,1)));
    
    switch myBaselineFunction.Algo
        case 'PolyFit'
            n = myBaselineFunction.pm1;

            spltprf = splitProfiles(prf, options.size4BslCor);
            for jj = 1:numel(spltprf)

                if size(spltprf{jj}, 1) > options.HRMS.shiftbsl
                    [Z, BSLPTS] = doPF(spltprf{jj}, n);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;
                else
                    [Z, BSLPTS] = doPF_HRMS(spltprf{jj}, options.HRMS.k, options.HRMS.n);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;

                end
            end

        case 'ArPLS'
            lambda = myBaselineFunction.pm1;
            ratio = myBaselineFunction.pm2;

            spltprf = splitProfiles(prf, options.size4BslCor);
            for jj = 1:numel(spltprf)

                if size(spltprf{jj}, 1) > options.HRMS.shiftbsl
                    [Z, BSLPTS] = doArPLS(spltprf{jj}(:, 2), lambda, ratio);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;
                else
                    [Z, BSLPTS] = doPF_HRMS(spltprf{jj}, options.HRMS.k, options.HRMS.n);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;

                end
            end

        case 'ArPLS2'
            lambda = myBaselineFunction.pm1;

            spltprf = splitProfiles(prf, options.size4BslCor);
            for jj = 1:numel(spltprf)

                if size(spltprf{jj}, 1) > options.HRMS.shiftbsl
                    [Z, BSLPTS] = doArPLS2(spltprf{jj}(:, 2), lambda);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;
                else
                    [Z, BSLPTS] = doPF_HRMS(spltprf{jj}, options.HRMS.k, options.HRMS.n);
                    [~, Ida] = intersect(prf(:,1), spltprf{jj}(:, 1));
                    z(Ida) = Z;
                    bslPts(Ida) = BSLPTS;
                    Noise(myCount, 1:5) = ...
                        [Spectra.("m/z axis")(IdN(ii)),...
                        spltprf{jj}(1, 1), spltprf{jj}(end, 1),...
                        max(infoDataset.minNoise, 2*std(spltprf{jj}(BSLPTS,2)-  Z(BSLPTS)), 'omitnan'),...
                        mean(Z(BSLPTS), 'omitnan')];
                    myCount = myCount +1;

                end
            end
    end

    Yc           = prf(:,2) - z;
    Yc(isnan(Yc) | Yc < 0) = 0;
    corProf(ii, :) = Yc;
    
end
corProf(corProf < 0) = 0;
corProf(~isfinite(corProf)) = 0;
try close(h), catch, end


%% 7. Recording new dataset and saving
% Create new dat file
ii = 1;
while 1
    newDataset = ['Dataset', num2str(ii)];
    if ~any(strcmp(newDataset, obj.Datasets.Name)), break; end
    ii = ii + 1;
end
mzInterval = [inf 0];
mkdir(fullfile(obj.Path2Fin, newDataset));
mkdir(fullfile(obj.Path2Fin, newDataset, 'Scans'));
newProfiles = [];
newSpectra = Spectra.("m/z axis"); newSpectra(:, 7) = 0;

%TIP, BPP and more
if ~options.ParallelMe, h = waitbar(0,'processing scans'); end
CurSeq = zeros(numel(newSpectra(:,1)), 2);
for ii = 1:numel(Profiles.Scan_Index)
    if ~options.ParallelMe,waitbar(ii/numel(Profiles.Scan_Index)); end

    ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);
    
    if isempty(myMS)
         MS = Spectra.("m/z axis"); MS(:, 2) = 0;
        
    else
        [~, Id1] = intersect(Spectra.("m/z axis"), myMS(:,1));
        if numel(Id1) < numel(myMS(:,1)), error(""); end
        MS = Spectra.("m/z axis"); MS(Id1, 2) = myMS(:,2);
    
    end
    
    MS(IndBsl, 2) = corProf(:, ii); % baselie correction
    MS(IndNoise, 2) = 0;            % Noise correction
    
    % Filter spikes if needed
    if ~isempty(MS) && options.RemSpks
        spkSz = options.SpkSz;
        MS    = spikesRemoval(MS, spkSz );
    end

    %Record MZSpectra
    IdNnz = MS(:, 2) > 0;
    newSpectra(:,2) = newSpectra(:,2) + double(IdNnz);
    IdStart = CurSeq(:, 2) == 0 & IdNnz;

    CurSeq(IdStart, 1) = ii;
    CurSeq(IdNnz, 2) = CurSeq(IdNnz, 2) + 1;

    % Check Length
    IdMatch = newSpectra(:,3) < CurSeq(:, 2) & ~IdNnz;
    newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
    newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
    CurSeq(~IdNnz, 2) = 0;
    newSpectra(:,5) = newSpectra(:,5) + MS(:, 2);
    newSpectra(:,6) = max([newSpectra(:,6), MS(:, 2)], [], 2);

    % reduced trailing zero in excess
    if ~isempty(MS)
        provMat      = [MS(2:end, 2); 0];
        provMat(:,2) = MS(:, 2);
        provMat(:,3) = [0; MS(1:end-1, 2)];
        MS         = MS(sum(provMat, 2) > 0, :);
    end

    if min(MS(:,1)) < mzInterval(1), mzInterval(1) = min(MS(:,1)); end
    if max(MS(:,1)) > mzInterval(2), mzInterval(2) = max(MS(:,1)); end

    if isempty(MS) | size(MS, 1) == 1 | sum(MS(:, 2)) == 0
        newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
            Profiles.Scan_Time(ii), ...
            Profiles.Injection_Time(ii), 0, 0, 0, ...
            nnz(MS(:,2)), size(MS), nan];
        
    else
        
        M = ChrMoment(MS);
        newProfiles(ii, :) = [Profiles.Scan_Index(ii), ...
            Profiles.Scan_Time(ii), ...
            Profiles.Injection_Time(ii), max(MS(:,2)), sum(MS(:,2)), ...
            min(MS(MS(:, 2) ~= 0, 2)), nnz(MS(:,2)),...
            size(MS), M(2)];
    end

    fileName = fullfile(obj.Path2Fin, newDataset, 'Scans', ['Scan#', ...
        num2str(Profiles.Scan_Index(ii, 1)), '.dat']);
    [fidWriteDat, errmsg]  = fopen(fileName, 'wb');
    fwrite(fidWriteDat, MS(:), "double");
    fclose(fidWriteDat);
end
try close(h); catch, end
IdMatch = newSpectra(:,3) < CurSeq(:, 2);
newSpectra(IdMatch, 3) = CurSeq(IdMatch, 2);
newSpectra(IdMatch, 4) = CurSeq(IdMatch, 1);
newSpectra(:, 7) = newSpectra(:,5)./newSpectra(:,2);

provMat = [newSpectra(2:end, 5); 0];
provMat(:,2) = newSpectra(:, 5);
provMat(:,3) = [0; newSpectra(1:end-1, 2)];
newSpectra = newSpectra(sum(provMat, 2) > 0, :);
short_mzAxis =  newSpectra(:,1);

options.baseline = [myBaselineFunction.Algo, ':', num2str(myBaselineFunction.pm1), ':', num2str(myBaselineFunction.pm2)];
infoDataset.dateofcreation = datetime("now");
infoDataset.AdditionalInformation = AdditionalInformation;
infoDataset.AdditionalInformation.nonEndingErrors{1} = '';
infoDataset.AdditionalInformation.ComputingTime = toc(tStart);
infoDataset.Label = obj.Datasets.Labels{end};
infoDataset.mzInterval =  mzInterval;
infoDataset.Options = options;
infoDataset.Profiles.size = size(newProfiles);
infoDataset.Spectra.size = size(newSpectra);
infoDataset.NoiseVector = Noise;

if infoDataset.minNoise == 1
    infoDataset.minNoise = mean(Noise(:, 5) < prctile(Noise(:, 5), 5), 'omitnan');
end

obj.Datasets = [obj.Datasets; ...
    {newDataset, '', true, true, true, true, false, ...
    ['Dataset', num2str(dts)], obj.Datasets.Labels{Id2Dataset}, ...
    options,  'Baseline_correction', {}}];

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
    function options = checkVarargin(dts, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
                
        
        % 2- Default parameters
        options.RemSpks       = true;
        options.SpkSz         = 3;
        options.XLim          = [0 inf];
        options.size4BslCor   = []; 
        options.baseline      = {};
        options.parentDataset = dts;
        options.size4Noise    = [];
        options.wdw4corr        = 2;
        options.ParallelMe      = false;
        options.Spk4noise       = 3;
        options.Sig2Nois        = 3;
        options.HRMS.shiftbsl   = 50;
        options.HRMS.n          = 1;
        options.HRMS.k          = 10;



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
        
        tgtIx = input('XLim');
        if ~isempty(tgtIx)
            XmM = varargin{tgtIx +1};
            options.XLim(1) = min(XmM);
            options.XLim(2) = max(XmM);
            
        end
        
        tgtIx = input('Ln4BslCor');
        if ~isempty(tgtIx)
            options.size4BslCor = varargin{tgtIx +1};
        end
        
        tgtIx = input('Baseline');
        if ~isempty(tgtIx)
            options.baseline = varargin{tgtIx +1};
        end
        
        tgtIx = input('Ln4Noise');
        if ~isempty(tgtIx)
            options.size4Noise = varargin{tgtIx +1};
        end
    end
end
