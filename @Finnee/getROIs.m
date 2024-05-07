%% DESCRIPTION

%
%% INPUT PARAMETERS
%
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [myROIs, obj] = getROIs(obj, dts, options, method)

% METHOD will be used at a later time when multiple methods will be
% available

if nargin == 2
    options.method.name             = 'default';
    options.method.min_MZ_Len       = 5; % min_MZ_peak_Length    
    options.method.min_Chrom_Len    = 7; % min_Chrom_peak_Length
    options.method.sgfd             = 'medium';
    options.method.WindowMZ         = 2;
    options.method.WindowTime       = 3;
    options.method.min_Sig2minNoise = 3;
    options.method.maxSize          = [3000 5000];
    options.method.myThreshs        = [0, 0.5, 1, 10];



    options.filterROIs.do           = true;
    options.filterROIs.method       = 'Random Forest';
    options.filterROIs.minPerLabel  = 75;
    options.filterROIs.LabelData    = [];
    options.filterROIs.model        = {};
    options.filterROIs.FilteredData = [];
end

AdditionalInformation = struct();
tStart = tic;

%% INTRODUCTION : LOAD THE NEEDED DATA
% Match the entered target dataset (dts) to the correct indice
Id2Dataset = strcmp(['Dataset', num2str(dts)], obj.Datasets.Name); 

% Load the dataset data (info dataset, profiles, spectra and axis)
infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;

fidRead = fopen(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'Profiles.dat'), 'rb');
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

infoDataset = load(fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset}, 'infoDataset.mat'));
infoDataset = infoDataset.infoDataset;
mzAxis = Spectra.("m/z axis");

IntThres = options.method.min_Sig2minNoise*infoDataset.minNoise;

%% FIND THE LIMITS FOR EACH ROI RECURSIVELY
% First loop - Take the non_zeros elements spectra and split based on
% zeroes values
XY = [Spectra.("m/z axis"), Spectra.Non_Zeros_Elements];
myLim = find(XY(:,2) == 0);
myLim(1:end-1, 2) = myLim(2:end);
myLim(:,3) = myLim(:, 2) - myLim(:,1);
myLim(myLim(:, 3) <= 5, :) = [];
myLim(:,3) = 1;
myLim(:,4) = inf;
myLim(:,5) = 0;
Tgt = [];

for ii = 1:numel(Profiles(:,1))
    ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);
    fMS = Spectra.("m/z axis");

    nnz_XIP =zeros(size(myLim(:, 5)));

    if ~isempty(myMS)

        [~, Id1] = intersect(fMS, myMS(:,1));
        if numel(Id1) < numel(myMS(:,1)), error(""); end
        fMS(Id1, 2) = myMS(:,2);
    else

        fMS(:, 2) = 0;
    end

    Id2Add = find(myLim(:,3) <= ii & myLim(:,4) >= ii);

    for jj = 1:numel(Id2Add)
        JJ = Id2Add(jj);
        nnz_XIP(JJ) = nnz(fMS(myLim(JJ,1):myLim(JJ,2), 2));
    end

    IDC = nnz_XIP > 0;

    myLim(IDC & myLim(:, 5) == 0, 3) = ii-1;
    myLim(IDC & myLim(:, 5) == 0, 5) = 1;
    if any(~IDC & myLim(:, 5) > options.method.min_Chrom_Len)
        IDCUT = find(~IDC & myLim(:, 5) > options.method.min_Chrom_Len);
        data2add = myLim(IDCUT, :);
        data2add(:, 4) = ii;
        Tgt = [Tgt; data2add];
    end

    myLim(~IDC, 5) = 0;
    myLim(~IDC, 3) = ii;
    myLim(IDC, 5) = myLim(IDC, 5) + 1;
end
IDCUT = find(myLim(:,5) > options.method.min_Chrom_Len);
data2add = myLim(IDCUT, :);
data2add(:, 4) = ii;
data2add(:, 3) = data2add(:, 3)-1;
Tgt = [Tgt; data2add];
XScans{size(Tgt, 1)} = [];
XProfiles{size(Tgt, 1)} = [];

% Make Extraxted profiles XProfiles and Extracted scans Xscans from
% previous intervals
for ii = 1:numel(Profiles(:,1))
    ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
        'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
    fidRead = fopen(ScanName, 'rb');
    myMS = fread(fidRead, ...
        [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
        "double");
    fclose(fidRead);
    fMS = Spectra.("m/z axis");

    if ~isempty(myMS)

        [~, Id1] = intersect(fMS, myMS(:,1));
        if numel(Id1) < numel(myMS(:,1)), error(""); end
        fMS(Id1, 2) = myMS(:,2);
    else
        fMS(Id1, 2) = 0;
    end

    Id2Add = find(Tgt(:,3) <= ii & Tgt(:,4) >= ii);

    for jj = 1:numel(Id2Add)
        JJ = Id2Add(jj);
        if isempty(XScans{JJ})
            XScans{JJ} = fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
            XScans{JJ}(:,2) = 0 + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
            XProfiles{JJ}(1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
            XProfiles{JJ}(1, 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
        else
            XScans{JJ}(:,1) = XScans{JJ}(:,1) + fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
            XScans{JJ}(:,2) = XScans{JJ}(:,2) + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
            XProfiles{JJ}(end+1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
            XProfiles{JJ}(end  , 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));

        end
    end

end

%% REDO THE SPLIT AS LONG AS SMALLEST SPLITS ARE DETECTED
% Check XProfiles and Xsans, removed trailing and leading zeros and
% update intervals split and redo when needed
myList = table();
while 1
    if isempty(Tgt), break; end
    Lim2Del = false(numel(XScans), 1);
    for ii = 1:numel(XScans)

        % 1. remove leading and trailing zeros
        % 1.1 Spectra (m/z profile)
        xyS = [];
        xyS(:, 2:3) = XScans{ii};
        xyS(:, 1) = (Tgt(ii,1):Tgt(ii,2))';

        if xyS(1, 2) == 0
            is = find(xyS(:,2) > 0, 1, 'first') - 1;
        else
            is = 1;
        end

        if xyS(end, 2) == 0
            ie = find(xyS(:,2) > 0, 1, 'last') + 1;
        else
            ie = length(xyS(:,2));
        end

        xyS = xyS(is:ie, :); 
        XScans{ii} =  xyS(:, 2:3);
        Tgt(ii,1:2) = [xyS(1, 1) xyS(end, 1)];

        % 1.2 Profile
        xyP = [];

        xyP(:, 2:3) = XProfiles{ii};
        cAxis = (Tgt(ii,3):Tgt(ii,4))';
        cAxis(cAxis <= 0) = [];
        xyP(:, 1) = cAxis';

        if xyP(1, 2) == 0
            is = find(xyP(:,2) > 0, 1, 'first') - 1;
        else
            is = 1;
        end

        if xyP(end, 2) == 0
            ie = find(xyP(:,2) > 0, 1, 'last') + 1;
        else
            ie = length(xyP(:,2));
        end

        xyP = xyP(is:ie, :); 
        XProfiles{ii} =  xyP(:, 2:3);
        Tgt(ii, 3:4) = [xyP(1, 1) xyP(end, 1)];


        % 2. Check if valid via size
        if ~any(XScans{ii}(:, 2) >= options.method.min_Chrom_Len) ...
                | ~any(XProfiles{ii}(:, 2) >= options.method.min_MZ_Len)

            Lim2Del(ii) = true;
            continue
        end

        % 3. Check if scan is splitted by zeros
       
        if any(xyS(2:end-1, 2) == 0)
                Lim2Del(ii) = false;
                continue
        end

        % 4.2. FFind if zeros intensities split
        if any(xyP(2:end-1, 2) == 0)
                Lim2Del(ii) = false;
                continue
        end
        
        % If within this increment we arrive here, profile and scans are
        % correct and scan abd profile are not spli within our assumptions.
        % Data are recorded 
        myList = [myList; {Tgt(ii, :), {XScans{ii}}, {XProfiles{ii}}}];
        Lim2Del(ii) = true;
    end
    
    Tgt(Lim2Del, :) = [];
    XScans(Lim2Del) = [];
    XProfiles(Lim2Del) = [];

    newLim = [];
    TargetScans = [];

    for ii = 1:numel(XScans)
        
        splitedMe = false;
        % 1. Make profile & Spectra
        xyS = [];
        xyS(:, 2:3) = XScans{ii};
        xyS(:, 1) = (Tgt(ii,1):Tgt(ii,2))';

        xyP = [];
        xyP(:, 2:3) = XProfiles{ii};
        xyP(:, 1) = (Tgt(ii,3):Tgt(ii,4))';

        % Find and split in mz dim
        is = find(xyS(2:end-1, 2) == 0, 1, 'first');
        if ~isempty(is)

            newLim(end+1, 1:5) = [Tgt(ii,1) xyS(is+1, 1) Tgt(ii,3:4 ) 0];
            newLim(end+1, 1:5) = [xyS(is+1, 1) Tgt(ii,2:4 ) 0];
            TargetScans = unique([TargetScans; (Tgt(ii,3):Tgt(ii,4))']);

            splitedMe = true;
            continue
        end

        % Find and split in time dim
        is = find(xyP(2:end-1, 2) == 0, 1, 'first');
        if ~isempty(is)
            newLim(end+1, 1:5) = [Tgt(ii, 1:2) xyP(is+1, 1) Tgt(ii, 4) 0];
            newLim(end+1, 1:5) = [Tgt(ii, 1:3) xyP(is+1, 1) 0];
            TargetScans = unique([TargetScans; (Tgt(ii,3):Tgt(ii,4))']);

            splitedMe = true;
            continue
        end

        error("WHAT THE HECK")
    end

    if isempty(newLim), break; end
    Tgt = newLim; clear newLim XScans XProfiles;
    TargetScans(TargetScans==0) = [];

    XScans{size(Tgt, 1)} = [];
    XProfiles{size(Tgt, 1)} = [];

    % Make Extraxted profiles XProfiles and Extracted scans Xscans
    for II = 1:numel(TargetScans)
        ii = TargetScans(II);
        ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
            'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
        fidRead = fopen(ScanName, 'rb');
        myMS = fread(fidRead, ...
            [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
            "double");
        fclose(fidRead);
        fMS = Spectra.("m/z axis");

        if ~isempty(myMS)

            [~, Id1] = intersect(fMS, myMS(:,1));
            if numel(Id1) < numel(myMS(:,1)), error(""); end
            fMS(Id1, 2) = myMS(:,2);
        else
            fMS(Id1, 2) = 0;
        end

        Id2Add = find(Tgt(:,3) <= ii & Tgt(:,4) >= ii);

        for jj = 1:numel(Id2Add)
            JJ = Id2Add(jj);
            if isempty(XScans{JJ})
                XScans{JJ} = fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                XScans{JJ}(:,2) = 0 + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                XProfiles{JJ}(1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                XProfiles{JJ}(1, 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
            else
                XScans{JJ}(:,1) = XScans{JJ}(:,1) + fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                XScans{JJ}(:,2) = XScans{JJ}(:,2) + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                XProfiles{JJ}(end+1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                XProfiles{JJ}(end  , 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));

            end
        end
    end
end

%% CHECK FOR SIZE
IdT = 1;
DarkNoise = min(infoDataset.NoiseVector(:,end));
while 1
    SM = myList.Var1;
    SM(:, end+1) = (SM(:,2) - SM(:, 1)).*(SM(:, 4)-SM(:, 3));
    IdX = find(SM(:, 6) > options.method.maxSize(1));

    if isempty(IdX), break; end
    IdT = IdT+1; if IdT > numel(options.method.myThreshs), break, end
    ItThres = DarkNoise*options.method.myThreshs(IdT);

    Tgt = SM(IdX, 1:5);
    XScans{size(Tgt, 1)} = [];
    XProfiles{size(Tgt, 1)} = [];
    myList(IdX, :) = [];

    newLim = [];
    TargetScans = [];

    for ii = 1:numel(XScans)
        TargetScans = unique([TargetScans; (Tgt(ii,3):Tgt(ii,4))']);
    end

    % Make Extracted profiles XProfiles and Extracted scans Xscans
    for II = 1:numel(TargetScans)
        ii = TargetScans(II);
        ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
            'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
        fidRead = fopen(ScanName, 'rb');
        myMS = fread(fidRead, ...
            [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
            "double");
        fclose(fidRead);
        fMS = Spectra.("m/z axis");
        myMS(myMS(:,2) <= ItThres, 2) = 0;

        if ~isempty(myMS)

            [~, Id1] = intersect(fMS, myMS(:,1));
            if numel(Id1) < numel(myMS(:,1)), error(""); end
            fMS(Id1, 2) = myMS(:,2);
        else
            fMS(Id1, 2) = 0;
        end

        Id2Add = find(Tgt(:,3) <= ii & Tgt(:,4) >= ii);

        for jj = 1:numel(Id2Add)
            JJ = Id2Add(jj);
            if isempty(XScans{JJ})
                XScans{JJ} = fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                XScans{JJ}(:,2) = 0 + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                XProfiles{JJ}(1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                XProfiles{JJ}(1, 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));

            else
                XScans{JJ}(:,1) = XScans{JJ}(:,1) + fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                XScans{JJ}(:,2) = XScans{JJ}(:,2) + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                XProfiles{JJ}(end+1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                XProfiles{JJ}(end  , 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));

            end
        end
    end

    newList = table();
    while 1
        if isempty(Tgt), break; end
        Lim2Del = false(numel(XScans), 1);
        for ii = 1:numel(XScans)

            % 1. remove leading and trailing zeros
            % 1.1 Spectra (m/z profile)
            xyS = [];
            xyS(:, 2:3) =  XScans{ii};

            if numel(XScans{ii}) < 2

                Lim2Del(ii) = true;
                continue
            end

            xyS(:, 1) = (Tgt(ii,1):Tgt(ii,2))';

            if xyS(1, 2) == 0
                is = find(xyS(:,2) > 0, 1, 'first') - 1;
            else
                is = 1;
            end

            if xyS(end, 2) == 0
                ie = find(xyS(:,2) > 0, 1, 'last') + 1;
            else
                ie = length(xyS(:,2));
            end

            xyS = xyS(is:ie, :);

            if numel(xyS) < 9

                Lim2Del(ii) = true;
                continue
            end

            try
            XScans{ii} =  xyS(:, 2:3);
            Tgt(ii,1:2) = [xyS(1, 1) xyS(end, 1)];
            catch
                disp("pp")
            end

            % 1.2 Profile
            xyP = [];
            
            xyP(:, 2:3) = XProfiles{ii};

            if numel(XProfiles{ii}) < 2

                Lim2Del(ii) = true;
                continue
            end

            cAxis = (Tgt(ii,3):Tgt(ii,4))';

            cAxis(cAxis == 0) = [];
            xyP(:, 1) = cAxis';

            if xyP(1, 2) == 0
                is = find(xyP(:,2) > 0, 1, 'first') - 1;
            else
                is = 1;
            end

            if xyP(end, 2) == 0
                ie = find(xyP(:,2) > 0, 1, 'last') + 1;
            else
                ie = length(xyP(:,2));
            end

            xyP = xyP(is:ie, :);

            if numel(xyP) < 9

                Lim2Del(ii) = true;
                continue
            end

            XProfiles{ii} =  xyP(:, 2:3);



            Tgt(ii, 3:4) = [xyP(1, 1) xyP(end, 1)];


            % 2. Check if valid via size
            if ~any(XScans{ii}(:, 2) >= options.method.min_Chrom_Len) ...
                    | ~any(XProfiles{ii}(:, 2) >= options.method.min_MZ_Len)

                Lim2Del(ii) = true;
                continue
            end

            % 3. Check if scan is splitted by zeros

            if any(xyS(2:end-1, 2) == 0)
                Lim2Del(ii) = false;
                continue
            end

            % 4.2. FFind if zeros intensities split
            if any(xyP(2:end-1, 2) == 0)
                Lim2Del(ii) = false;
                continue
            end

            % If within this increment we arrive here, profile and scans are
            % correct and scan and profile are not split within our assumptions.
            % Data are recorded
            newList = [newList; {Tgt(ii, :), {XScans{ii}}, {XProfiles{ii}}}];
            Lim2Del(ii) = true;
        end

        Tgt(Lim2Del, :) = [];
        XScans(Lim2Del) = [];
        XProfiles(Lim2Del) = [];

        newLim = [];
        TargetScans = [];

        for ii = 1:numel(XScans)

            splitedMe = false;
            % 1. Make profile & Spectra
            xyS = [];
            xyS(:, 2:3) = XScans{ii};
            xyS(:, 1) = (Tgt(ii,1):Tgt(ii,2))';

            xyP = [];
            xyP(:, 2:3) = XProfiles{ii};
            xyP(:, 1) = (Tgt(ii,3):Tgt(ii,4))';

            % Find and split in mz dim
            is = find(xyS(2:end-1, 2) == 0, 1, 'first');
            if ~isempty(is)

                newLim(end+1, 1:5) = [Tgt(ii,1) xyS(is+1, 1) Tgt(ii,3:4 ) 0];
                newLim(end+1, 1:5) = [xyS(is+1, 1) Tgt(ii,2:4 ) 0];
                TargetScans = unique([TargetScans; (Tgt(ii,3):Tgt(ii,4))']);

                splitedMe = true;
                continue
            end

            % Find and split in time dim
            is = find(xyP(2:end-1, 2) == 0, 1, 'first');
            if ~isempty(is)
                newLim(end+1, 1:5) = [Tgt(ii, 1:2) xyP(is+1, 1) Tgt(ii, 4) 0];
                newLim(end+1, 1:5) = [Tgt(ii, 1:3) xyP(is+1, 1) 0];
                TargetScans = unique([TargetScans; (Tgt(ii,3):Tgt(ii,4))']);

                splitedMe = true;
                continue
            end

            error("WHAT THE HECK")
        end

        if isempty(newLim), break; end
        Tgt = newLim; clear newLim XScans XProfiles;
        TargetScans(TargetScans == 0) = [];

        XScans{size(Tgt, 1)} = [];
        XProfiles{size(Tgt, 1)} = [];

        % Make Extracted profiles XProfiles and Extracted scans Xscans
        for II = 1:numel(TargetScans)
            ii = TargetScans(II);
            ScanName = fullfile(obj.Path2Fin, obj.Datasets.Name{Id2Dataset},...
                'Scans', ['Scan#', num2str(Profiles.Scan_Index(ii)), '.dat']);
            fidRead = fopen(ScanName, 'rb');
            myMS = fread(fidRead, ...
                [Profiles.Length_MS_Scan_Profile(ii) Profiles.Column_MS_Scan_Profile(ii)],...
                "double");
            fclose(fidRead);
            fMS = Spectra.("m/z axis");
            myMS(myMS(:,2) <= ItThres, 2) = 0;


            if ~isempty(myMS)

                [~, Id1] = intersect(fMS, myMS(:,1));
                if numel(Id1) < numel(myMS(:,1)), error(""); end
                fMS(Id1, 2) = myMS(:,2);
            else
                fMS(Id1, 2) = 0;
            end

            Id2Add = find(Tgt(:,3) <= ii & Tgt(:,4) >= ii);

            for jj = 1:numel(Id2Add)
                JJ = Id2Add(jj);
                if isempty(XScans{JJ})
                    XScans{JJ} = fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                    XScans{JJ}(:,2) = 0 + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                    XProfiles{JJ}(1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                    XProfiles{JJ}(1, 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                else
                    XScans{JJ}(:,1) = XScans{JJ}(:,1) + fMS(Tgt(JJ,1):Tgt(JJ,2), 2);
                    XScans{JJ}(:,2) = XScans{JJ}(:,2) + double(fMS(Tgt(JJ,1):Tgt(JJ,2), 2) ~= 0);
                    XProfiles{JJ}(end+1, 1) = sum(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));
                    XProfiles{JJ}(end  , 2) = nnz(fMS(Tgt(JJ,1):Tgt(JJ,2), 2));

                end
            end
        end
    end

    myList = [myList; newList];
end

SM = myList.Var1;
SM(:, end+1) = (SM(:,2) - SM(:, 1)).*(SM(:, 4)-SM(:, 3));
IdX = SM(:, 6) > options.method.maxSize(2);
myList(IdX, :) = [];

%% MAKE AND SAVE THE ROIS STRUCTURES
variable_names_types = [["ID", "double"]; ...
    ["TargetTime", "double"]; ...
    ["TimeMin", "double"]; ...
    ["TimeMax", "double"]; ...
    ["Targetmz", "double"]; ...
    ["mzMin", "double"]; ...
    ["mzMax", "double"]];

% Make Target table with the limits of each ROIs
Target = table('Size', [0, size(variable_names_types,1)],...
    'VariableNames', variable_names_types(:,1),...
    'VariableTypes', variable_names_types(:,2));

for ii = 1:size(myList, 1)

    mSpectra = mzAxis(myList.Var1(ii, 1):myList.Var1(ii, 2));
    mSpectra(:, 2) = myList.Var2{ii}(:,1);

    if  myList.Var1(ii, 3) == 0
        Profile = Profiles.Scan_Time(1:myList.Var1(ii, 4));
    else
        Profile = Profiles.Scan_Time(myList.Var1(ii, 3):myList.Var1(ii, 4));
    end
    Profile(:, 2) = myList.Var3{ii}(:,1);

    %TODO: normalisation, modify to calculate area, not sum everywehere
    mSpectra(:, 2) = mSpectra(:, 2)*mean(diff(Profile(:,1)));
    mSpectra(:, 3) = 1:size(mSpectra, 1);
    Profile(:, 2) = Profile(:, 2)*mean(diff(mSpectra(:,1)));
    Profile(:, 3) = 1:size(Profile, 1);

    Target = [Target; {nan, ...
        nan, Profile(1, 1), Profile(end, 1),...
        nan, mSpectra(1, 1), mSpectra(end, 1)}];
end


Target = sortrows(Target,"TimeMin","ascend");
Target.ID = (1:size(Target, 1))';
myROIs = obj.mkMnROI(find(Id2Dataset), Target, ...
    [options.method.min_MZ_Len options.method.min_Chrom_Len IntThres]);

if options.filterROIs.do
    DesVar = myROIs.DescriptiveVars;
    Id2ROI = DesVar.IdROI;
    X = [DesVar.nnz, DesVar.CV, DesVar.Kurtosis, ...
        DesVar.ROIQual1, DesVar.ROIQual2, DesVar.LocalMaxDens, DesVar.LocalMax];

    % TAG THE ROIs
    if isempty(options.filterROIs.model)
        Y = nan(size(Id2ROI));
        Path2ROIs = myROIs.Path2ROIs;

        while 1

            rdi = randi(size(DesVar, 1));
            while ~isnan(Y(rdi))

                rdi = randi(size(DesVar, 1));
            end

            Id2X = Id2ROI(rdi);
            Tgt = myROIs.Targets(myROIs.Targets.ID == Id2X, :);

            fileName = fullfile(myROIs.Path2ROIs, ['ROI#', ...
                num2str(Id2X), '.dat']);
            [fidReadDat, errmsg]  = fopen(fileName, 'rb');
            Data = reshape(fread(fidReadDat, "double"), Tgt.sizeROI);
            fclose(fidReadDat);
            fig = figure;
            ax = uiaxes(fig);
            mesh(ax, Data);
            uiwait(fig)

            answer = questdlg('Would you like to?', ...
            	'answer  Menu', ...
            	'keep', 'discard', 'keep');

            switch answer
                case 'keep'
                    Y(rdi) = 1;

                case 'discard'
                    Y(rdi) = 0;
            end

            msgfig = msgbox (sprintf('Total number of ROis: %.0f, ROI label as good: %.0f/%.0f, ROI label as bad: %.0f/%.0f',...
                numel(Y), sum(Y == 1, 'omitnan'), options.filterROIs.minPerLabel, sum(Y == 0, 'omitnan'), options.filterROIs.minPerLabel));
            uiwait(msgfig, 1)
            try close(msgfig); catch, end

            if sum(Y == 0, 'omitnan') > options.filterROIs.minPerLabel ...
                    && sum(Y == 1, 'omitnan') > options.filterROIs.minPerLabel
                break

            end
        end

        Id2Rem = isnan(Y);
        Y_model = Y(~Id2Rem); X_model = X(~Id2Rem, :);
        options.filterROIs.LabelData = {X_model; Y_model};
        c = cvpartition(size(Y_model, 1),'Holdout', 0.33);
        B = TreeBagger(25, X_model(c.training, :), Y_model(c.training));
        predChar1 = B.predict(X_model(c.test, :));
        pd = str2double(predChar1);
        cl2 = Y_model(c.test);
        hpc = plotconfusion(cl2', pd');
        options.filterROIs.model = B;
    end

    B = options.filterROIs.model;
    options.filterROIs.FilteredData = [Id2ROI, str2double(B.predict(X))];
    myROIs.Filter = options.filterROIs;
    save(fullfile(myROIs.myPath, 'myROIs.mat'), 'myROIs')
end

AdditionalInformation.nonEndingErrors{1} = '';
AdditionalInformation.ComputingTime = toc(tStart);
AdditionalInformation.CreatedObject = myROIs;
newLineTable.TypeOfAction = 'getROIs';
newLineTable.Options4Creations = {options};
newLineTable.DateOfCreation = datetime;
newLineTable.Output = fullfile(myROIs.myPath, 'myROIs.mat');
newLineTable.AdditionalInformation = AdditionalInformation;
newLineTable = struct2table(newLineTable, 'AsArray', true);

if isempty(obj.Datasets.SecondaryActions{Id2Dataset})
    obj.Datasets.SecondaryActions{Id2Dataset} = newLineTable;

else
    obj.Datasets.SecondaryActions{Id2Dataset} = [obj.Datasets.SecondaryActions{Id2Dataset} ...
        ; newLineTable];

end
myFinnee = obj; %#ok<*NASGU>
save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')

end
