%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analyzeMe(obj, ScanRate)
% Function dealing in the plotting of the profile

option.baseline.type = "PF_HRMS";
option.baseline.pm1  = 6;
option.baseline.pm2  = 2;


switch ScanRate
    case 'low'
        wdw_lm     = 2;
        wdw_SG     = 3;
        wdw_Noise  = 5;
        min_len    = 7;



    case 'med'
        wdw_lm     = 3;
        wdw_SG     = 7;
        wdw_Noise  = 7;
        min_len    = 10;

    case 'high'
        wdw_lm     = 5;
        wdw_SG     = 15;
        wdw_Noise  = 9;
        min_len    = 15;

end

SpikeNoise = [];
MySplits = {};

% SPLIT DATA;

XY = obj.Data;
VectorIsNul = find(XY(:, 2) == 0);

%%% Find chrom portion
IdX = find(diff(VectorIsNul) > min_len+1);
for ii = 1:numel(IdX)
    MySplits{end+1} =  XY(VectorIsNul(IdX(ii)):VectorIsNul(IdX(ii)+1), :);

end

%%% Find spikes
IdX = find(diff(VectorIsNul) > 1 & diff(VectorIsNul) <= 3);
for ii = 1:numel(IdX)
    SpikeNoise(end+1) =  max(XY(VectorIsNul(IdX(ii)):VectorIsNul(IdX(ii)+1), 2));

end

SzDN = numel(SpikeNoise);
while 1
    SpikeNoise(isoutlier(SpikeNoise)) = [];
    if SzDN == numel(SpikeNoise)
        break;

    else
        SzDN = numel(SpikeNoise);

    end

end 
SpikeNoise = mean(SpikeNoise) + 2*std(SpikeNoise);

% CHECK IF BASELINE CORRECTION IS NEEDED
if numel(MySplits) == 1 & height(MySplits{1}) == height(XY)
    doBaseline = true;

else
    doBaseline = false;

end

% DO FOR EACH SPLIT SPECTRA, SMOOTHING AND NOISE
for ii = 1:numel(MySplits)
    cXY = MySplits{ii};
    disp("pp")

    FitMe = [];
    %%% Measure noise and filter
    kk = 1;
    
    for jj = 1:height(cXY)-wdw_Noise+1
        [p, S, mu] = polyfit(cXY(jj:jj+wdw_Noise-1, 1), cXY(jj:jj+wdw_Noise-1, 2), 2);
        y = polyval(p, cXY(jj:jj+wdw_Noise-1, 1), [], mu);
        FitMe(kk, jj:jj+wdw_Noise-1) = y;
        if kk == wdw_Noise
            kk = 1;

        else
            kk = kk +1;

        end
    end
    FitMe(FitMe == 0) = NaN;
    cXY(:, end+1) =  mean(FitMe, 1, 'omitnan');
    MySplits{ii} = cXY;

    NoiseVector = cXY(:, 2)- cXY(:, 5);
    SzDN = numel(NoiseVector);

    while 1
        NoiseVector(isoutlier(NoiseVector)) = [];
        if SzDN == numel(NoiseVector)
            break;

        else
            SzDN = numel(NoiseVector);

        end
    end
    Noise(ii) = mean(NoiseVector) + 2*std(NoiseVector);
end

[~,g] = sgolay(2, wdw_SG);
XY_filt = filter2(g(:, 1), XY(:, 2), 'same');
NoiseVector = (XY(:, 2)-XY_filt).^2;
[z, w] = doPF_HRMS(XY, option.baseline.pm1, 3);

% Definition
XY = obj.AxisX.XY;
% foY = obj.AxisY.fo;
% 
% fig = figure('Name', obj.FigureTitle);
% dcm_obj = datacursormode(fig);
% set(dcm_obj,'UpdateFcn',@myupdatefcn);
% c = uicontextmenu; 
% uimenu(c, 'Label','ExportObj','Callback',@doExport);
% fig.UIContextMenu = c;
% 
% % Main plot
% infoX = obj.AxisX;
% infoY = obj.AxisY;
% 
% switch upper(obj.TraceType)
%     case {'EMP', '', 'OTR'} 
%         plot(obj.Data(:,1), obj.Data(:,2));
%         title(obj.Title, 'Interpreter', 'none');
%         xlabel([infoX.Label, ' / ', infoX.Unit]);
%         ylabel([infoY.Label, ' / ', infoY.Unit]);
% 
% end
% 
%     function doExport(~, ~)
% 
%         assignin('base', 'currentTrace', obj)
%     end
% 
%     function txt = myupdatefcn(~,event_obj)
%         % New data cursor update function
% 
%         pos = get(event_obj,'Position');
%         xString = [infoX.Label, ' = ', num2str(pos(1), foX), ' ', infoX.Unit];
%         yString = [infoY.Label, ' = ', num2str(pos(2), foY), ' ', infoY.Unit];
%         txt = {xString, yString};
%     end
end