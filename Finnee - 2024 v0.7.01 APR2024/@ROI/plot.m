%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot(obj)
% Function dealing in the plotting of the profile

% Definition
foX = obj.AxisX.fo;
foY = obj.AxisY.fo;
foZ = obj.AxisZ.fo;

fig = figure('Name', obj.Title);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@myupdatefcn);
c = uicontextmenu;
uimenu(c, 'Label','ExportObj','Callback',@doExport);
fig.UIContextMenu = c;

% Main plot
infoX = obj.AxisX;
infoY = obj.AxisY;
infoZ = obj.AxisZ;

[X, Y] = meshgrid(obj.AxisX.Data(:,1), obj.AxisY.Data(:,1));
surf(Y, X,  obj.Data)
title(obj.Title, 'Interpreter', 'none');
ylabel([infoX.Label, ' / ', infoX.Units]);
xlabel([infoY.Label, ' / ', infoY.Units]);
zlabel([infoZ.Label, ' / ', infoZ.Units]);

    function doExport(~, ~)

        assignin('base', 'currentTrace', obj)
    end

    function txt = myupdatefcn(~,event_obj)
        % New data cursor update function

        pos = get(event_obj,'Position');
        xString = [infoX.Label, ' = ', num2str(pos(2), foX), ' ', infoX.Units];
        yString = [infoY.Label, ' = ', num2str(pos(1), foY), ' ', infoY.Units];
        zString = [infoZ.Label, ' = ', num2str(pos(3), foZ), ' ', infoZ.Units];
        txt = {xString, yString, zString };
    end
end