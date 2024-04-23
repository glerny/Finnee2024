%% DESCRIPTION

function LocMax = LocalMaxima3D(XY, DXY, smoothMe)

%% CORE OF THE FUNCTION
% Finding non nul local maxima

if nargin == 2
    smoothMe = 'None';
end

switch smoothMe
    case 'None'

    case 'SmoothLow'
        h = sgsdf_2d(-1:1, -1:1, 1, 1, 0);
        XY = filter2(h, XY, 'same');

    case 'SmoothMed'
        h = sgsdf_2d(-2:2, -2:2, 1, 1, 0);
        XY = filter2(h, XY, 'same');

    case 'SmoothHigh'
        h = sgsdf_2d(-5:5, -5:5, 2, 2, 0);
        XY = filter2(h, XY, 'same');

    otherwise
        error('Smooth paramter not recognised')
end

LocMax = [];
for ii = 2:size(XY, 1)-1
    for jj = 2:size(XY, 2)-1
        Int = [max(1, ii-DXY(1)) min(ii+DXY(1), size(XY, 1)) ...
            max(1, jj-DXY(2)) min(jj+DXY(2), size(XY, 2))];
        if XY(ii, jj) == max(XY(Int(1):Int(2), Int(3):Int(4)), [], 'all', 'omitnan') && ...
                XY(ii, jj) > 0
            LocMax(end+1, :) = [XY(ii, jj), ii, jj];

        end
    end
end
LocMax = sortrows(LocMax, -1);

