%% DESCRIPTION
% PROFILE is the class that deals with two-dimensional representations.
% Those can either be electropherograms, chromatograms, MS spectra or
% others.
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef ROI

    properties
        Title       % Title of the Trace
        FigureTitle % Additional info (normally where the data come from'.
        % Will be displayed in the figure Title
        AxisX       % Information about AxisX (label, unit, precision)
        AxisY       % Information about AxisY (label, unit, precision)
        AxisZ       % Information about AxisY (label, unit, precision)
        Noise       % Noise vector in the m/z dimension
        BlankNoise  % minimum noise
        Data        % The Data
        AditionalInformation % Structure
    end

    properties (Dependent)
        InfoROI      % Get back the data of the Axis

    end

    methods

        function obj = ROI(container, index, target)
            % Creator method.

            if isa(container, 'struct')
                if ~isfield(container, 'Path2ROIs')
                    error("The peak list table should contained a linked to the region of interests")
                end

                obj.AditionalInformation.fullPath =  container.Path2ROIs;

                if nargin  == 3

                    obj.AditionalInformation.target = target;
                end

                iTgt = find(container.Targets.ID == index);

                if isempty(iTgt)
                    return

                end
                % load axis
                myFinnee = container.FinneeFile;
                dts = container.Dataset;
                Tgt = container.Targets(iTgt, :);

                obj.AxisX = myFinnee.Datasets.Labels{dts}.AxisX;
                obj.AxisX.Data = container.AxisX(container.AxisX >= Tgt.TimeMin & container.AxisX <= Tgt.TimeMax);

                obj.AxisY = myFinnee.Datasets.Labels{dts}.AxisY;
                obj.AxisY.Data = container.AxisY(container.AxisY >= Tgt.mzMin & container.AxisY <= Tgt.mzMax);

                obj.AxisZ = myFinnee.Datasets.Labels{dts}.AxisZ;

                minNoise = container.minNoise;
                NoiseVector = container.NoiseVector;
                obj.Noise = ones(size(obj.AxisY.Data)) * minNoise;
                for ii = 1:numel(obj.AxisX.Data)
                    cNoise = ones(size(obj.AxisY.Data)) * minNoise;

                    try
                    Id2Dn = find(NoiseVector(:, 2) <= obj.AxisX.Data(ii) & NoiseVector(:, 3) >= obj.AxisX.Data(ii));

                    catch
                        Id2Dn = [];
                    end
                    
                    if ~isempty(Id2Dn)
                        [C, Id2, Id1] = intersect(obj.AxisY.Data (:,1), NoiseVector(Id2Dn, 1));

                        if ~isempty(Id2)
                            cNoise(Id2) = NoiseVector(Id2Dn(Id1), 4);
                            obj.Noise = max(obj.Noise, cNoise);
                        end
                    end
                end

                obj.BlankNoise = minNoise;


                obj.Data = [];
                fileName = fullfile(container.Path2ROIs, ['ROI#', ...
                    num2str(index), '.dat']);
                [fidReadDat, errmsg]  = fopen(fileName, 'rb');

                DT = fread(fidReadDat, "double");

                if ~isempty(DT)
                    obj.Data = reshape(DT, Tgt.sizeROI);

                end
                
                fclose(fidReadDat);


                obj.Title = [myFinnee.FileID, ' | ',  myFinnee.Datasets.Name{dts}, ' | ', 'ROI #', num2str(index)];
                obj.FigureTitle{1} = [myFinnee.FileID, ' | ',  myFinnee.Datasets.Name{dts}, ' | ', 'ROI #', num2str(index)];
                

            else
                error("This container is not recognised")
            end
        end

        function InfoROI = get.InfoROI(obj)
            % Pass backward InfoAxis

            InfoROI.Title     = obj.Title;
            InfoROI.FT        = obj.FigureTitle;
            InfoROI.TT        = obj.TraceType;
            InfoROI.AxisX     = obj.AxisX;
            InfoROI.AxisY     = obj.AxisY;
            InfoROI.P2Fin     = obj.Path2Fin;
            InfoROI.Loc       = obj.DataStorage;
            InfoROI.AdiPrm    = obj.AdiParam;
            InfoROI.Precision = obj.Precision;
        end
    end

end

