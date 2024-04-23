%% DESCRIPTION
% TRACE is the class that deals with two-dimensional representations. 
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


classdef Spectrum
    
    properties
        Title       % Title of the Trace
        FigureTitle % Additional info (normally where the data come from'. 
        % Will be displayed in the figure Title 
        TraceType   % information about the typr of trace: 
        AxisX       % Information about AxisX (label, unit precision)
        AxisY       % Information about AxisY (label, unit precision)
        Data        % The Data
        AditionalInformation % Structure
        Precision   % precision of
    end
    
    properties (Dependent)
        InfoTrc      % Get back the data of the Axis
    end
    
    methods
        
        function obj = Spectrum(infoTrc, data2write)
            % Creator method.
            
            if nargin == 0
                obj.Title = '';
                obj.FigureTitle = '';
                obj.TraceType = 'EMP';
                obj.AxisX = Axis;
                obj.AxisY = Axis;
                obj.Data = [];
                obj.AditionalInformation  = struct();
                
            else
                obj.Title      	= infoTrc.Title;
                obj.FigureTitle = infoTrc.FT;
                obj.TraceType   = infoTrc.TT;
                obj.AxisX       = infoTrc.AxisX;
                obj.AxisY       = infoTrc.AxisY;
                obj.AditionalInformation    = infoTrc.AdiPrm;
                
                % control of data2write
                if nargin == 1
                    obj.Data = [];
                    obj.TraceType   = 'EMP';
               
                end

                if isempty(data2write)
                    obj.Data = [];
                    obj.TraceType   = 'EMP'; 

                else
                    obj.Data = data2write; 

                end
            end
                            
        end
        
        
        function InfoTrc = get.InfoTrc(obj)
            % Pass backward InfoAxis
            
            InfoTrc.Title     = obj.Title;
            InfoTrc.FT        = obj.FigureTitle;
            InfoTrc.TT        = obj.TraceType;
            InfoTrc.AxisX     = obj.AxisX;
            InfoTrc.AxisY     = obj.AxisY;
            InfoTrc.P2Fin     = obj.Path2Fin;
            InfoTrc.Loc       = obj.DataStorage;
            InfoTrc.AdiPrm    = obj.AdiParam;
            InfoTrc.Precision = obj.Precision;
        end
    end
    
end

