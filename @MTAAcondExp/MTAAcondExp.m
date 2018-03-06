classdef MTAAcondExp < hgsetget 
    properties        
        path
        filename       
        Trial
        Data
        parameters        
        meta
        data
        hash
        tag
        ext        
    end
    
    methods
        % function plot(Pfs,unit,mask)
        % MAIN args
        % Trial
        % Data
        % Spk
        % State
        % parameters
        %   .binDims
        %   .smoothingWeights
        %   .randomizationMethod
        %
            
        % PARAMETERS struct:
        %    .state
        %    .binDims
        %    .smoothingWeights
        %    .randomizationMethod
        %    .numIterations
        %    
            
        % HASH inputs
        %    Trial.filebase
        %    Data.hash
        %    Spk.hash
        %    State.hash
        %    parameters
        %
        
        function Obj = MTAAcondExp(Trial,varargin)
            
            Obj.path = load('MTAPaths.mat');
            switch class(obj),
                case 'MTATrial',
                    Trial = Obj;
                    
                    
                    if xyzp.isempty,
                        try,
                            xyz = preproc_xyz(Trial,'trb');
                        catch err,
                            disp(err)
                            xyz = preproc_xyz(Trial);
                            trackingMarker = 'hcom';
                        end
                        xyz.data = sq(xyz(:,trackingMarker,1:numel(binDims)));
                    else
                        xyz = xyzp;
                    end
                    
                    Pfs.path = Trial.spath;
                    Pfs.tag  = tag;
                    
                    Pfs.session.sessionName = TrialName;
                    Pfs.session.trialName   = TrialName;
                    Pfs.session.mazeName    = MazeName;

                    if ischar(states),
                        Pfs.parameters.states = states;
                        pfsState = Trial.stc{states,xyz.sampleRate}.copy;
                    elseif isa(states,'MTAData'),
                        pfsState = states.copy;
                        pfsState.resample(xyz);
                        Pfs.parameters.states = pfsState.label;                       
                    end
                    
                    Pfs.parameters.type   = type;
                    Pfs.parameters.spkShuffle = spkShuffle;
                    Pfs.parameters.posShuffle = posShuffle;
                    Pfs.parameters.numIter  = numIter;
                    if isempty(SmoothingWeights)&&numel(binDims)==numel(SmoothingWeights)
                        SmoothingWeights = Nbin./30;
                    end
                    Pfs.parameters.smoothingWeights   = SmoothingWeights;
                    Pfs.parameters.binDims = binDims;
                    Pfs.parameters.bootstrap = bootstrap;
                    
                    Pfs.adata.trackingMarker = trackingMarker;
                    Pfs.adata.bins = [];
                    if isempty(boundaryLimits),boundaryLimits = Trial.maze.boundaries(1:numel(type),:);end
                    Pfs.adata.binSizes = round(abs(diff(boundaryLimits(1:numel(binDims),:),1,2))./binDims');
                    
                    Pfs.data =struct( 'clu',        [],...
                        'elClu',      [],...
                        'el',         [],...
                        'maxRate',    [],...
                        'maxRateInd', [],...
                        'maxRatePos', [],...
                        'rateMap',    [],...
                        'meanRate',   [],...
                        'si',         [],...
                        'spar',       []);
                    Pfs.update_filename(Trial,tag);
                    Pfs.ext = 'pfs';
            
        end

    end
    
end