classdef MTAAcondExp < hgsetget 
    properties        
        path
        filename       
        hash
        tag
        ext
        Trial
        Data
        spkOpts
        parameters        
        meta
        data
    end
    
    methods
        function Obj = MTAAcondExp(Trial,varargin)

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
            
        % HASH struct:
        %    State.hash
        %    Data.hash
        %    Spk.hash
        %    parameters
            
        end

    end
    
end