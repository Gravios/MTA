classdef MTAAPfs < MTAAnalysis

    properties (Abstract)
        path
        ext
        parameters
        mdata
        data
    end

    methods

        function Anal = MTAAPfs(Session_argin, parameters, funcHandle, varargin)
            [SessionName,MazeName,TrialName] = DefaultArgs(Session_args,{'jg05-20120317','cof','all'});
            Session = MTASession(SessionName,[],MazeName)
            Trial = MTATrial(SessionName,{{'Spk',S} 
            
            [units,states,overwrite,name,trialName,type,spk_shuffle,pos_shuffle,numBSiterations,smooth,nbins,mazeName,constrained_to_maze]=DefaultArgs(varargin,{[],'walk',0,[],'all','xy','n',0,1,0.03,50,'cof',1});
            % Use parameters to populate database            
            Anal.genPath;
            Anal.ext = 'pfs';
            
            % Pre-process data
            
            
        end
    end
    
end