classdef MTAApfs < MTAAnalysis

    properties (Abstract)
        path
        ext
        parameters
        mdata
        data
    end

    methods

        function Pfs = MTAApfs(Session_args, parameters, funcHandle, varargin)
            prop = properties('MTAAPfs');
            for i = 1:length(prop)
                PlaceField.(prop{i})=[];
            end
            
            if isempty(Session_args),return,end
            
            if isa(Session_args,'MTATrial'),
                Session = Session_args;
                SessionName = Session.name;
                MazeName    = Session.Maze.name;
                TrialName   = Session.trialName;
            elseif iscell(Session_args)
                [SessionName,MazeName,TrialName] = DefaultArgs(Session_args,{'jg05-20120317','cof','all'});
                Session = MTASession(SessionName,[],MazeName);
                Session = MTATrial(Session,{{'Spk',S}},TrialName);
            end
            
            Pfs.ext = 'pfs';
            
            [units,states,overwrite,name,type,spkShuffle,posShuffle,numIter,smooth,binDims,constrained_to_maze]=DefaultArgs(varargin,{[],'walk',0,'all','xy','n',0,1,0.03,[20,20,40],1});
            % Use parameters to populate database            
            Pfs.genPath;

            Pfs.mdata.sessionName = SessionName;
            Pfs.mdata.trialName   = TrialName;
            Pfs.mdata.mazeName    = MazeName;
            Pfs.mdata.name        = name;
            

            pfsState = Session.stc{states,s.xyz.sampleRate};
            
            
            Pfs.parameters.states = pfsState.label;
            Pfs.parameters.type   = type;
            Pfs.parameters.spkShuffle = spkShuffle;
            Pfs.parameters.posShuffle = posShuffle;
            Pfs.parameters.numIter  = numIter;
            Pfs.parameters.smooth   = smooth;
            Pfs.parameters.binDims = binDims;
            Pfs.parameters.constrained_to_maze = constrained_to_maze;


            
            %% Load State specific periods
            pfsState = Session.stc{states,s.xyz.sampleRate};
            if isempty(name)
                pf_tmpfile = [Session.spath.analysis Session.filebase ...
                    '.pfs.' type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    spkShuffle num2str(posShuffle) 'bs' num2str(numIter) 'sm' num2str(smooth*100) 'bn' num2str(nbins) '.mat'];
            else
                pf_tmpfile = [Session.spath.analysis Session.filebase '.Pfs.' name '.mat'];
            end

            
            
            % Pre-process data

        end
        
        function Pfs = genPath(Pfs)
        end
    end
    
end