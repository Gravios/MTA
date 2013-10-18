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
            
            [units,states,overwrite,tag,type,spkShuffle,posShuffle,numIter,smooth,binDims,constrained_to_maze]=DefaultArgs(varargin,{[],'walk',0,'all','xy','n',0,1,0.03,[20,20,40],1});
            % Use parameters to populate database            
            
            Pfs.path = Session.spath;
            Pfs.tag        = tag;
            
            Pfs.session.name = SessionName;
            Pfs.session.trialName   = TrialName;
            Pfs.session.mazeName    = MazeName;

            pfsState = Session.stc{states,s.xyz.sampleRate};

            Pfs.parameters.states = pfsState.label;
            Pfs.parameters.type   = type;
            Pfs.parameters.spkShuffle = spkShuffle;
            Pfs.parameters.posShuffle = posShuffle;
            Pfs.parameters.numIter  = numIter;
            Pfs.parameters.smooth   = smooth;
            Pfs.parameters.binDims = binDims;
            Pfs.parameters.constrained_to_maze = constrained_to_maze;

            Pfs.adata.bins = [];
                        
            Pfs.updateFileName(tag);
            
            % load existing data
            if exist(pf_tmpfile,'file')&&overwrite~=1,
                load(pf_tmpfile);
                return
            end
            Session.spk.create(Session.xyz.sampleRaet,pfsState.label,units);
            
        end
        

        function Pfs = updateFilename(Pfs,Session,varargin)
            Pfs.tag= varargin{1};

            pfsState = Session.stc{states,s.xyz.sampleRate};
            if isempty(Pfs.tag)
                binDimTag = num2str(binDims);
                binDimTag(isspace(binDimTag)) = '_';
                pf_tmpfile = fullfile(Pfs.path, [Session.filebase ...
                    '.pfs.' type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    spkShuffle num2str(posShuffle) 'bs' num2str(numIter) 'sm' num2str(smooth*100) 'bd' binDimTag '.mat']);
            else
                pf_tmpfile = fullfile(Session.spath.analysis, [Session.filebase '.Pfs.' Pfs.tag '.mat']);
            end

    end
    
end