classdef MTAApfs %< MTAAnalysis

    properties 
        path
        session
        tag
        ext
        parameters
        mdata
        adata
        data
    end

    methods

        function Pfs = MTAApfs(Session_args, parameters, varargin)     
            %Pfs = MTAAnalysis();
            prop = properties('MTAAPfs');
            for i = 1:length(prop)
                Pfs.(prop{i})=[];
            end
            
            %if isempty(Session_args),return,end
            
            %if isa(Session_args,'MTATrial'),
                Session = Session_args;
                SessionName = Session.name;
                MazeName    = Session.maze.name;
                TrialName   = Session.trialName;
%             elseif iscell(Session_args)
%                 [SessionName,MazeName,TrialName] = DefaultArgs(Session_args,{'jg05-20120317','cof','all'});
%                 Session = MTASession(SessionName,[],MazeName);
%                 Session = MTATrial(Session,{{'Spk',S}},TrialName);
%             end
            
            Pfs.ext = 'pfs';
            
            [units,states,overwrite,tag,type,spkShuffle,posShuffle,numIter,smooth,binDims]=DefaultArgs(varargin,{[],'walk',0,[],'xy','n',0,1,0.03,[20,20,40]});
            % Use parameters to populate database            
            
            Pfs.path = Session.spath;
            Pfs.tag  = tag;
            
            Pfs.session.name = SessionName;
            Pfs.session.trialName   = TrialName;
            Pfs.session.mazeName    = MazeName;

            pfsState = Session.stc{states,Session.xyz.sampleRate};

            Pfs.parameters.states = pfsState.label;
            Pfs.parameters.type   = type;
            Pfs.parameters.spkShuffle = spkShuffle;
            Pfs.parameters.posShuffle = posShuffle;
            Pfs.parameters.numIter  = numIter;
            Pfs.parameters.smooth   = smooth;
            Pfs.parameters.binDims = binDims;
            Pfs.parameters.constrained_to_maze = constrained_to_maze;

            Pfs.adata.trackingMarker = Session.trackingMarker;
            Pfs.adata.bins = [];
            Pfs.adata.binSizes = [];
                        
            Pfs.updateFilename(Session,tag);
            
            numUnits = numel(units);
            
            %% load existing data
            epftmp = exist(pf_tmpfile,'file');
            if epftmp&&overwrite~=1,
                load(pf_tmpfile);
                % Load specific units
                Pfs.data = Pfs.data(ismember(Pfs.data.clu,units));
                return
            elseif ~epftmp
                %% Instantiate Pfs Data Variables if DeNovoCalc
                dind = 1:numUnits;
                Pfs.data =struct('maxRate',    zeros(5,numUnits),...
                                 'maxRateInd', zeros(5,numUnits,numel(binDims)),...
                                 'maxRatePos', zeros(5,numUnits,numel(binDims)),...
                                 'rateMap',    zeros(prod(binSizes),numUnits,numIter),...
                                 'stdMap',     zeros(prod(binSizes),numUnits));
            elseif epftmp&&overwrite==1,
                load(pf_tmpfile);
                
                oldUnitInds = find(ismember(Pfs.data.clu,units));
                numOldUnits = numel(oldUnitInds);
                numNewUnits = numUnits - numOldUnits;
                tnumUnits = numNewUnits + numel(Pfs.data.clu);
                if numNewUnits>0,
                    newdata =struct( 'clu',        zeros(1,numNewUnits),...
                                     'elClu',      zeros(1,numNewUnits),...
                                     'el',         zeros(1,numNewUnits),...
                                     'maxRate',    zeros(5,numNewUnits),...
                                     'maxRateInd', zeros(5,numNewUnits,numel(binDims)),...
                                     'maxRatePos', zeros(5,numNewUnits,numel(binDims)),...
                                     'rateMap',    zeros(prod(binSizes),numUnits,numIter),...
                                     'stdMap',     zeros(prod(binSizes),numUnits));
                end
                field = fieldnames(newdata);
                for f = 1:numel(field);
                    Pfs.data = cat(2,Pfs.data.(field(f)),newdata.(field(f)));
                end
                
                dind = [oldUnitInds(:);[tnumUnits-numNewUnits:tnumUnits]'];
                
            end
            
            %% load Units into spk object;
            Session.spk.create(Session,Session.xyz.sampleRate,pfsState.label,units);
            

        end
        

        function Pfs = updateFilename(Pfs,Session,varargin)
            Pfs.tag= varargin{1};

            pfsState = Session.stc{Pfs.parameters.states,Session.xyz.sampleRate};
            if isempty(Pfs.tag)
                binDimTag = num2str(Pfs.parameters.binDims);
                binDimTag(isspace(binDimTag)&isspace(circshift(binDimTag',1)'))=[];
                binDimTag(isspace(binDimTag)) = '_';
                pf_tmpfile = fullfile(Pfs.path, [Session.filebase ...
                    '.pfs.' Pfs.parameters.type '.' Session.trackingMarker '.' pfsState.label '.' ...
                    Pfs.parameters.spkShuffle 'ps' num2str(Pfs.parameters.posShuffle) 'bs' num2str(Pfs.parameters.numIter) 'sm' num2str(Pfs.parameters.smooth*100) 'bd' binDimTag '.mat']);
            else
                pf_tmpfile = fullfile(Session.spath.analysis, [Session.filebase '.Pfs.' Pfs.tag '.mat']);
            end
        end

    end
    
end