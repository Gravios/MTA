classdef MTATrial < MTASession
% MTATrial(name,varargin) - Is a subclass of MTASession structured to organize 
%   subsets of a parent session to aid in the analysis of neural and spatial data.
%
%   Ext: *.trl.mat
%   Saved Trial object contain only basic information required to load subsets of
%   a full session
% 
%   name - string: Same name as the directory of the session
%   
%   -or- 
%
%   name - MTASession: The Session object from which data will be copied
%
%
%   varargin:
%     [preLoadedFields,trialName,new_xyzPeriods,overwrite,mazeName]
%
%     trialName:       string, designation of trial, the Trial representing 
%                      the full Session has the default name 'all'
%
%     new_xyzPeriods:  numericArray, new Periods defined in the refrence 
%                      frame of the MTASession objects xyzPeriods
%
%     overwrite:       boolean, flag to overwrite saved Trials
%
%     mazeName:        string, 3-4 letter name of the testing arena 
%                      (e.g. 'rof' := rectangular open-field)
%
%
%---------------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Trial,
%     Trial = MTATrial(name,trialName);
%     
%     Create new Trial,
%     Trial = MTATrial(name,trialName,new_xyzPeriods,overwrite,mazeName);
%
%---------------------------------------------------------------------------------------------------------
%     examples:
%       load saved Trial,
%         Trial = MTATrial('jg05-20120309','all');
%
%       Create New Trial from a subset of the total session
%         Trial = MTATrial('jg05-20120309','crt1',[1,10000;11000,21000],1,'rof','normal')
%              
%---------------------------------------------------------------------------------------------------------

    methods 
        function Trial = MTATrial(Session,varargin)
            [trialName,sync,overwrite,mazeName] = DefaultArgs(varargin,{'all',[],0,'cof'});
            if ~isa(Session,'MTASession'),
                Session = MTASession(Session,mazeName);
            end
            if strcmp(trialName,'nil'),overwrite = 1;end
            
            Trial = Trial@MTASession(Session,mazeName);
            
            Trial.trialName = trialName;
            Trial.filebase = [Trial.name '.' Trial.maze.name '.' Trial.trialName];
            Trial.stc.updateFilename(Trial.filebase);
            
            if exist(fullfile(Trial.spath, [Trial.filebase '.trl.mat']),'file')&&~overwrite
                ds = load(fullfile(Trial.spath, [Trial.filebase '.trl.mat']));
                Trial.sync = ds.sync;
                if isfield(ds,'stcmode'),
                    if ~isempty(ds.stcmode)&&exist(fullfile(Trial.spath, [Trial.filebase '.stc.' ds.stcmode '.mat']),'file')
                        Trial.stc.updateFilename([Trial.filebase '.stc.' ds.stcmode '.mat']);
                        Trial.stc.load;
                    end
                end
            else
                switch class(sync)
                    case 'MTADepoch'
                        Trial.sync = sync;
                    case 'double'
                        msg.message = 'The provided synchronization periods are empty.';
                        msg.identifier = 'MTATrial:MTAtrial:EmptySync';
                        msg.stack = dbstack;
                        if isempty(sync),error(msg),end
                        
                    otherwise
                        error('sync format not recognized')
            end
            
            Trial.trackingMarker = Session.trackingMarker;
            Trial.resync(Trial.xyz);
            Trial.stc.updateSync(Trial.sync);
            props = properties(Trial);
            for p = 1:numel(props),
                prop = Trial.(props{p});
                if isa(prop,'MTAData'),
                    if prop.isempty||strcmp(props{p},'xyz'),
                        continue,
                    else
                        Trial.resync(prop);
                    end
                elseif isa(prop,'MTAStateCollection')
                    prop.sync = Trial.sync;
                    if ~prop.isempty
                        for s = numel(prop(:)),
                            Trial.sync.resync(prop{s});
                        end
                    end
                end
            end
        end
        

        function save(Trial)
            stcmode = [];
            if ~isempty(Trial.stc)
                stcmode = Trial.stc.mode;
            end
            trialName = Trial.trialName;
            sync = Trial.sync;
            save(fullfile(Trial.spath,[Trial.filebase '.trl.mat']),'trialName','sync','stcmode');            
        end

    end
end
