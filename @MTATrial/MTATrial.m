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
%     [trialName,mazeName,overwrite,sync]
%
%     trialName:       string,       designation of trial, the Trial representing 
%                                    the full Session has the default name 'all'
%
%     mazeName:        string,       3-4 letter name of the testing arena 
%                                    (e.g. 'rof' := rectangular open-field)
%
%     overwrite:       boolean,      flag to overwrite saved Trials
%
%     sync:            MTADepoch,    New sync object who data field
%                                    contains new synchronization periods
%                      numericArray, New periods to select from the
%                                    MTASession object in secords
%
%
%---------------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Trial,
%     Trial = MTATrial(name,trialName);
%     
%     Create new Trial,
%     Trial = MTATrial(name,trialName,mazeName,overwrite,sync);
%
%---------------------------------------------------------------------------------------------------------
%     examples:
%       load saved Trial,
%         Trial = MTATrial('jg05-20120309','all');
%
%       Create New Trial from a subset of the total session
%         Trial = MTATrial('jg05-20120309','crt1','rof',false,[1,1000;1100,2100])
%              
%---------------------------------------------------------------------------------------------------------

    methods 
        function Trial = MTATrial(Session,varargin)
            [trialName,mazeName,overwrite,sync] = DefaultArgs(varargin,{'all','cof',0,[]});
            if ~isa(Session,'MTASession'),
                Session = MTASession(Session,mazeName);
            end
            if strcmp(trialName,'nil'),overwrite = 1;end
            
            Trial = Trial@MTASession(Session,mazeName);
            
            Trial.trialName = trialName;
            Trial.filebase = [Trial.name '.' Trial.maze.name '.' Trial.trialName];
            
            if exist(fullfile(Trial.spath, [Trial.filebase '.trl.mat']),'file')&&~overwrite
                
                ds = load(fullfile(Trial.spath, [Trial.filebase '.trl.mat']));
                Trial.sync = ds.sync;
                try
                 if isfield(ds,'stcmode'),
                     if ~isempty(ds.stcmode),
                         Trial.stc.updateFilename([Trial.filebase,'.stc.' ds.stcmode '.mat']);
                         Trial.stc.updateMode(ds.stcmode);
                         if exist(Trial.stc.fpath,'file')
                             Trial.stc.load;
                         else
                             Trial.stc.updateMode('default');                  
                             Trial.stc.load;
                         end                    
                     end
                 end
                end             
            else
                
                msg.message = 'The provided synchronization periods are empty.';
                msg.identifier = 'MTATrial:MTAtrial:EmptySync';
                Trial.stc.updateFilename([Trial.filebase,'.stc.default.mat']);                
                switch class(sync)
                    
                  case 'MTADepoch'
                    if sync.isempty,
                        msg.stack = dbstack;
                        error(msg);
                    else
                        Trial.sync = sync.copy;
                    end
                    
                  case 'double'

                        if isempty(sync),
                            msg.stack = dbstack;
                            error(msg),
                        else
                            Trial.sync.resample(1);
                            Trial.sync = MTADepoch(Trial.spath,[Trial.filebase '.sync.mat'],sync,1,Trial.sync.sync,0);
                        end
                        
                    otherwise
                        error('sync format not recognized')
                end
                
            end
            
            Trial.trackingMarker = Session.trackingMarker;
            Trial.stc.updateSync(Trial.sync.copy);
            Trial.stc.updatePath(Trial.spath);
            props = properties(Trial);
            for p = 1:numel(props),
                if strcmp(props{p},'sync'),continue,end
                prop = Trial.(props{p});
                                
                if isa(prop,'MTAData'),
                    Trial.resync(prop);
                elseif isa(prop,'MTAStateCollection')
                    if ~prop.isempty,
                        for s = 1:numel(prop.states(:)),
                            Trial.resync(prop.states{s});
                        end
                    end
                end
                
            end
        end
        

        

    end
    
    methods (Static)
        function Trial = validate(Trial)
            
            if isa(Trial,'MTATrial'),
                return;
            elseif ischar(Trial),
                Trial = MTATrial(Trial);
            elseif iscell(Trial),
                Trial = MTATrial(Trial{:});
            elseif isstruct(Trial),
                Trial = MTATrial(Trial.sessionName,...
                                 Trial.trialName,...
                                 Trial.mazeName);
            else
                error('MTA:validate_trial: unrecognized format');
            end

        end
    end

end
