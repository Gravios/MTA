classdef MTATrial < MTASession
% MTATrial(name,varargin)
% 
% A subclass of MTASession structured to organize subsets of a parent session
% to aid in the analysis of neural and spatial data.
% 
% Inputs ---------------------------------------------------------------------------------
%
%  name - string: Same name as the directory of the session
%   
%  -or- 
%
%  name - MTASession: The Session object from which data will be copied
%
%
%  varargin: [trialName,mazeName,overwrite,sync]
%
%    trialName:       string,       designation of trial, the Trial representing 
%                                   the full Session has the default name 'all'
%
%    mazeName:        string,       3-4 letter name of the testing arena 
%                                   (e.g. 'rof' := rectangular open-field)
%
%    overwrite:       boolean,      flag to overwrite saved Trials
%
%    sync:            MTADepoch,    New sync object who data field
%                                   contains new synchronization periods
%                     numericArray, New periods to select from the
%                                   MTASession object in secords
%
%-- General Use --------------------------------------------------------------------------
%
%     Load from saved Trial,
%     Trial = MTATrial(name,trialName);
%
%     Parse and load by filebase
%     Trial = MTATrial.validate(filebase);
%     
%     Create new Trial,
%     Trial = MTATrial(name,trialName,mazeName,overwrite,sync);
%
%-- Examples -----------------------------------------------------------------------------
%
%     load saved Trial,
%       Trial = MTATrial('jg05-20120309','cof','all');
%
%     Parse and load saved Trial,
%       Trial = MTATrial.validate('jg05-20120309.cof.all');
%
%     Create New Trial from a subset of the total session
%       Trial = MTATrial('jg05-20120309','rof','crt1',false,[1,1000;1100,2100])
%              
%-- Notes --------------------------------------------------------------------------------
%
%  Ext: *.trl.mat
%
%  Saved Trial object contain only basic information required to load subsets of a full 
%  session
%
%
    methods 
        function Trial = MTATrial(Session,varargin)
            [mazeName,trialName,overwrite,sync] = DefaultArgs(varargin,{'cof','all',0,[]});
            if ~isa(Session,'MTASession'),
                Session = MTASession(Session,mazeName);
            end
            if strcmp(trialName,'nil'),overwrite = 1;end
            
            Trial = Trial@MTASession(Session,mazeName);
            
            Trial.trialName = trialName;
            Trial.filebase = [Trial.name '.' Trial.maze.name '.' Trial.trialName];

            % LOAD ? Trial file            
            if ~overwrite && exist(fullfile(Trial.spath, [Trial.filebase '.trl.mat']),'file') 

                ds = load(fullfile(Trial.spath, [Trial.filebase '.trl.mat']));
                Trial.sync = ds.sync;
                try
                    % LOAD ? saved state collection by mode
                    if isfield(ds,'stcmode')
                        if ~isempty(ds.stcmode)
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
                    % LOAD ? subject
                    if isfield(ds,'subject')
                        if ~isempty(ds,'subject')
                            % TODO 
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
            

            Trial.stc.updateSync(Trial.sync.copy);
            Trial.stc.updatePath(Trial.spath);
            props = properties(Trial);
            for p = 1:numel(props),
                if strcmp(props{p},'sync'),continue,end
                prop = Trial.(props{p});
                                
                if isa(prop,'MTAData'),
                    prop.resync(Trial);
                elseif isa(prop,'MTAStateCollection')
                    if ~prop.isempty,
                        for s = 1:numel(prop.states(:)),
                            prop.states{s}.resync(Trial);
                        end
                    end
                end
                
            end
        end
        

        

    end
    
    methods (Static)
        function Trial = validate(Trial)
        %function Trial = validate(Trial)
        %
        % Validate MTATrial object. If Trial exists as a
        % non-MTATrial object attempt to load Trial from disk based
        % on the class and contents of Trial.
        %            
        % acceptable input types: char, struct, cell
        %
        %  char:
        %  
        %    Session name
        %    Trial = MTATrial.validate('jg05-20120317');
        %
        %    Trial filebase
        %    Trial = MTATrial.validate('jg05-20120317.cof.all');
        %
        %  struct:
        %    sessionList = get_session_list('jg05');
        %    Trial = MTATrial.validate(sessionList(1));
        %
        %  cell:
        %    Session name
        %    Trial = {'jg05-20120317'};
        %
        %    Trial filebase elements
        %    Trial = {'jg05-20120317','cof','all'};
        %
        %    Full Trial Constructor
        %    Trial = {'jg05-20120309','crt1','rof',false,[1,1000;1100,2100]};
        %
        %
        %    Trial = MTATrial.validate(Trial);
        %
            
            if isa(Trial,'MTATrial'),
                return
            elseif isa(Trial,'MTASession'),
                Trial = MTATrial(Trial);
                return
                
            elseif ischar(Trial),
                pat =['(?<sessionName>[a-z_A-Z]+\d{2,4}[-]\d{8,8}[a-z_A-Z]{0,1})\.'...
                      '(?<mazeName>\w+)\.'...
                      '(?<trialName>\w+)'];
                tok = regexp(Trial,pat,'names');
                if ~isempty(tok),
                    Trial = MTATrial.validate(tok);
                else
                    Trial = MTATrial(Trial);                                    
                end
                return

                
            elseif iscell(Trial) && numel(Trial)==1,
                if isa(Trial{1},'MTASession')
                    Trial = Trial{1};
                elseif ischar(Trial{1})
                    pat =['(?<sessionName>[a-z_A-Z]+\d{2,8}[-]\d{8,8}[a-z_A-Z]{0,1})\.'...
                          '(?<mazeName>\w+)\.'...
                          '(?<trialName>\w+)'];
                    tok = regexp(Trial,pat,'names');
                    if ~isempty(tok),
                        Trial = MTATrial.validate(tok);
                    else
                        Trial = MTATrial(Trial);                                    
                    end
                end    
                return

                
            elseif isstruct(Trial),
                stcMode = '';
                if isfield(Trial,'stcMode'),
                    stcMode = Trial.stcMode;
                end
                meta = [];
                if isfield(Trial,'subject')
                    meta = Trial.subject;
                end
                

                try
                    Trial = MTATrial(Trial.sessionName,...
                                     Trial.mazeName,...
                                     Trial.trialName);
                catch err
                    disp(err)
                    disp('MTATrial:validate - Attempting to construct trial.')
                    QuickTrialSetup(Trial);
                    try
                        Trial = MTATrial(Trial.sessionName,...
                                         Trial.mazeName,...
                                         Trial.trialName);
                    catch err
                        disp(err)
                        disp('MTATrial:validate - Failed to construct trial.')
                    end
                    

                end
                
                if ~isempty(stcMode),
                    try,
                        Trial.load('stc',stcMode);
                    catch err
                        disp(err);
                    end
                end
                
                Trial.meta = meta;

            else
                error('MTA:validate_trial: unrecognized format');
            end

        end
    end

end
