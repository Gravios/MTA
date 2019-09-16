classdef MTASession < hgsetget
% MTASession(name,varargin) 
% Data structure to organize the analysis of neural and spatial data.
%
%   name - string: Same name as the directory of the session
%
%   varargin:
%     [mazeName, overwrite, TTLValue, xyzSystem, ephySystem, xyzSampleRate]
%
%     mazeName:        string, 3-4 letter name of the testing arena 
%                      (e.g. 'rof' := rectangular open-field)
%
%     overwrite:       boolean, flag to overwrite saved Sessions
%
%     TTLValue:        string, used to synchronize position and electrophysiological data
%
%     dataLogger:       string, name/id of the system(s) to record subject
%
%     xyzSampleRate:   numeric, samples per second for xyz tracking
%                      Vicon(M): 119.881035 Hz
%                      Vicon(V): 199.997752 Hz or 149.9974321 Hz
%                      Optitrack: duno
%
%-------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Session,
%     Session = MTASession(name,mazeName);
%     
%     Create new session,
%     Session = MTASession(name,mazeName,overwrite,TTLValue,xyzSystem,ephySystem);
%
%-------------------------------------------------------------------------------------------------
%     examples:
%       load saved session,
%         Session = MTASession('jg05-20120309','rof');
%
%       Create New Session
%         Session = MTASession('jg05-20120309','rof',1,'0x0040',{'nlx','vicon'},119.881035);
%   
%-------------------------------------------------------------------------------------------------

    properties (SetAccess = public)

        %filebase - string: full file head in the format name.(Maze.name).trialName
        filebase = '';

        %spath - struct: same as path but with Session name appended to the end
        spath
        
        %path - struct: holds all paths of the constructed data tree created by MTAConfiguration.m
        path
        
        %name - string: name of the directory holding session information
        name = '';

        %trialName - string: designation of trial the full Session has the default name 'all'
        trialName = '';
        
        %par - struct: contains parameter information regarding the recording systems, units ect...
        parameters
        
        %sampleRate - double: Sample Rate of electrophysiological recording system
        sampleRate     
        
        %Maze - MTAMaze: Object containing all maze information
        maze

        %Model - MTAModel: Object contianing all marker information
        model          

        %sync - MTADepoch: Loading periods relative to the primary recording system
        sync

        %xyz - MTADxyz: (Time,Marker,Dimension) XYZ position of each marker 
        xyz 

        %ang - MTADang: (Time,Marker,Marker,Dimension) marker to marker angles 
        ang 

        %stc - MTAStateCollection: Stores the periods of events
        stc
        
        %spk - MTASpk: Stores neurons' action potential timing and features
        spk

        %Fet - MTAFet: Object containing behavioral features
        fet

        %lfp - MTADlfp: (Time,channel) local field potential
        lfp

        %ufr - MTADufr: (Time,ClusterId) unit firing rates with lfpSampleRate
        ufr = [];      

        %nq - struct: clustering quality and spike waveform characteristics
        nq = {};

        %fbr - MTADfbr: (Time,ClusterId) fiber photometery with lfpSampleRate
        fbr = MTADfbr([]);

    end
    
    properties (SetAccess = protected)
        % parametersHash - String: used to detect changes to parameter structure before saving to file
        parametersHash = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz';
        
    end
    

    methods

        function Session = MTASession(varargin)
        % Session Constructor
            
            % DEFARGS ----------------------------------------------------------------------------
            defargs = struct('name',                      [],                                  ...
                             'mazeName',                  'cof',                               ...
                             'overwrite',                 false,                               ...
                             'TTLValue',                  '0x0040',                            ...
                             'dataLoggers',               {{'nlx','vicon'}},                   ...
                             'xyzSampleRate',             []                                   ...
            );

            [name,mazeName,overwrite,TTLValue,dataLoggers,xyzSampleRate] = ...
                DefaultArgs(varargin,defargs,'--struct');
            %-------------------------------------------------------------------------------------
            % LOAD paths from configuration folder, setup in MTAstartup->MTAConfiguration
            Session.path = load('MTAPaths.mat');
            
            
            if isa(name,'MTASession'),
            % COPY MTASession object from 'name' to 'Session'
                prop = properties('MTASession');
                for i = 1:length(prop)
                    if ismethod(name.(prop{i}),'copy')
                        Session.(prop{i})=name.(prop{i}).copy;
                    else
                    Session.(prop{i})=name.(prop{i});
                    end
                end
            elseif isempty(name), return;                 
            else
            % CREATE new session object
                Session.name = name;
                Session.spath = fullfile(Session.path.project,Session.name);
                Session.trialName = 'all';
                Session.maze = MTAMaze(mazeName);
                Session.filebase = [Session.name '.' Session.maze.name '.' Session.trialName];
                if exist(fullfile(Session.spath,[Session.filebase, '.ses.mat']),'file')...
                        &&~overwrite
                    Session = Session.load();
                    Session = Session.updatePaths();
                    %Session.xyz.load(Session.sync);
                elseif overwrite
                    varargin(1:3) = [];
                    warning(['Overwriting Session: ' ...
                             fullfile(Session.spath, [Session.filebase, '.ses.mat'])])
                    Session.create(varargin{:});
                else
                    varargin(1:3) = [];                    
                    warning(['Session with maze, ' ...
                             Session.maze.name ', does not exist: ' ...
                             'attempting to creating session']);
                    Session.create(varargin{:});
                end  
            end
        end%function MTASession
    end%methods

    methods (Static)
        function Session = validate(Session)
            
            if isa(Session,'MTATrial'),
                Session = MTASession.validate(Session.filebase);
            elseif isa(Session,'MTASession'),
                return;
            elseif ischar(Session),
                pat =['(?<sessionName>[a-z_A-Z]+\d{2,4}[-]\d{8,8}[a-z_A-Z]{0,1})\.'...
                      '(?<mazeName>\w+)\.'...
                      '(?<trialName>\w+)'];
                tok = regexp(Session,pat,'names');
                if ~isempty(tok),
                    Session = tok;
                    Session = MTASession.validate(Session);
                else
                    Session = MTASession(Session);                                    
                end
                
            elseif isstruct(Session),
                stcMode = '';
                if isfield(Session,'stcMode'),
                    stcMode = Session.stcMode;
                end

                Session = MTASession(Session.sessionName,...
                                     Session.mazeName,...
                                     false);

                if ~isempty(stcMode),
                    try,
                        Session.load('stc',stcMode);                        
                    catch err
                        disp(err);
                    end

                end

            else
                error('MTA:validate_trial: unrecognized format');
            end

        end%function validate
    end%methods(static)
    
end%class MTASession

