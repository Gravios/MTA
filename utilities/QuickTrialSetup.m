function varargout = QuickTrialSetup(Session,varargin)
%function QuickTrialSetup(Session,varargin)
%
% Session:        MTASession,char,struct
%                 Basis for trial setup. Either the MTASession
%                 ojbect, a string containing the Session name or a
%                 list of structs containing trial constructor
%                 info. (see get_session_list)
%
% Variables:
%   trialName:      string, defarg  - 'all' 
%                   Name of the new Trial
%
%   mazelName:      string, defarg  - 'cof' 
%                   Name of the maze
%
%   offsets:        matrix, defarg  - [0,0]  
%                   number of seconds to skip/clip from begining/end of
%                   each motion tracking trial
%                   see MTA:MTADepoch for more details
%
%   dropSyncInd,    matrix, defarg  - []    
%                   Indicies of the Motion Tracking segments to ignore, in
%                   order of the Session.xyz.sync.data
%
%   includeSyncInd,    matrix, defarg  - []    
%                   Indicies of the Motion Tracking segments to include, in
%                   order of the Session.xyz.sync.data
%
%   autolabel:          Logical, defarg  - true
%                   Try to do some automatic behavior labeling
%
%   debug:          Logical, defarg  - false
%                   Display some diagnostic plots
%
%   host: ???
%
%   local: ???
%
%   overwrite:      Logical, defarg - false
%                   Overwrite current version of Trial
%


% DEFARGS ----------------------------------------------------------------------
HostConf = load('MTAConf');

defargs = struct( 'trialName',  'all',                                       ...
                   'mazeName',  'cof',                                       ...
                    'offsets',  [0,0],                                       ...
                'dropSyncInd',  [],                                          ...
             'includeSyncInd',  [],                                          ...
                  'autolabel',  true,                                        ...
                      'debug',  false,                                       ...
                       'host',  HostConf.host_server,                        ...
                      'local',  false,                                       ...
                  'overwrite',  false                                        ...
);%-----------------------------------------------------------------------------

[trialName, mazeName, offsets, dropSyncInd, includeSyncInd,...
 autolabel, debug, host,local, overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');


% MAIN -------------------------------------------------------------------------
% Try to grab a session list if one exists otherwise pass on the 

if ischar(Session)
    Sessions = get_session_list_v2(Session);
elseif isstruct(Session)
    Sessions = Session;
end

if isstruct(Sessions),
    for s = Sessions,
        if local,
            MTAstartup(s.project,host,HostConf.data_server);
        else
            MTAstartup(s.project,host,s.hostServer);
        end
        
        Session = MTASession.validate(s);
        
        xsync = Session.xyz.sync.copy;
        xsync = xsync + s.offsets;
        if isfield(s,'includeSyncInd'),
            if ~isempty(s.includeSyncInd)
                dropSyncInd = ~ismember(1:xsync.size(1),s.includeSyncInd);
            end
        end
        
        xsync.data(dropSyncInd,:) = [];
        Trial = MTATrial(Session,...
                         s.mazeName,...
                         s.trialName,...                         
                         overwrite,...
                        xsync);
        Trial.save();
    end
elseif (isa(Sessions,'MTASession')&&~isa(Sessions,'MTATrial'))||ischar(Sessions),
    Session = Sessions;
    if ischar(Session),
        Session = MTASession(Session,mazeName); 
    end
    xsync = Session.xyz.sync.copy;
    xsync = xsync+offsets;
    if ~isempty(includeSyncInd)
        dropSyncInd = ~ismember(1:xsync.size(1),includeSyncInd);
    end
    xsync.data(dropSyncInd,:) = [];
    Trial = MTATrial(Session,Session.maze.name,trialName,true,xsync);
    Trial.save;
end



assert(offsets(:,1)>=0&offsets(:,2)<=0,'MTA:utilities:QuickTrialSetup:offsets, see help QuickTrialSetup for offsets specifications');

% Run labelBhv if all required markers are present
% labelBhv functions only on Sessions with the H5B4(H0B9) model
rmarkers = {'spine_lower','pelvis_root', 'spine_middle', 'spine_upper',...
    'head_back',  'head_left',   'head_front',  'head_right'};

if sum(ismember(Trial.model.ml,rmarkers))==numel(rmarkers)&&autolabel
    Trial.stc.updateMode('auto_wbhr');
    Trial = labelBhv(Trial,Trial.stc);
    %Trial = labelAuxBhv(Trial,Trial.stc);
    Trial.save;
end

if nargout==1,varargout{1} = Trial;end

if debug,
    Trial.load('xyz');
    Session.load('xyz');
    figure, 
    subplot(1,2,1);
    plot(Session.xyz(round(offsets(1)*Session.xyz.sampleRate):end,'head_front',3));
    hold on,plot(Trial.xyz(:,'head_front',3),'r');
    Lines(Trial.stc{'r'}(:),[],'m');
    xlabel('samples');
    ylabel('Height');
    legend('Session','Trial','Location','NorthEastOutside');
    title([Trial.filebase , 'Rearing Events']);
    subplot(1,2,2);
    plot(Trial.xyz(:,7,1),Trial.xyz(:,7,2),'.');
end

% END MAIN -------------------------------------------------------------------------