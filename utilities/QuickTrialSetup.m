function varargout = QuickTrialSetup(Session,varargin)
%function QuickTrialSetup(Session,varargin)
%
% Variables:
%   trialName:      string, defarg  - 'all' 
%                   Name of the new Trial
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
%   autolabel:          Logical, defarg  - true
%                   Try to do some automatic behavior labeling
%
%   debug:          Logical, defarg  - false 
%                   Display some diagnostic plots
%
%

[trialName,offsets,dropSyncInd,autolabel,debug] = DefaultArgs(varargin,{'all',[0,0],[],true,false});
xsync = Session.xyz.sync.copy;
xsync = xsync+offsets;
xsync.data(dropSyncInd,:) = [];
Trial = MTATrial(Session,trialName,[],true,xsync);
Trial.save;

assert(offsets(:,1)>=0&offsets(:,2)<=0,'MTA:utilities:QuickTrialSetup:offsets, see help QuickTrialSetup for offsets specifications');

%% Run labelBhv if all required markers are present
% labelBhv functions only on Sessions with the H5B4(H0B9) model
rmarkers = {'spine_lower','pelvis_root', 'spine_middle', 'spine_upper',...
              'head_back',  'head_left',   'head_front',  'head_right',...
                'head_top'};

if sum(ismember(Trial.model.ml,rmarkers))==numel(rmarkers)&&autolabel
    Trial.stc.updateMode('auto_wbhr');
    Trial = labelBhv(Trial,Trial.stc);
    Trial = labelAuxBhv(Trial,Trial.stc);
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



