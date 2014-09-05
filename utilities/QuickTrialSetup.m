function Trial = QuickTrialSetup(Session,varargin)
%function QuickTrialSetup(Session,varargin)
%[trialName,offsets,dropSyncInd,debug] = DefaultArgs(varargin,{'all',[0,0],[],false});

[trialName,offsets,dropSyncInd,debug] = DefaultArgs(varargin,{'all',[0,0],[],false});
xsync = Session.xyz.sync.copy;
xsync = xsync+offsets;
xsync.data(dropSyncInd,:) = [];
Trial = MTATrial(Session,trialName,[],true,xsync);
Trial.save;

%% Run labelBhv if all required markers are present
% labelBhv functions only on Sessions with the H5B4(H0B9) model
rmarkers = {'spine_lower','pelvis_root', 'spine_middle', 'spine_upper',...
              'head_back',  'head_left',   'head_front',  'head_right',...
                'head_top'};

if sum(ismember(Trial.model.ml,rmarkers))==numel(rmarkers)
    Trial.stc.updateMode('auto_wbhr');
    Trial = labelBhv(Trial,Trial.stc);
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



