function QuickTrialSetup(Session,varargin)
[trialName,offsets,dropSyncInd] = DefaultArgs(varargin,{'all',[8,0],[]});
xsync = Session.xyz.sync.copy;
xsync = xsync+offsets;
xsync.data(dropSyncInd,:) = [];
xsync.resample(1);
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
