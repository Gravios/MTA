function req20160310_1_preproc(Trial,varargin);
% 1. Preproc features
defargs = {...
    ... sampleRate
        12,...
    ...
    ...states 
       {'walk','rear','turn','pause','groom','sit'},...
    ...
    ...stcMode
       'hand_labeled_rev3_jg',...
    ...train
       false,...
    ...RefTrial
       {'jg05-20120317','all','cof'},...
    ...
    ...RefStcMode
       'hand_labeled_rev3_jg',...
    ...
    ...nNeurons 
       200,...
    ...
    ...nIter
       100,...
    ...
    ...rndMethod
       'WSBNT'...
};


nIter    = 100;
rndMethod = 'WSBNT';


[sampleRate,states,stcMode,train,RefTrial,RefStcMode] = DefaultArgs(varargin,defargs);

if train
    %Train Parm
    Trial = MTATrial.validate();
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = [];
    rMean = []; rStd = [];
else
    %Test Parm
    Trial = MTATrial.validate(Trial);    
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;

    RefTrial = MTATrial.validate(RefTrial);
    RefTrial.load('stc',RefStcMode);

    rfet = fet_all(RefTrial,sampleRate,[]);
    rfet.data = [rfet.data,rfet.data.^2];
    rafet = rfet.copy;
    for sh = 1:rfet.size(2)-1;
        rfet.data = [rfet.data,circshift(rafet.data',-sh)'.*rafet.data];
    end
    [~,raMean,raStd] = unity(rafet);
    [~,rMean,rStd] = unity(rfet);
    clear('rfet')
end

% LOAD all features
try
    afet = Trial.fet.copy;
    afet.label = 'fet_all';
    afet.load;
catch
    afet = fet_all(Trial,sampleRate,RefTrial);
    afet.save;
end

try
    fet = Trial.fet.copy;
    fet.label = 'fet_all_quad';
    fet.load;
catch
    fet.data = [fet.data,fet.data.^2];
    afet = fet.copy;
    for sh = 1:fet.size(2)-1;
        fet.data = [fet.data,circshift(afet.data',-sh)'.*afet.data];
    end
    fet.save;
end


% NORMALIZE features
if ~isempty(rMean)&&~isempty(rStd),
    afet = unity(afet,[],raMean,raStd);
    fet = unity(fet,[],rMean,rStd);
else
    afet = afet.unity;
    fet  =  fet.unity;
end

if train
    [tstc,~,tfet] = resample_whole_state_bootstrap_noisy_trim(stc,fet,states);
    tstc.states{end}.data = [1,tfet.size(1)];    
    [stateOrd,fetInds,miAstates] = select_features_hmi(Trial,tstc,tfet,states,false);
    save(fullfile(Trial.path.data,'analysis','req20160310_1_preproc.mat'),...
         'stateOrd','fetInds','miAstates');
else
    load(fullfile(Trial.path.data,'analysis','req20160310_1_preproc.mat'));
end

save(fullfile(Trial.spath,'req20160310_1_preproc.mat'));