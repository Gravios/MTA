function req20160310_1_preproc(Trial,varargin);
%function req20160310_1_preproc(Trial,varargin);
% 1. Preprocess features and select state order
%
% varargin:
%
%  sampleRate,12,                                            
%  states,    {{'walk','rear','turn','pause','groom','sit'}},
%  stcMode,   'hand_labeled_rev3_jg',                        
%  train,     true,                                          
%  RefTrial,  'jg05-20120317.all.cof',                       
%  RefStcMode,'hand_labeled_rev3_jg',                        
%  nNeurons,  200,                                           
%  nIter,     100,                                           
%  rndMethod, 'WSBNT'                                        
%  
% Out:
%  None
%
%

% Default Arguments
defargs = ...
  struct('sampleRate',12,                                            ...
         'states',    {{'walk','rear','turn','pause','groom','sit'}},...
         'stcMode',   'hand_labeled_rev3_jg',                        ...
         'train',     true,                                          ...
         'RefTrial',  'jg05-20120317.all.cof',                       ...
         'RefStcMode','hand_labeled_rev3_jg',                        ...
         'nNeurons',  200,                                           ...
         'nIter',     100,                                           ...
         'rndMethod', 'WSBNT'                                        ...
);

[sampleRate,states,stcMode,train,RefTrial,...
 RefStcMode,nNeurons,nIter,rndMethod] = DefaultArgs(varargin,defargs);

if train
    %Train Parm
    Trial = MTATrial.validate(Trial);
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = [];
    rMean = []; rStd = [];
    fet = fet_all(Trial,sampleRate);        
else
    %Test Parm
    Trial = MTATrial.validate(Trial);    
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;

    RefTrial = MTATrial.validate(RefTrial);
    RefTrial.load('stc',RefStcMode);

    rfet = fet_all(RefTrial,sampleRate,[]);
    rafet = rfet.copy;
    rfet.data = [rfet.data,rfet.data.^2];
    for sh = 1:rfet.size(2)-1;
        rfet.data = [rfet.data,circshift(rafet.data',-sh)'.*rafet.data];
    end
    [~,raMean,raStd] = unity(rafet);
    [~,rMean,rStd] = unity(rfet);
    clear('rfet')
    fet = fet_all(Trial,sampleRate,RefTrial);    
end

% LOAD all features


afet = fet.copy;
fet.data = [afet.data,afet.data.^2];
qfet = fet.copy;
for sh = 1:qfet.size(2)-1;
    fet.data = [fet.data,circshift(qfet.data',-sh)'.*qfet.data];
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

save(fullfile(Trial.spath,'req20160310_1_preproc-tfet.mat'),...
     'states','fetInds','stateOrd','tfet','tstc','fet','-v7.3');
save(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'),...
     'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod','-v7.3');
