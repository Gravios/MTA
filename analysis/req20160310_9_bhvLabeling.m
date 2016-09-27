function req20160310_9_bhvLabeling(Trial)


rlist = get_session_list('training_hand_labeled');
slist = {'hand_labeled_jg';'hand_labeled_Ed'};



Trial = MTATrial.validate(Trial);
RefTrial = MTATrial.validate('jg05-20120317.cof.all');
RefTrial.load('stc','hand_labeled_rev3_jg');

%'states','fetInds','stateOrd','afet','nNeurons','nIter','rndMethod'
ds = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
afet = ds.afet;
ds = load(fullfile(RefTrial.spath,'req20160310_1_preproc-afet.mat'));

bs = load(fullfile(RefTrial.spath,'req20160310_5_genfigs.mat'));
ts = load(fullfile(RefTrial.spath,'req20160310_8_genOptfigs.mat'));

opn = struct;

for s = 1:5,
    sbind = bs.bFetInds{s};
    gStates = ds.states(cellfun(@isempty,...
                         regexp(ds.states,...
                                ['(^',strjoin(ds.stateOrd(1:s-1),'$)|(^'),'$)'])...
                         )...
                 );

    f = numel(ts.bfets{s});

    sub_fet = afet.copy;
    sub_fet.data = afet(:,sbind(1:f));
    sub_fet.label = [afet.label '-req20160310-optfet' ds.stateOrd{s} '-' num2str(f)];
    sub_fet.key = 'x';
    sub_fet.updateFilename(Trial);
    
    model = ['MTAC_BATCH-' ds.stateOrd{s} sub_fet.label ...
             '_SR_'  num2str(sub_fet.sampleRate) ...
             '_NORM_' num2str(0) ...             
             '_REF_' RefTrial.filebase ...
             '_STC_' RefTrial.stc.mode ...
             '_NN_'  num2str(ds.nNeurons)...
             '_NI_'  num2str(ds.nIter)...
             '_'     'NN_multiPN'...
             '_'     'RAND_' ds.rndMethod];

    [opn(s).Stc,opn(s).d_state,opn(s).labelingStats, ...
     opn(s).labelingStatsMulti,opn(s).model,opn(s).p_state] = ...
        bhv_nn_multi_patternnet(Trial,gStates,Trial.stc,sub_fet,...
                                [], model,ds.nNeurons,ds.nIter,...
                                ds.rndMethod,'targetState',ds.stateOrd{s});
end


tmat = zeros([size(opn(s).d_state,1),6]);
for s = 1:5,
    [~,sind] = max(ButFilter(opn(s).d_state,3,1.4/(Trial.xyz.sampleRate*0.5),'low'),[],2);
    %[~,sind] = max(opn(s).d_state,[],2);
    if s ~= 5,
        tmat(sind==2,s) = 7-s;
    else
        tmat(sind==2,s) = 7-s;
        tmat(sind==1,s+1) = 6-s;
    end
end
[~,sind] = max(tmat,[],2);

smat = zeros(size(tmat'));
smat([1:6:prod(size(smat))]-1+sind')=1;

smat = MTADxyz('data',smat','sampleRate',Trial.xyz.sampleRate);
ymat = MTADxyz('data',stc2mat(Trial.stc,Trial.load('xyz'),ds.stateOrd),'sampleRate',Trial.xyz.sampleRate);

ind = Trial.stc{'a'}.cast('TimeSeries');
ind = ind.data==1&any(ymat.data,2)&any(smat.data,2);

tcm = confmat(ymat(ind,:),smat(ind,:)); % DEP: netlab
confusionMatrix= round(tcm./Trial.xyz.sampleRate,2);
precision= round(diag(tcm)./sum(tcm,2),4).*100;
sensitivity= round(diag(tcm)'./sum(tcm),4).*100;
accuracy = sum(diag(tcm))/sum(tcm(:));

% Create new StateCollection ... well copy
tfb = Trial.filebase;
tfb(tfb=='.')='-';
Stc = Trial.stc.copy;
Stc.updateMode(['hmis_REF_' tfb]);
Stc.states = {};


% Populate Stc object with the new states
for i = 1:numel(ds.stateOrd),
    sts = ThreshCross(smat(:,i)==1,0.5,0);
    if ~isempty(sts),
        sts = bsxfun(@plus,sts,[1,0]);
    end
    
    Stc.addState(Trial.spath,...
             Trial.filebase,...
             sts,...
             Trial.xyz.sampleRate,...
             Trial.xyz.sync.copy,...
             Trial.xyz.origin,...
             ds.stateOrd{i},...
             [Trial.stc{ds.stateOrd{i}}.key],...
             'TimePeriods');
end



%save(fullfile(Trial.spath,[mfilename,num2str(s),'.mat']),'sbind','accum_acc','accum_pre','accum_sen','-v7.3');