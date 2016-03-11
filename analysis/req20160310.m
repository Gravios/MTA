

train = true;
states = {'walk','rear','turn','pause','groom','sit'};
target = 'sit';
sampleRate = 12;
gStates = states;


if train
    %Train Parm
    stcMode = 'hand_labeled_rev3_jg';
    Trial = MTATrial.validate('jg05-20120317');
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = [];
    rMean = []; rStd = [];
else
    %Test Parm
    %Trial = MTATrial.validate('Ed03-20140625');
    stcMode = 'hand_labeled_rev1_Ed';
    Trial = MTATrial.validate({'Ed05-20140529','all','ont'});    
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = MTATrial.validate('jg05-20120317');
    RefTrial.load('stc','hand_labeled_rev3_jg');
    rfet = fet_all(RefTrial,sampleRate,[]);
    rfet.data = [rfet.data,rfet.data.^2];
    rafet = rfet.copy;
    for sh = 1:rfet.size(2)-1;
        rfet.data = [rfet.data,circshift(rafet.data',-sh)'.*rafet.data];
    end
    [~,rMean,rStd] = unity(rfet);
    clear('rfet')
end


% LOAD all features
fet = fet_all(Trial,sampleRate,RefTrial);
fet.data = [fet.data,fet.data.^2];
afet = fet.copy;
for sh = 1:fet.size(2)-1;
    fet.data = [fet.data,circshift(afet.data',-sh)'.*afet.data];
end

% NORMALIZE features
if ~isempty(rMean)&&~isempty(rStd),
    fet = unity(fet,[],rMean,rStd);
else
    fet = fet.unity;
end


if train,
    [tstc,~,tfet] = resample_whole_state_bootstrap_noisy_trim(stc,fet,states);
    [stateOrd,fetInds] = select_features_hmi(Trial,tstc,tfet,states,true);
    mpn = struct;
    for s = stateOrd,
        tstates = {[strjoin({gStates{find(cellfun(@isempty,regexp(gStates,gStates{sind})))}},'+'),'&gper'],...
                   gStates{sind}})
        
        [tsnem] = mta_tsne(Trial,tfet,12,tstc,tstates,5,2,80);

        sfet = fet.copy;
        sfet.data = tfet(:,fetInds{s});

        [mpn.Stc,mpn.d_state,mpn.labelingStats, ...
         mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
            bhv_nn_muli_patternnet(Trial,states,stcMode,sfet,[],[],200,100,'WSBNT');
        
        gStates(stateOrd)=[];
    end
end


%'tsnem','mpn'