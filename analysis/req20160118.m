

Trial = MTATrial('jg05-20120317');
clear('mod');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev2_jg';
mod.featureSet = 'fet_tsne_rev10';
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'whole_state_bootstrap';
argin = struct2varargin(mod);
[stc,d_state,ls,lsm,model] = bhv_nn_multi_patternnet(Trial,argin{:});


Trial = MTATrial('Ed03-20140625');
clear('mod');
mod.states     = {'walk','rear','turn','pause','groom','sit'};
mod.stcMode    = 'hand_labeled_rev1_Ed';
mod.featureSet = 'fet_tsne_rev5';
mod.sampleRate = 12;
mod.nNeurons   = 100;
mod.randomizationMethod = 'whole_state_bootstrap';
argin = struct2varargin(mod);
[stc,d_state,ls,lsm,model] = bhv_nn_multi_patternnet(Trial,argin{:});




% Ed 
slist =     {'hand_labeled_Ed'; 'hand_labeled_jg'; 'hand_labeled_Ed'; 'hand_labeled_jg'};
refTrial =  {  'Ed03-20140625';   'jg05-20120317';   'jg05-20120317';   'Ed03-20140625'};
fetSet  = 'fet_tsne_rev5';


for sli = 1:numel(slist),
    for rti = 1:numel(refTrial),
        SesList = SessionList(slist{sli});
        states  = {'walk','rear','turn','pause','groom','sit'};
        model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];

        stc = {}; d_state = {}; ls = {}; lsm = {};
        for s = SesList
            Trial = MTATrial(s.sessionName,s.trialName,s.mazeName);
            clear('mod');
            mod.states     = states;
            mod.stcMode    = s.stcMode;
            mod.featureSet = fetSet;
            mod.model      = model;
            mod.sampleRate = 12;
            mod.nNeurons   = 100;
            argin = struct2varargin(mod);
            [stc{end+1},d_state{end+1},ls{end+1},lsm{end+1}] = bhv_nn_multi_patternnet(Trial,argin{:});
            stc{end}.save(1);
        end

        save(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']),...
             '-v7.3','slist','model','fetSet','states','stc','d_state','ls');

    end
end

%Plot Over all accuracies 


prop = 'accuracy';
figure,plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'.')
xlim([0,5])
ylim([50,100])

prop = 'sensitivity';
prop = 'precision';
figure,plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                     repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'.')
xlim([0,7])
hax = gca;
hax.XTickLabelMode = 'manual';
hax.XTickLabel = cat(2,{''},states,{''});
ylabel(prop)


figure,hold on,
for sli = 1:numel(slist),
    for rti = 1:numel(refTrial),    
        prop = 'accuracy';
        load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));
        ,plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'.')
    end
end