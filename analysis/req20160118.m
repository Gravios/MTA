

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
slist =     {'hand_labeled_Ed'; 'hand_labeled_jg'};
refTrial =  {  'Ed03-20140625';   'jg05-20120317'};
fetSet  = 'fet_tsne_rev10';


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

sli = 1;
rti = 1;
fetSet  = 'fet_tsne_rev5';
SesList = SessionList(slist{sli});
SesList = {SesList(:).sessionName};
model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];
load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));

prop = 'accuracy';
figure,plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
xlim([0,5])
ylim([50,100])
ylabel(prop);
title({['Training Set: ' refTrial{rti}],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',1:numel(SesList));
set(gca,'XTickLabel',SesList);
set(gca,'XTickLabelRotation',90);
pause(.1)

prop = 'precision';
figure,plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'d-')
xlim([0,7])
hax = gca;
hax.XTickLabelMode = 'manual';
hax.XTickLabel = cat(2,{''},states,{''});
ylabel(prop)
title({['Training Set: ' refTrial{rti}],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
legend(SesList)
pause(.1)




prop = 'sensitivity';
figure,plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                     repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'d-')
xlim([0,7])
hax = gca;
hax.XTickLabelMode = 'manual';
hax.XTickLabel = cat(2,{''},states,{''});
ylabel(prop)
title({['Training Set: ' refTrial{rti}],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
legend(SesList)
pause(.1)


figure,hold on,
slist =     {'hand_labeled_Ed'; 'hand_labeled_jg'};
refTrial =  {  'Ed03-20140625';   'jg05-20120317'};
fetSet  = 'fet_tsne_rev5';
k = 1;
htl = {};
set(gcf,'interpreter','none')
for sli = 1:numel(slist),
    for rti = 1:numel(refTrial),    
        prop = 'accuracy';
        model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];
        ds = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));
        scat = scatter(k*ones([numel(ds.ls),1]),cell2mat(cellfun(@subsref,ds.ls, ...
                             repmat({substruct('.',prop)},[1,numel(ds.ls)]),'uniformoutput',false))'.*100,10);

        scat.MarkerFaceColor = scat.CData;
        k=k+1;
        htl = cat(2,htl,{['{' slist{sli} ' - ' refTrial{rti} '}']});
    end
end

xlim([0.5,k-0.5]);
ylim([50,100])

set(gca,'XTickLabelMode','manual');
set(gca,'XTick',1:k-1);
set(gca,'XTickLabel',{});


set(gca,'XTickLabel',htl);
set(gca,'XTickLabelRotation',90);

%Plot Over all accuracies 

