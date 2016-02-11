
%Trian bhv_nn_multi_patternnet and plot 
rlist = SessionList('training_hand_labeled');
slist = {'hand_labeled_jg';'hand_labeled_Ed'};
fetSet  = 'fet_tsne_rev15';
sampleRate = 12;
nNeurons = 100;
nIter = 100;
states = {'walk','rear','turn','pause','groom','sit'};
rndMethod = 'WSBNT';
norm = true;
mref = true;


for s = rlist
    clear('mod');
    mod.states     = states;
    mod.stcMode    = s.stcMode;
    mod.featureSet = fetSet;
    mod.sampleRate = sampleRate;
    mod.nNeurons   = nNeurons;
    mod.nIter      = nIter;
    mod.randomizationMethod = rndMethod;
    mod.normalize = norm;
    argin = struct2varargin(mod);
    bhv_nn_multi_patternnet(MTATrial.validate(s),argin{:});
end



% Ed 
for sli = 1:numel(slist),
    for rti = 1:numel(rlist),
        SesList = SessionList(slist{sli});
        model = ['MTAC_BATCH-' fetSet ...
                 '_SR_' num2str(sampleRate) ...
                 '_NORM_' num2str(norm) ...
                 '_REF_' rlist(rti).sessionName, '.' rlist(rti).mazeName '.' rlist(rti).trialName ...
                 '_STC_' rlist(rti).stcMode ...
                 '_NN_' num2str(nNeurons) ...
                 '_NI_' num2str(nIter) ...
                 '_NN_multiPN_RAND_' rndMethod];

        stc = {}; d_state = {};p_state = {}; ls = {}; lsm = {};mdl = {};
        for s = SesList
            Trial = MTATrial(s.sessionName,s.trialName,s.mazeName);
            Trial.load('stc',s.stcMode);
            clear('mod');
            mod.states     = states;
            mod.stcMode    = s.stcMode;
            mod.featureSet = fetSet;
            mod.model      = model;
            mod.sampleRate = sampleRate;
            mod.nNeurons   = nNeurons;
            mod.nIter      = nIter;
            mod.map2reference = mref;
            mod.normalize = norm;
            argin = struct2varargin(mod);
            [stc{end+1},d_state{end+1},ls{end+1},lsm{end+1},mdl{end+1},p_state{end+1}] = bhv_nn_multi_patternnet(Trial,argin{:});
            stc{end}.save(1);
        end

        if mod.map2reference,
            mapped = '-map2ref';
        else
            mapped = '';
        end
        
        save(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,mapped,'.mat']),...
             '-v7.3','slist','rlist','nNeurons','nIter','sampleRate','model','fetSet','rndMethod',...
                     'states','stc','d_state','p_state','ls','mdl');

    end
end
%Plot Over all accuracies 


% Visuallize the Accuracy, Precision, and Sensitivity of model with
% respect to a second set of labels (Usually hand labeled)


sli = 1;
rti = 1;
SesList = SessionList(slist{sli});
SesList = {SesList(:).sessionName};

model = ['MTAC_BATCH-' fetSet ...
         '_SR_' num2str(sampleRate) ...
         '_NORM_' num2str(norm) ...         
         '_REF_' rlist(rti).sessionName, '.' rlist(rti).mazeName '.' rlist(rti).trialName ...
         '_STC_' rlist(rti).stcMode ...
         '_NN_' num2str(nNeurons) ...
         '_NI_' num2str(nIter) ...         
         '_NN_multiPN_RAND_' rndMethod];

load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,mapped,'.mat']))
slist = {'hand_labeled_jg';'hand_labeled_Ed'};


%figure
hfig = figure(3923992),clf
set(hfig,'Position',       [40         40        1500         500])

prop = 'accuracy';
subplot(131);
plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
xlim([0,5])
ylim([0,100])
ylabel(prop);
title({['Training Set: ' rlist(rti).sessionName],...
       ['Randomization: ' rndMethod ],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',1:numel(SesList));
set(gca,'XTickLabel',SesList);
set(gca,'XTickLabelRotation',90);
pause(.1)

prop = 'precision';
subplot(132);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',4,6)','d-')
xlim([0,7])
hax = gca;
hax.XTickLabelMode = 'manual';
hax.XTickLabel = cat(2,{''},states,{''});
ylim([0,100])
ylabel(prop)
title({['Training Set: ' rlist(rti).sessionName],...
       ['Randomization: ' rndMethod ],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
legend(SesList,'location','southwest');
pause(.1)

prop = 'sensitivity';
subplot(133);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                     repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',6,4),'d-')
xlim([0,7])
hax = gca;
hax.XTickLabelMode = 'manual';
hax.XTickLabel = cat(2,{''},states,{''});
ylim([0,100])
ylabel(prop)
title({['Training Set: ' rlist(rti).sessionName],...
       ['Randomization: ' rndMethod ],...
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
legend(SesList,'location','southwest');
pause(.1)






figure,hold on,
fetSet  = 'fet_tsne_rev13';
k = 1;
htl = {};
set(gcf,'interpreter','none')
for sli = 1:numel(slist),
    for rti = 1:numel(refTrial),    
        prop = 'accuracy';

        model = ['MTAC_BATCH-' fetSet ...
                 '_SR_' num2str(sampleRate) ...
                 '_NORM_' num2str(norm) ...
                 '_REF_' rlist(rti).sessionName, '.' ...
                         rlist(rti).mazeName '.' ...
                         rlist(rti).trialName ...
                 '_STC_' rlist(rti).stcMode ...
                 '_NN_' num2str(nNeurons) ...
                 '_NI_' num2str(nIter) ...
                 '_NN_multiPN_RAND_' rndMethod];

        ds = load(fullfile(MTASession().path.data,'analysis',...
                           [slist{sli},'-',model,mapped,'.mat']));

        scat = scatter(k*ones([numel(ds.ls),1]),cell2mat(cellfun(@subsref,ds.ls, ...
                             repmat({substruct('.',prop)},[1,numel(ds.ls)]),'uniformoutput',false))'.*100,10);

        scat.MarkerFaceColor = scat.CData;
        k=k+1;
        htl = cat(2,htl,{['{' slist{sli} ' - ' rlist(rti).sessionName '}']});
    end
end
pp
xlim([0.5,k-0.5]);
ylim([50,100])

set(gca,'XTickLabelMode','manual');
set(gca,'XTick',1:k-1);
set(gca,'XTickLabel',{});


set(gca,'XTickLabel',htl);
set(gca,'XTickLabelRotation',90);






