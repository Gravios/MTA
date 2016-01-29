
%Trian bhv_nn_multi_patternnet and plot 

rlist = SessionList('training_hand_labeled');
fetSet  = 'fet_tsne_rev12';
sampleRate = 12;
nNN = 10;
states = {'walk','rear','turn','pause','groom','sit'};
rndMethod = 'WSB';

for s = rlist'
    clear('mod');
    mod.states     = states;
    mod.stcMode    = s.stcMode;
    mod.featureSet = fetSet;
    mod.sampleRate = sampleRate;
    mod.nNeurons   = nNN;
    mod.randomizationMethod = rndMethod
    argin = struct2varargin(mod);
    bhv_nn_multi_patternnet(MTATrial(s),argin{:});
end



% Ed 
slist =     {     'hand_labeled_Ed';      'hand_labeled_jg'};
for sli = 1:numel(slist),
    for rti = 1:numel(refTrial),
        SesList = SessionList(slist{sli});
        model = ['MTAC_BATCH-' fetSet '_SR_' num2str(sampleRate) '_REF_' rlist(rti).name, ...
                 '.cof.all_STC_' rlist(rti).stcMode '_NN_' num2str(nNN) '_NN_multiPN_RAND_'...
                rndMethod];

        stc = {}; d_state = {}; ls = {}; lsm = {};
        for s = SesList
            Trial = MTATrial(s.sessionName,s.trialName,s.mazeName);
            Trial.load('stc',s.stcMode);
            clear('mod');
            mod.states     = states;
            mod.stcMode    = rlist(rti).stcMode;
            mod.featureSet = fetSet;
            mod.model      = model;
            mod.sampleRate = sampleRate;
            mod.nNeurons   = nNN;
            mod.map2reference = true;
            argin = struct2varargin(mod);
            [stc{end+1},d_state{end+1},ls{end+1},lsm{end+1}] = bhv_nn_multi_patternnet(Trial,argin{:});
            stc{end}.save(1);
        end

        if mod.map2reference,
            mapped = '-map2ref';
        else
            mapped = '';
        end
        
        save(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,mapped,'.mat']),...
             '-v7.3','slist','rlist','nNN','sampleRate','model','fetSet','rndMethod',...
                     'states','stc','d_state','ls');

    end
end
%Plot Over all accuracies 


% Visuallize the Accuracy, Precision, and Sensitivity of model with
% respect to a second set of labels (Usually hand labeled)
slist =     {     'hand_labeled_Ed';      'hand_labeled_jg'};
refTrial =  {       'Ed03-20140625';        'jg05-20120317'};
trnStcMode ={'hand_labeled_rev1_Ed', 'hand_labeled_rev2_jg'};

sli = 2;
rti = 2;
fetSet  = 'fet_tsne_rev11';
SesList = SessionList(slist{sli});
SesList = {SesList(:).sessionName};
model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];
load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));

prop = 'accuracy';
figure,plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
xlim([0,5])
ylim([20,100])
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
ylim([20,100])
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
ylim([20,100])
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



slist =     {     'hand_labeled_Ed';      'hand_labeled_jg'};
refTrial =  {       'Ed03-20140625';        'jg05-20120317'};
trnStcMode ={'hand_labeled_rev1_Ed', 'hand_labeled_rev2_jg'};

fetSet = 'fet_tsne_rev12';
sli = 1;SesList = SessionList(slist{sli});
states  = {'walk','rear','turn','pause','groom','sit'};
stc = {}; d_state = {}; ls = {}; lsm = {};

s = SesList(2);
for rti = 1:2,

    model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} ...
             '.cof.all_STC_' trnStcMode{rti} '_NN_100_NN_multiPN_RAND_WSB'];

    Trial = MTATrial(s.sessionName,s.trialName,s.mazeName);
    Trial.load('stc',s.stcMode);
    clear('mod');
    mod.states     = states;
    mod.stcMode    = trnStcMode{rti};
    mod.featureSet = fetSet;
    mod.model      = model;
    mod.sampleRate = 12;
    mod.nNeurons   = 100;
    mod.map2reference = false;
    argin = struct2varargin(mod);
    [~,d_state{end+1},ls{end+1}] = bhv_nn_multi_patternnet(Trial,argin{:});
end
 



sli = 1;
rti = 1;
states  = {'walk','rear','turn','pause','groom','sit'};
fetSet  = 'fet_tsne_rev10';
SesList = SessionList(slist{sli});
SesList = {SesList(:).sessionName};
model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];
ds = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));

sli = 1;
rti = 2;
fetSet  = 'fet_tsne_rev10';
SesList = SessionList(slist{sli});

model = ['MTAC_BATCH-' fetSet '_SR_12_REF_' refTrial{rti} '.cof.all_NN_100_NN_multiPN_RAND_WSB'];
os = load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']));

Trial = MTATrial(SesList(1).sessionName,...
                 SesList(1).trialName,...
                 SesList(1).mazeName);
sts = stc2mat(Trial.load('stc',SesList(1).stcMode),Trial.load('xyz'),states);


figure,sp = [];
sp(end+1) = subplot(411),imagesc(ds.d_state{1}')
sp(end+1) = subplot(412),imagesc(os.d_state{1}')
sp(end+1) = subplot(413),imagesc(os.d_state{1}'+ds.d_state{1}')
sp(end+1) = subplot(414),imagesc(~~sts')
linkaxes(sp,'xy');

msl = (os.d_state{1}+ds.d_state{1})>190;
sum(all(~~sts==msl,2)&&~(all(msl==0|~sts,2)))

sum(any(~~sts,2))