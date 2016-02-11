
%% Model Parameters
rlist = SessionList('training_hand_labeled');
slist = {'hand_labeled_jg';'hand_labeled_Ed'};
fetSet  = 'fet_tsne_rev13';
sampleRate = 12;
nNeurons = 100;
states = {'walk','rear','turn','pause','groom','sit'};
norm = true;
mref = true;


%% Train Neural Network
for  s = rlist
    Trial = MTATrial.validate(s);
    Trial.load('stc',s.stcMode);
    states = {'walk','rear','turn','pause','groom','sit'};
    featureSet = 'fet_tsne_rev13';
    sampleRate = 12;
    nNeurons = 100;
    model = ['MTAC_' featureSet ...
             '_SR_'  num2str(sampleRate) ...
             '_NORM_'  num2str(norm)...
             '_REF_' Trial.name ...
             '_STC_' s.stcMode ...
             '_NN_'  num2str(nNeurons) ];

    features = feval(featureSet,Trial,sampleRate);
    if norm, [~,rm,rs] = unity(features); end
    %features.data =  features.data./5;

    bhv_nn (Trial,true,states,features,model,'nNeurons',nNeurons);
end



%% Label and Compare to Hand Labeled Data
for rti = 1:numel(rlist),
    refTrial = MTATrial.validate(rlist(rti));
    refTrial.load('stc',rlist(rti).stcMode);
    refFet = feval(featureSet,refTrial,sampleRate);
    if norm, [~,rm,rs] = unity(refFet); end

    for sli = 1:numel(slist),

        SesList = SessionList(slist{sli});

        model = ['MTAC_'  featureSet ...
                 '_SR_'   num2str(sampleRate) ...
                 '_NORM_' num2str(norm)...
                 '_REF_'  refTrial.filebase ...
                 '_STC_'  refTrial.stc.mode ...
                 '_NN_'   num2str(nNeurons) ];
        ls = {};Stc = {};d_state = {};
        for s = SesList,
            Trial = MTATrial.validate(s);
            Trial.load('stc',s.stcMode);
            StcHL = Trial.stc.copy;

            xyz = Trial.load('xyz');
            features = feval(featureSet,Trial,sampleRate);
            features.map_to_reference_session(Trial,refTrial);
            if norm, features.unity([],rm,rs); end        
            %features.data =  features.data./5;
            % Create state matrix (N x k) N=samples, k=states, Domain:Boolean
            [Stc{end+1},d_state{end+1}] = bhv_nn (Trial,false,states,features,model);
            shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
            ysm = MTADxyz('data',double(0<stc2mat(Stc{end},  xyz,states)),'sampleRate',xyz.sampleRate); 

            aper = Trial.stc{'a'}.cast('TimeSeries');
            aper.resample(ysm);
            ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

            tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlabe
            ls{end+1}.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
            ls{end}.precision = round(diag(tcm)./sum(tcm,2),4).*100;
            ls{end}.confusionMatrix = round(tcm./xyz.sampleRate,2);
            ls{end}.accuracy = round(sum(diag(tcm))/sum(tcm(:)),4);
            ls{end}.acm = cat(2,ls{end}.confusionMatrix,ls{end}.precision);
            ls{end}.acm = cat(1,ls{end}.acm,[ls{end}.sensitivity,ls{end}.accuracy]);


        end
        save(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']),...
             '-v7.3','slist','rlist','nNeurons','sampleRate','model','fetSet',...
             'states','Stc','d_state','ls');
    
    end

end




sli = 2;
rti = 2;
SesList = SessionList(slist{sli});
SesList = {SesList(:).sessionName};

model = ['MTAC_'  featureSet ...
         '_SR_'   num2str(sampleRate) ...
         '_NORM_' num2str(norm)...
         '_REF_'  rlist(rti).sessionName, '.' rlist(rti).mazeName '.' rlist(rti).trialName ...
         '_STC_'  rlist(rti).stcMode ...
         '_NN_'   num2str(nNeurons) ];

load(fullfile(MTASession().path.data,'analysis',[slist{sli},'-',model,'.mat']))


%figure
hfig = figure(23093020),clf
set(hfig,'Position',       [40         40        1500         500])
prop = 'accuracy';
subplot(131);
plot(cell2mat(cellfun(@subsref,ls, ...
                             repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
xlim([0,5])
ylim([20,100])
ylabel(prop);
title({['Training Set: ' rlist(rti).sessionName],...
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
       ['Labeling Set: ' slist{sli}],...
       ['Feature  Set: ' fetSet]});
legend(SesList,'location','southwest');
pause(.1)



saveas('')













rcomp = Trial.stc{'r'}&Stc{4}{'r'};

figure,plot(xyz(:,7,3))
Lines(Trial.stc{'r'}(:),[],'r');
Lines(Stc{4}{'r'}(:),[],'g');

                    
                    
                    
                    
                    
%% New heuristics
td = [];
for s = Stc{'w'}.data',
    td(end+1) = sqrt(sum((xyz(s(2),1,:)-xyz(s(1),1,:)).^2,3));
end                    



td = [];
for s = Stc{'r'}.data',
    td(end+1) = median(xyz(s(1):s(2),5,3));
end                    
th = [];
for s = StcHL{'r'}.data',
    th(end+1) = median(xyz(s(1):s(2),5,3));
end                    
figure,hold on,
plot(diff(StcHL{'r'}.data,1,2)/120,th','.b')
plot(diff(Stc{'r'}.data,1,2)/120,td','.r')


z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rl =  z.segs(Stc{'r'}(:,1)-60,120,nan);

z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rh =  z.segs(StcHL{'r'}(:,1)-60,120,nan);

z = Trial.xyz.copy;
z.data = xyz(:,5,3);
rh =  z.segs(StcHL{'r'}(:,1)-120,240,nan);


c = convn(mean(rh(90:150,:),2)',rh(:,2)','full');
figure,plot(c)
figure,plot(diff(c))

figure,sp = [];
sp(end+1) = subplot(211);
imagesc(shl.data')
sp(end+1) = subplot(212);
imagesc(ysm.data')
linkaxes(sp,'xy');

s = 6;
xlim(Stc{'r'}(s,:)+[-60,60])
StcHL{'r'}(s,:)
     


f = 4
figure,hold on
eds = linspace(0,170,100);
hs = bar(eds,histc(tfet(:,f),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='r';
hs = bar(eds,histc(trfet(:,f),eds),'histc');hs.FaceAlpha=.5;hs.FaceColor='c';

Trial = MTATrial('jg05-20120317');
tfet = fet_tsne(Trial,20,false);
ofet = fet_tsne(Trial,20,false);
features = fet20151007(Trial,20,false);

figure,hold on
plot(tfet(:,1))
plot([ofet(:,1),features(:,1)])



%% Stop here this hnn is shit
%% Heirarchical nn labeling 
Trial = MTATrial('jg05-20120317');

Trial.load('stc','hand_labeled_rev2');
features = fet_tsne(Trial,20,false);

nNeurons = 10;
states = {'rear'};
model = 'MTAC_fet_tsne_REFjg0520120317_NN_rear';



% Train Model
      bhv_nn (Trial,true,states,features,model,false,true);

Trial = MTATrial('Ed01-20140707');
StcHL = Trial.load('stc','hand_labeled_rev1');
%features = fet20151007(Trial,20,false);
features = fet_tsne(Trial,20,false);

% Label Trial
Stc = bhv_nn (Trial,false,states,features,model,false,true);


xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,{states{1},['gper-' states{1}]})),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,{states{1},'other'})),'sampleRate',xyz.sampleRate); 

ind = any(shl.data,2)&any(ysm.data,2);

tcm = confmat(shl(ind,:),ysm(ind,:)); % DEP: netlab
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
precision = round(diag(tcm)./sum(tcm,2),4).*100;
cm = round(tcm./xyz.sampleRate,2);



xyz = Trial.load('xyz');

hvar = xyz.copy;
hvar.filter('ButFilter',3,2.4,'low');
hvar = hvar.vel(1,[1,2]);

hvar = features;
hind = 6;

%hind = 1;
%%
eds = linspace(-3,2,100);
figure, hold on
ind = Trial.stc{'a-w-n'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'c';
hs.FaceAlpha = .5;
ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(hvar(ind,hind)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .5;




f = 2; s ='n' ;
figure,
eds = 50:.5:120;
eds = -3:.05:2;
subplot(311),
bar(eds,hist(rfet(rStc{s},f),eds),'histc'),Lines(mean(rfet(rStc{s},f)),[],'r')
subplot(312),
bar(eds,hist(features(StcHL{s},f),eds),'histc'),Lines(mean(features(StcHL{s},f)),[],'r')
subplot(313),
bar(eds,hist(nfet(StcHL{s},f),eds),'histc'),Lines(mean(nfet(StcHL{s},f)),[],'r')
