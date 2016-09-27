
subjectName = 'jg05';

sessionList = get_session_list(subjectName);

Trial = MTATrial.validate('jg05-20120317');
Trial.load('stc','hl_3_jg_r');

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
ang = create(MTADang,Trial,xyz);

ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',xyz);
[State, hmm, decode] = gausshmm(ang(ind.data,'head_back','head_front',2),2);

headPitchState = zeros([size(xyz,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;
headPitchState(ind.data==1,2) = State'~=stsInd;

figure,hold on
plot(ang(headPitchState(:,1)==1,'head_back','head_front',2));
plot(headPitchState);


figure,hold on
eds = linspace(-pi/2,pi/2,200);
ind = headPitchState(:,1)==1;
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
ind = headPitchState(:,2)==1;
hs = bar(eds,histc(ang(ind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'r';
hs.EdgeColor = 'r';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;


% HMM head pitch during locomotion and pause states
angThreshHMM = mean([prctile(ang(headPitchState(:,1)==1,'head_back','head_front',2),2),prctile(ang(headPitchState(:,2)==1,'head_back','head_front',2),98)]);



% MEAN rhm ang threshold
compute_session_rhm_distribution(Trial,[],'loc+pause','hl_3_jg_r');

afig = hgload(fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig'));
rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','compute_session_rhm_distribution-hangle'));
rhm_distrb = rhm_distrb(1);
figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
mrhmp(mrhmp==0) = nan;
rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
    +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));






[rhm,fs] = fet_rhm(Trial,[],'mtchglong');
trhm = rhm.copy;
rhm.data  = log10(rhm.data);
rhm.data(rhm<-9) = nan;
rhm.data(nniz(rhm.data))=nan;
vel = xyz.vel(1,[1,2]);
vel.resample(rhm);
vnn = nniz(vel);
rhm.data = (rhm.data-repmat(nanmean(rhm(vnn,:)),[rhm.size(1),1]))...
           ./repmat(nanstd(rhm(vnn,:)),[rhm.size(1),1]);

rhmp = rhm.copy;
rhmp.data = nanmean(rhm(:,6<fs&fs<12),2);

ind = Trial.stc{'x+p'};
ind.cast('TimeSeries',rhmp);
ind.data = ind.data&nniz(rhmp);
[State, hmm, decode] = gausshmm(rhmp(ind.data),2);



headPitchState = zeros([size(rhmp,1),2]);
[~,stsInd] = max(cat(1,hmm.state.Mu));
headPitchState(ind.data==1,1) = State'==stsInd;
headPitchState(ind.data==1,2) = State'~=stsInd;


figure,hold on
plot(rhmp.data);
plot(headPitchState(:,1));
plot(headPitchState(:,2));

figure,hold on
plot(ang(headPitchState(:,1)==1,'head_back','head_front',2));
plot(headPitchState);

hps = rhmp.copy;
hps.data = headPitchState+eps;
hps.resample(xyz);

figure,hold on
eds = linspace(-pi/2,pi/2,200);
sind = hps(:,1)>0.5;
hs = bar(eds,histc(ang(sind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
sind = hps(:,2)>0.5;
hs = bar(eds,histc(ang(sind,'head_back','head_front',2),eds),'histc');
hs.FaceColor = 'r';
hs.EdgeColor = 'r';
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;


% HMM head pitch during locomotion and pause states
rhmThreshHMM = mean([prctile(ang(hps(:,1)>0.5,'head_back','head_front',2),98),prctile(ang(hps(:,2)>0.5,'head_back','head_front',2),2)]);
