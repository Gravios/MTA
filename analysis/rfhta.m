% rfhta.m
% saves an average of the rearing features?

slist = {'jg05-20120309','jg04-20120128','jg04-20120130','jg05-20120315'};
%slist = {'er01-20110721','jg04-20120210','jg04-20120211','jg04-20120212','jg04-20120213'};

Session = MTASession(slist{1});

rper = Session.Bhv.getState('rear').state;
[~,rfet] = rear(Session,'com');

fet = GetSegs(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),3),round(rper-2.*Session.xyzSampleRate),round(4*Session.xyzSampleRate),0);
hfet = reshape(fet,round(4*Session.xyzSampleRate),size(rper,1),2);

fet = GetSegs(rfet,round(rper-2.*Session.xyzSampleRate),round(4*Session.xyzSampleRate),0);
afet = reshape(fet,round(4*Session.xyzSampleRate),size(rper,1),2);


for si=2:length(slist),
    Session = MTASession(slist{si});

    rper = Session.Bhv.getState('rear').state;
    [~,rfet] = rear(Session,'com');

    fet = GetSegs(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),3),round(rper-2.*Session.xyzSampleRate),round(4*Session.xyzSampleRate),0);
    hfet = cat(2,hfet,reshape(fet,round(4*Session.xyzSampleRate),size(rper,1),2));

    fet = GetSegs(rfet,round(rper-2.*Session.xyzSampleRate),round(4*Session.xyzSampleRate),0);
    afet = cat(2,afet,reshape(fet,round(4*Session.xyzSampleRate),size(rper,1),2));
end

t = -2:4/size(afet,1):2;
t = [t(1:size(afet,1)/2),t((size(afet,1)/2+2):(size(afet,1)+1))];

save('/data/homes/gravio/data/analysis/RearingFetAve.mat','hfet','afet','t')



%% start here - Figures

load('/data/homes/gravio/data/analysis/RearingFetAve.mat')
%% Rear onset ave feature
% $$$ figure,hold on
% $$$ errorbar(t,sq(mean(afet(:,:,1),2)),sq(std(afet(:,:,1),[],2)))
% $$$ plot(t,sq(mean(afet(:,:,1),2)),'r')
% $$$ 
% $$$ %% Rear offset ave feature
% $$$ figure,hold on,
% $$$ errorbar(t,sq(mean(afet(:,:,2),2)),sq(std(afet(:,:,2),[],2)))
% $$$ plot(t,sq(mean(afet(:,:,2),2)),'r')


%load('/data/homes/gravio/data/analysis/RearingStats_ca3_dpf.mat')
load('/data/homes/gravio/data/analysis/RearingStats_ca1_dpf.mat')


%% each ration gives slightly different results
rratio = sq(mean(arfccg(tbin>=-2000&tbin<=-500,:,:,:,:))./mean(arfccg(tbin>=-500&tbin<=1000,:,:,:,:)));
%rratio = sq(mean(arfccg(tbin>=-1000&tbin<=0,:,:,:,:))./mean(arfccg(tbin>=0&tbin<=1000,:,:,:,:)));

%% Rear onset ratio 
figure,
plot(sq(rratio(1,1,:)),sq(rratio(1,2,:)),'.')

%% Rear offset ratio
figure,
plot(sq(rratio(2,1,:)),sq(rratio(2,2,:)),'.')




[~,sind] = sort(sq(sum(arfccg(38:48,1,1,sind),1)-sum(arfccg(53:63,1,2,sind),1)));
[~,sind] = sort(sq(sum(sq(arfccg(50:65,1,1,:)),1)));
figure,
subplot(121),imagesc(tbin(2:end-1),1:size(arfccg,4),sq(arfccg(25:75,1,1,sind)-arfccg(25:75,1,2,sind))'),colorbar
subplot(122),imagesc(tbin(2:end-1),1:size(arfccg,4),unity(sq(arfccg(25:75,1,1,sind)))'),colorbar
