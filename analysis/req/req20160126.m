

Trial = MTATrial('jg05-20120317');
fxyz = Trial.load('xyz');
fxyz.resample(30);
fxyz.filter('ButFilter',7,10);
fvxy = fxyz.vel([],[1,2]);
fang = create(MTADang,Trial,fxyz);

%figure,plot(fvxy(:,1));

%fet = fet_tsne_rev5(Trial,20);


ind = Trial.stc{'a'}.cast('TimeSeries');
ind.resample(30);

ifet = fvxy(~~ind.data,1);

[tState, thmm, decode] = gausshmm(ifet,3);

state = zeros(size(ind));
state(~~ind.data) = tState;
figure;sp=[];
sp(end+1)=subplot(211);plot(state);
Lines(Trial.stc{'w',30}(:),[],'c');
sp(end+1)=subplot(212);plot(fvxy(:,1));
linkaxes(sp,'x');

ifet = fang(~~ind.data,3,4,2);

[atState, athmm, adecode] = gausshmm(ifet,3);

astate = zeros(size(ind));
astate(~~ind.data) = atState;
figure;sp=[];
sp(end+1)=subplot(211);plot(astate);
Lines(Trial.stc{'r',30}(:),[],'c');
sp(end+1)=subplot(212);plot(fang(:,3,4,2));
linkaxes(sp,'x');



State = {};

papo = parpool(3);

%for fet

tState = cell([1,3]);
thmm   = cell([1,3]);
decode = cell([1,3]);
for itr = 1:3,    
    [tState{itr}, thmm{itr}, decode{itr}] = gausshmm(ifet,itr+1);
end


ifet = fet(nniz(fet),2);
[State, hmm, decode] = gausshmm(ifet,4);

ifet = fet(nniz(fet),3);
[State, hmm, decode] = gausshmm(ifet,4);

ifet = fet(nniz(fet),4);
[State, hmm, decode] = gausshmm(ifet,4);


% Peak Picking
Trial = MTATrial('jg05-20120317');
fxyz = Trial.load('xyz');
fxyz.filter('ButFilter',3,20);
fvxy = fxyz.vel([],[1,2]);
fang = create(MTADang,Trial,fxyz);

fz = Trial.load('xyz');
fz.filter('ButFilter',3,[1,10],'bandpass');
[e,in] = findpeaks(abs(diff(fz(:,1,3))),'MinPeakHeight',0.5);
figure,plot(abs(diff(fz(:,1,3)))),hold on,plot(in,e,'.');

ipd =[];
wp = in(e>20);
for p = 1:numel(wp)-1;
    %ipd(p) = sqrt(sum([fxyz(wp(p+1),1,[1,2])-fxyz(wp(p),1,[1,2])].^2,3));
    ipd(p) = sqrt(sum([fxyz(wp(p+1),1,3)-fxyz(wp(p),1,3)].^2,3));
end
figure,hist(log10(ipd),100)



