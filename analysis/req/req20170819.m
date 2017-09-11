% MjgER2016 exploration 

% Label the sessions used in the analysis

label_bhv_msnn('BHV_S4H5');
label_ripples(Trial);


Trial = MTATrial.validate('jg05-20120310.cof.all');


stc = Trial.load('stc','msnn_ppsvd');
xyz = Trial.load('xyz');



% COMPUTE ccg betwee behavior transiton and theta trasitions
tper = stc{'theta',xyz.sampleRate};

wper = get_state_transitions(stc,Trial,{'pause','walk'},[],xyz);
wper = get_state_transitions(stc,Trial,{'walk','pause'},[],xyz);

dind = wper(2:end,1)-wper(1:end-1,2);
wper([1000;dind]<240,:) = [];

[tccg,t] = CCG([tper(:,1);wper(:,1)],...
               [ones([size(tper,1),1]);2*ones([size(wper,1),1])],...
               20,30,xyz.sampleRate,[1,2],'count');
figure,bar(t,tccg(:,2,1))


% COMPUTE ccg betwee behavior transiton and spw and 
sper = stc{'spw',xyz.sampleRate};
wper = get_state_transitions(stc,Trial,{'pause','sit'},[],xyz);

dind = wper(2:end,1)-wper(1:end-1,2);
wper([1000;dind]<240,:) = [];
[sccg,t] = CCG([sper(:,1);wper(:,1)],...
               [ones([size(sper,1),1]);2*ones([size(wper,1),1])],...
               20,30,xyz.sampleRate,[1,2],'count');
figure,bar(t,sccg(:,2,1))



% COMPUTE spw rate for each behavior
srate = [];
states = {'walk','rear','turn','pause','groom','sit','shake'};
for s = 1:numel(states)
    bper = [stc{states{s}}];
    bdur = sum((diff(bper.data,1,2)+1)./bper.sampleRate);
    ss = [stc{[states{s},'&spw'],xyz.sampleRate}];
    srate(s) = size(ss,1)/bdur;
end


% RELATE spw to slow respiration
ang = create(MTADang,Trial,xyz);
lang = MTADfet.encapsulate(Trial,...
                           ButFilter(unity(ang(:,2,4,3)),3,[0.7,1.5]./(ang.sampleRate.*0.5),'bandpass'),...
                           ang.sampleRate,...
                           'Slow Respiration',...
                           'sresp',...
                           'x');

%figure,plot(lang.data), plotSTC(stc);

plang = lang.phase([0.7,3]);
plm = LocalMinima(plang.data,30,-pi+0.2);

%figure,plot(lang.data), plotSTC(stc); Lines(plm,[],'r');

sper = stc{'spw',xyz.sampleRate};
sind = sper(2:end,1)-sper(1:end-1,2);
sper.data([1000;sind]<360,:) = [];


plm(~WithinRanges(plm,stc{'s'}(:,:)))=[];
plms = plm(lang(plm)>-0.008&lang(plm)<-0.001);
plms = plms(1:3:end);


%figure,rose(plang(sper(1:1:end,1)))

%figure,plot(lang.data), plotSTC(stc); Lines(plms,[],'r');

[pccg,t] = CCG([sper(:,1);plms],...
               [ones([size(sper,1),1]);2*ones([size(plms,1),1])],...
               10,60,xyz.sampleRate,[1,2],'count');
figure,
subplot(311);
bar(t,pccg(:,1,1));
subplot(312);
bar(t,pccg(:,2,1));
subplot(313);
bar(t,pccg(:,2,2));


% EXAMINE coherence between lfp and rhm
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',72);
rhm = fet_rhm(Trial);
lfp.resample(rhm);

prhm = rhm.phase([6,14]);
plfp = lfp.phase([6,14]);

plfpm = LocalMinima(plfp.data,5,-pi+0.2);
plfpm(~WithinRanges(plfpm,[stc{'w'}(:,:)]))=[];



frhm = prhm.copy;
flfp = plfp.copy;

frhm.data = circ_dist(circshift(prhm.data,-1),prhm.data)./(2.*pi).*prhm.sampleRate;
flfp.data = circ_dist(circshift(plfp.data,-1),plfp.data)./(2.*pi).*prhm.sampleRate;

ind = stc{'p'};
figure,plot(flfp(ind),frhm(ind),'.');



