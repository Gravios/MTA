% MjgER2016 exploration 

% Label the sessions used in the analysis
label_bhv_msnn('BHV_S4H5');



Trial = MTATrial.validate('jg05-20120310.cof.all');
stc = Trial.load('stc','msnn');
xyz = Trial.load('xyz');

tper = Trial.stc{'theta',xyz.sampleRate};
wper = Trial.stc{'walk',xyz.sampleRate};
wper = get_state_transitions(Trial.stc,Trial,{'pause','walk'},[],xyz);


[tccg,t] = CCG([tper(:,1);wper(:,1)],...
               [ones([size(tper,1),1]);2*ones([size(wper,1),1])],...
               20,20,xyz.sampleRate,[1,2],'count');


figure,bar(t,tccg(:,1,2))

lfp = Trial.load('lfp',72);
rhm = fet_rhm(Trial);
lfp.resample(rhm);

prhm = rhm.phase([6,14]);
plfp = lfp.phase([6,14]);

frhm = prhm.copy;
flfp = plfp.copy;

frhm.data = circ_dist(circshift(prhm.data,-1),prhm.data)./(2.*pi).*prhm.sampleRate;
flfp.data = circ_dist(circshift(plfp.data,-1),plfp.data)./(2.*pi).*prhm.sampleRate;

ind = stc{'p'};
figure,plot(flfp(ind),frhm(ind),'.');



