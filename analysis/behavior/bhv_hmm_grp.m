


trialList = get_session_list('jg05');

fet = [];
tfs = [];
for t = trialList,
    Trial = MTATrial(t.name,t.trialName,t.mazeName);
    ifet = fet_lgr_rev3(Trial,30);
    tfs = cat(2,tfs,ifet.size(1));
    fet = cat(1,fet,ifet.data);
end

hmm_states = 64;
nind = nniz(fet);
log_inds = [4,9:18];
fet(nind,log_inds) = log10(fet(nind,log_inds));

[State, hmm, decode] = gausshmm(fet(nind,:),hmm_states);
%[~, rearstate ] = max(cat(1,hmm.state.Mu));
 
figure,plot(decode.q_star)

det = zeros([size(fet,1),1]);
det(nind) = decode.q_star;


mdet = MTADxyz('data',det(end-tfs(4):end),'sampleRate',ifet.sampleRate);

figure,plot(mdet.data)
Lines(Trial.stc{'r',mdet.sampleRate}(:),[],'r');
Lines(Trial.stc{'w',mdet.sampleRate}(:),[],'b');
Lines(Trial.stc{'n',mdet.sampleRate}(:),[],'g');
Lines(Trial.stc{'m',mdet.sampleRate}(:),[],'m');


edgs = linspace(.5,36.5,37);
sts = 'rwnms';
figure
for s = 1:numel(sts),
    subplot(numel(sts),1,s);
    bar(edgs,histc(mdet(Trial.stc{sts(s)}),edgs),'histc');
    title(Trial.stc{sts(s)}.label)
    ylim([0,1e4])
end



figure,hist(fet(:,9),100)

figure,hist(log10(fet(:,11)),100)

edgs = linspace(-3,1,100);
ind = Trial.stc{'a-r'};
figure,hold on,
bar(edgs,histc(log10(fet(nind,1)),edgs),'histc')
ind = Trial.stc{'r'};
bar(edgs,histc(log10(fet(ind,11)),edgs),'histc')


edgs = linspace(-3,1,100);
i = 9;
mfet = MTADxyz('data',log10(fet(:,i)),'sampleRate',fet.sampleRate);
figure,sbar(Trial,mfet,'w',edgs)






% newlogs 4,9,10