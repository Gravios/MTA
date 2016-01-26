

Trial = MTATrial('jg05-20120317');

fet = fet_tsne_rev5(Trial,20);

ifet = fet(nniz(fet),1);

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


