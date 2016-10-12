
Trial = MTATrial('jg05-20120310');
stcMode = 'NN0317R';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
tstates = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);

Stc = Trial.load('stc',stcMode);

sum(diff([Stc{states{3},1}.data],1,2))
sum(diff([Stc{tstates{3},1}.data],1,2))

