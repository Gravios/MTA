
Trial = MTATrial('jg05-20120310');
stcMode = 'NN0317R';
pfStatsTag = '';

states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
tstates = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
Stc = Trial.load('stc',stcMode);


sum(diff([Stc{states{3},1}.data],1,2))

sum(diff([Stc{states{1},1}.data],1,2))
sum(diff([Stc{tstates{4},1}.data],1,2))+sum(diff([Stc{tstates{5},1}.data],1,2))
figure,plotSTC(Stc);


analysisFileName = fullfile(Trial.spath,[Trial.filebase,pfStatsTag,'_pfStats.mat']);




for sts = 1:numel(states),
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end
