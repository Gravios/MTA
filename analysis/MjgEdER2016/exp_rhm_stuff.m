
Trial = MTATrial('jg05-20120310');
Trial.load('stc','NN0317R');


ds = load(fullfile(Trial.spath,[Trial.filebase,'_bhvtheta_pfStats.mat']));


ds.pfmstats(2,ds.cluMap==10)

pfkbs = {};
for sts = 1:numel(ds.states),
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = ds.cluMap;
    defargs.states = ds.states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end

ufr  = Trial.ufr.copy;


label = ['IPu',num2str(unit),'p'pind];


unit = 1;
pind = 1;

figure
for s = 1:7

xinds = sq(ds.pfmstats(s,ds.cluMap==unit).patchRateInd(1,1,pind,1,:));
xinds = pfkbs{1}.adata.bins{1}(xinds(nniz(xinds)));

yinds = sq(ds.pfmstats(s,ds.cluMap==unit).patchRateInd(1,1,pind,2,:));
yinds = pfkbs{1}.adata.bins{2}(yinds(nniz(yinds)));

ppos = [xinds,yinds];
subplot(1,7,s);plot(ppos(:,1),ppos(:,2),'.'),xlim([-500,500]),ylim([-500,500])
end