sesList = get_session_list('MjgEdER2016_bhv');
stcMode = 'NN0317R';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);


ds = {};
for s = 1:7%numel(sesList),
    Trial = MTATrial.validate(sesList(s));
    ds{s} = load(fullfile(Trial.spath,[Trial.filebase,'_bhvtheta_pfStats.mat']));
end

Trials = arrayfun(@MTATrial.validate,sesList,'UniformOutput',false);

pfstats = cellfun(@(Trial) load(fullfile(Trial.spath,[Trial.filebase,'_bhvtheta_pfStats.mat'])),...
                  Trials,'UniformOutput',false);






% Characterize mean and maximum firing rates between states
ind = ':';
figure, hold on,
for s = [1,2,3,4,5]
    plot([ds{s}.pfkstats(1,ind).peakFR],[ds{s}.pfkstats(2,ind).peakFR],'.b')
end

ind = ismember(ds{6}.cluMap,ds{7}.cluMap);

figure,plot([ds{6}.pfkstats(4,ind).peakFR],[ds{7}.pfkstats(4,ind).peakFR],'.r')
hold on,plot([ds{6}.pfmstats(4,ind).peakFR],[ds{7}.pfmstats(4,ind).peakFR],'.g')

figure,plot([ds{6}.pfkstats(4,ind).patchMFR(:,],[ds{7}.pfkstats(4,ind).patchMFR],'.r')

figure,plot([ds{6}.pfkstats(,ind).peakFR],[ds{7}.pfkstats(4,ind).peakFR],'.r')


% Characterize the influence of rhythmic head motion (rhm) on unit firing rate

npfkbs = {};
for sts = 1:numel(ds.states),
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = ds.cluMap;
    defargs.states = ds.states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end



unit = ds.cluMap(9);
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

figure
for s = 1:7
    subplot(1,7,s);
    pfkbs{s}.plot(unit);
end




unit = cluMap(8);
pind = 1;

figure
for s = 1:7

    xinds = sq(pfmstats(s,cluMap==unit).patchRateInd(1,1,pind,1,:));
    xinds = pfkbs{1}.adata.bins{1}(xinds(nniz(xinds)));

    yinds = sq(pfmstats(s,cluMap==unit).patchRateInd(1,1,pind,2,:));
    yinds = pfkbs{1}.adata.bins{2}(yinds(nniz(yinds)));

    ppos = [xinds,yinds];
    subplot(1,7,s);plot(ppos(:,1),ppos(:,2),'.'),xlim([-500,500]),ylim([-500,500])
end





figure
for s = 1:7

    xinds = sq(pfkstats(s,cluMap==unit).patchRateInd(1,1,pind,1,:));
    xinds = pfkbs{1}.adata.bins{1}(xinds(nniz(xinds)));

    yinds = sq(pfkstats(s,cluMap==unit).patchRateInd(1,1,pind,2,:));
    yinds = pfkbs{1}.adata.bins{2}(yinds(nniz(yinds)));

    ppos = [xinds,yinds];
    subplot(1,7,s);plot(ppos(:,1),ppos(:,2),'.'),xlim([-500,500]),ylim([-500,500])
end


figure
for s = 1:7
    subplot(1,7,s);
    pfkbs{s}.plot(unit);
end



ufr  = Trial.ufr.copy;
label = ['IPu',num2str(unit),'p'pind];


ds.pfmstats.PatchRateInd