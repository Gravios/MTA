



t = 20;
Trial = Trials{t};
unitSubset = units{t};

pft = pfs_2d_theta(Trial,unitSubset);

bstates = {'theta','rear','hloc+hpause&theta','lloc+lpause&theta'};
bstatesLabels = {'theta','rear','highBhv','lowBhv'};
pfs = {};
for s = 1:numel(bstates),
    pfs{s} = MTAApfs(Trial,'tag',['HRZxTHP_',bstatesLabels{s}]);
end

dfs = req20180123_ver5(Trial,[],'13');

phzOrder = [9:16,1:8];

figure();
for unit = unitSubset,
    clf();
    subplot(1,6,1);
    plot(pft,unit,'mean','text');
    title(['pft unit:',num2str(unit)]);
    subplot(1,6,2);
    plot(dfs{1},unit,'mean','text',[],false);
    title(['pfb unit:',num2str(unit)]);    
    rmap = {};
    for s = 1:numel(bstates),        
        rmap{s} = plot(pfs{s},unit,1,'text',[],false);    
    end
    rmax = max(cell2mat(cf(@(m) max(reshape(m,[],1)), rmap)));
    for s = 1:numel(bstates),    
        subplot(1,6,s+2);
        imagescnan({pfs{s}.adata.bins{1},linspace(-2*pi,2*pi,32),repmat(rmap{s}(:,phzOrder),1,2)'},...
                   [0,rmax]);
        axis('xy');
        title(bstatesLabels{s});
    end
    waitforbuttonpress();
end