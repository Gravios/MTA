


%% figure Permutation Test part 1

Trial = MTATrial('jg05-20120310');
pfs = MTAApfsPerm(Trial,'numIter',1000);
units = pfs.data.clu;
pfr = MTAApfs(Trial,units,'rear');
pfw = MTAApfs(Trial,units,'walk');

hfig = figure(471);
set(hfig,'position',[2,402,1598,369]);
unit = units(1);
while unit~=-1,
unit = 65;
    mr = sort([pfr.data.rateMap(:,unit==units,1);pfw.data.rateMap(:,unit==units,1)],'descend');
    mr = mr(nniz(mr));
    mr = median(mr(1:20));
    %subplot(131),pfr.plot(unit,[],1,[0,mr]);
    %subplot(132),pfw.plot(unit,[],1,[0,mr]);
    subplot(133),pfs.plot(unit,'sigks',1);
    title(['unit: ' num2str(unit)])
    unit = figure_controls(hfig,unit,units);
end


%% figure Permutation Test part 1

Trial = MTATrial('jg05-20120310');

pfs = MTAApfsPerm(Trial,[],[],true,'numIter',1000);
dbstop in MTAApfsPerm at 228


rper = Trial.stc{'r'};
wper = Trial.stc{'w'};

rper.cast('TimeSeries');
wper.cast('TimeSeries');


durStateA = sum(rper.data);
durStateB = sum(rper.data);


pooledState = rper.data+(wper.data*2);
psi = find(pooledState);

% STATE LABEL CONTROL
w = 6000;
hfig = figure(3857),hold on
c = 'rb';
for i = psi(w:w+3000)'
p = patch([-.5,-.5,.5,.5]+i,[0,1,1,0],c(pooledState(i)));
set(p,'EdgeColor',c(pooledState(i)));
end
set(gca,'xticklabelmode','manual');
set(gca,'yticklabelmode','manual');
set(gca,'xticklabel',{});
set(gca,'yticklabel',{});
set(hfig,'position',[2,670,1598,101]);


% STATE LABEL RANDOMIZED
rind = randperm(length(psi));
permStateA = psi(rind(1:durStateA));
permStateB = psi(rind(durStateA+1:end));

pooledState = zeros(size(pooledState));
pooledState(permStateA) = 1;
pooledState(permStateB) = 2;
psir = find(pooledState);


w = 6000;
hfig = figure(3857),hold on
c = 'rb';
for i = psir(w:w+3000)'
p = patch([-.5,-.5,.5,.5]+i,[0,1,1,0],c(pooledState(i)));
set(p,'EdgeColor',c(pooledState(i)));
end
set(gca,'xticklabelmode','manual');
set(gca,'yticklabelmode','manual');
set(gca,'xticklabel',{});
set(gca,'yticklabel',{});
set(hfig,'position',[2,670,1598,101]);


dbstop in MTAApfsPerm at 228
pfs = MTAApfsPerm(Trial,[],[],true,'numIter',1000);
hfig = figure(23487);
set(hfig,'position',[2,402,1598,369]);
subplot(131),imagescnan({Pfs.adata.bins{1},Pfs.adata.bins{2}, ...
                     reshape(rateMapA,Pfs.adata.binSizes')'},[],[],true,[0,0,0]),axis xy
subplot(132),imagescnan({Pfs.adata.bins{1},Pfs.adata.bins{2}, ...
                     reshape(rateMapB,Pfs.adata.binSizes')'},[],[],true,[0,0,0]),axis xy 
subplot(133),imagescnan({Pfs.adata.bins{1},Pfs.adata.bins{2}, ...
                     reshape(rateMapA-rateMapB,Pfs.adata.binSizes')'},[],[],true,[0,0,0]),axis xy 

% STATE shuff distrib
[~,mind] = min(pfs.data.rateMap(:,pfs.data.clu==65,1));
figure,hist(sq(pfs.data.rateMap(mind,pfs.data.clu==65,:)),50)
% $$$ set(gca,'yticklabelmode','manual');
% $$$ set(gca,'yticklabel',{});