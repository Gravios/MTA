%% Primary Goal - establish a spatially binned correlations between
%% unit firing rate (ufr) and other behavioral (bhv) variables such
%% as speed (s), height (z), acceleration (a), head direction (ang)
%% etc ...

s = MTASession('jg05-20120310');
Trial = MTATrial(s,{'ang',{'ufr',s.xyzSampleRate},{'CluRes',s.xyzSampleRate}},'all');

bhvs = {'walk','rear'};


nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';


fet_corr = {{'ufr','vel'},...
            {'ufr','ang'},...
            {'ufr','dz'},...
            {'fufr','vel'},...
            {'fufr','ang'},...
            {'fufr','dz'}};


Trial.ufr = (Trial.ufr+min(fet.ufr{1}(:))).*20000;

for bhv = 1:length(bhvs), 
    fet.xyz{bhv} = SelectPeriods(sq(Trial.xyz(:,Trial.Model.gmi(Trial.trackingMarker),[1,2])),Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
    fet.ufr{bhv} = SelectPeriods(Trial.ufr,Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
    fet.fufr{bhv} = SelectPeriods(Filter0(gausswin(21)./sum(gausswin(21)),Trial.ufr),Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
    fet.vel{bhv} = SelectPeriods([0;Filter0(gausswin(21)./sum(gausswin(21)),Trial.vel(7,[1,2]))],Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
    fet.dz {bhv} = SelectPeriods([0;diff(Filter0(gausswin(11)./sum(gausswin(11)),sq(Trial.xyz(:,Trial.Model.gmi(Trial.trackingMarker),3))))],Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
    fet.ang{bhv} = SelectPeriods(Trial.ang(:,4,5,2),Trial.Bhv.getState(bhvs{bhv}).state,'c',1,1);
end
 

for bhv = 1:length(bhvs), 
    for unit = 1:size(Trial.map,1),
        if length(Trial.res(Trial.clu==unit))>9,
            for y = 1:length(ybins),
                for x = 1:length(xbins),
                    spkdist = sqrt(sum((fet.xyz{bhv}-repmat([xbins(x),ybins(y)],size(fet.xyz{bhv},1),1)).^2,2));
                    [dist,distInd] = sort(spkdist);                  
                    maxDistInd = find(dist<50,1,'last');
                    if length(distInd(1:maxDistInd))>10,
                    if ~isempty(find(fet.(fet_corr{f}{1}){bhv}(distInd(1:maxDistInd),unit)>0,1,'first')),keyboard,end
                        for f = 1:length(fet_corr)
                            [pfc{bhv,f}(y,x,unit),pfp{bhv,f}(y,x,unit)] = RankCorrelation(fet.(fet_corr{f}{1}){bhv}(distInd(1:maxDistInd),unit),fet.(fet_corr{f}{2}){bhv}(distInd(1:maxDistInd)));                            
                        end
                    else
                        for f = 1:length(fet_corr)
                            pfc{bhv,f}(y,x,unit) = nan;
                            pfp{bhv,f}(y,x,unit) = nan;
                        end
                    end
                end
            end
        end
    end
end


unit = 1;
f = 1;
while unit~=-1
clf,
subplot2(3,1,1,1)
imagescnan({xbins,ybins,pfc{1,f}(:,:,unit)},[],[],1),axis xy
subplot2(3,1,1,2)
pfr.plot(unit,1)
subplot2(3,1,1,3)
pfw.plot(unit,1)
title(num2str(unit))
unit = figure_controls(gcf,unit)
end



figure
for b = 1:length(bhvs),
for f = 1:length(fet_corr),
subplot2(length(bhvs),length(fet_corr),b,f)
imagescnan({ybins,xbins,pfc{b,f}(:,:,32).*(pfp{b,f}(:,:,32)<0.05)},[],[],1)
title([bhvs{b} ' ' fet_corr{f}{1} ' vs ' fet_corr{f}{2} ' unit: ' num2str(unit)])
end
end





xyoccw = xyocc(Trial,SelectPeriods(Trial.xyz(:,7,[1,2]),Trial.Bhv.getState('walk').state,'c',1,1),50);
xyoccr = xyocc(Trial,SelectPeriods(Trial.xyz(:,7,[1,2]),Trial.Bhv.getState('rear').state,'c',1,1),50);


[accg,tbin] = autoccg(Trial);



figure
unit=1;
while unit~=-1,
tpf =pfknnmd{1}(:,:,unit);
tpfr = pfknnmr{1}(:,:,unit);
tpf(tpf<=200) = 1;
tpf(tpf>200) = 0;
subplot2(3,4,2,2)
imagesc(xbins,ybins,tpf.*tpfr*20000),axis xy,colorbar

tpf =pfknnmd{2}(:,:,unit);
tpfr = pfknnmr{2}(:,:,unit);
tpf(tpf<=200) = 1;
tpf(tpf>200) = 0;
subplot2(3,4,1,2)
imagesc(xbins,ybins,tpf.*tpfr*20000),axis xy,colorbar

subplot2(3,4,2,1)
pfw.plot(unit,1)

subplot2(3,4,1,1)
pfr.plot(unit,1)

subplot2(3,4,2,3)
plot(spks{1}(unit).xyz(:,7,1),spks{1}(unit).xyz(:,7,2),'.')
xlim([-500,500])
ylim([-500,500])


subplot2(3,4,1,3)
plot(spks{2}(unit).xyz(:,7,1),spks{2}(unit).xyz(:,7,2),'.')
xlim([-500,500])
ylim([-500,500])

subplot2(3,4,2,4)
imagesc(xbins,ybins,xyoccw'),axis xy, colorbar

subplot2(3,4,1,4)
imagesc(xbins,ybins,xyoccr'),axis xy, colorbar

subplot2(3,4,3,1:4)
bar(tbin,accg(:,unit)),axis tight

title(num2str(unit))
unit = figure_controls(gcf,unit)
end


