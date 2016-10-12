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


for bhv = 1:length(bhvs),
% $$$ 
% $$$     spks{bhv} = Trial.spkFeatures({bhvs{bhv}},{'ufr','ang','xyz' ,'vel','acc'});
% $$$ 
% $$$     pfknnvd{bhv} = zeros(length(ybins),length(xbins),size(Trial.map,1));
% $$$     pfknnmd{bhv} = zeros(length(ybins),length(xbins),size(Trial.map,1));
% $$$ 
% $$$     pfknncr{bhv} = zeros(length(ybins),length(xbins),size(Trial.map,1));
 
    pfknnvr{bhv} = zeros(length(ybins),length(xbins),size(Trial.map,1));
    pfknnmr{bhv} = zeros(length(ybins),length(xbins),size(Trial.map,1));
% $$$ 
% $$$     nnn{bhv} = zeros(length(Trial.map),1);
% $$$     nss{bhv} = zeros(length(Trial.map),1); 

    for unit = 1:size(Trial.map,1),
% $$$         nss{bhv}(unit) = size(spks{bhv}(unit).ufr,1);
% $$$         nnn{bhv}(unit) = round(log10(nss{bhv}(unit)).^4);
        if nnn{bhv}(unit)*2<length(spks{bhv}(unit).ufr)&nnn{bhv}(unit)>9,
            for y = 1:length(ybins),
                for x = 1:length(xbins),
                    spkdist = sqrt(sum((sq(mean(spks{bhv}(unit).xyz(:,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]),2))-repmat([xbins(x),ybins(y)],size(spks{bhv}(unit).xyz,1),1)).^2,2));
                    [~,distInd] = sort(spkdist);

% $$$                     pfknnmd{bhv}(y,x,unit) = median(spkdist(distInd(1:nnn{bhv}(unit))));
% $$$                     pfknnvd{bhv}(y,x,unit) = var   (spkdist(distInd(1:nnn{bhv}(unit))));
% $$$ 
% $$$                     kncr{bhv}(y,x,unit) = corr(spks{bhv}(unit).ufr(distInd(1:nnn{bhv}(unit))),spks{bhv}(unit).xyz(distInd(1:nnn{bhv}(unit)),7,3));
                    %pfknnvr{bhv}(y,x,unit) =   (spks{bhv}(unit).ufr(distInd(1:nnn{bhv}(unit))));

                    pfknnmr{bhv}(y,x,unit) = median(spks{bhv}(unit).ufr(distInd(1:nnn{bhv}(unit))));
                    pfknnvr{bhv}(y,x,unit) = var   (spks{bhv}(unit).ufr(distInd(1:nnn{bhv}(unit))));

                end
            end
        end
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


