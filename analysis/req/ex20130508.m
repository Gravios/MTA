
%%Primary Goal: K-nearest neighbour of xy spike rates


Trial = MTATrial('jg05-20120310',{'ufr','CluRes'},'all');

sessionLength = Trial.xyzPeriods(end)-Trial.xyzPeriods(1)+1;
uind = round(linspace(round(Trial.lfpSampleRate/Trial.xyzSampleRate),size(Trial.ufr,1),sessionLength));
Trial.ufr = ufr(uind,:);


spks = Trial.spkFeatures({'walk'},{'xyz',{'ufr',Trial.xyzSampleRate}});


nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';
pfknnvr = zeros(length(ybins),length(xbins),size(Trial.map,1));
pfknnmr = zeros(length(ybins),length(xbins),size(Trial.map,1));
pfknnvd = zeros(length(ybins),length(xbins),size(Trial.map,1));
pfknnmd = zeros(length(ybins),length(xbins),size(Trial.map,1));

nnn = zeros(length(Trial.map),1);
nss = zeros(length(Trial.map),1); 
for unit = 1:size(Trial.map,1),
    nss(unit) = size(spks(unit).ufr,1);
    nnn(unit) = round(log10(nss(unit)).^4);
    if nnn(unit)*2<length(spks(unit).ufr)&nnn(unit)>9,
        for y = 1:length(ybins),
            for x = 1:length(xbins),
                spkdist = sqrt(sum((sq(mean(spks(unit).xyz(:,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]),2))-repmat([xbins(x),ybins(y)],size(spks(unit).xyz,1),1)).^2,2));
                [~,distInd] = sort(spkdist);
                pfknnmd(y,x,unit) = median(spkdist(distInd(1:nnn(unit))));
                pfknnvd(y,x,unit) = var   (spkdist(distInd(1:nnn(unit))));
                pfknnmr(y,x,unit) = median(spks(unit).ufr(distInd(1:nnn(unit))));
                pfknnvr(y,x,unit) = var   (spks(unit).ufr(distInd(1:nnn(unit))));
            end
        end
    end
end



%% Figures

[accg,atbin] = autoccg(Trial.name);
pfw = MTAPlaceField(Trial,[],{'walk'},1);
    spkpos_head = sq(mean(spks(unit).xyz(:,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]),2));

hfig = figure(2389);
unit = 1;
while unit~=-1
    clf
    subplot2(2,2,1,1);    pfw.plot(unit,1),
    subplot2(2,2,1,2);    bar(atbin,accg(:,unit)),axis tight
                          title(['u: ' num2str(unit) ' e: ' num2str(pfw.cluMap(unit,2)) ' c: ' num2str(pfw.cluMap(unit,3))]);
    subplot2(2,2,2,1);    imagesc(xbins,ybins,pfknn(:,:,unit)*20*1000),axis xy,colorbar
                          title([num2str(unit) ' k-nn : ' num2str(nnn(unit))])
    spkpos_head = sq(mean(spks(unit).xyz(:,Trial.Model.gmi({'head_back','head_left','head_front','head_right'}),[1,2]),2));
    subplot2(2,2,2,2);    plot(spkpos_head(:,1),spkpos_head(:,2),'.'),xlim([-500,500]),ylim([-500,500])
    unit = figure_controls(hfig,unit);
end

