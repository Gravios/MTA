configure_default_args();
EgoProCode2D_load_data(); 

pclr = cool(3);

binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;

hbaBinEdg = [-1.2,-0.2,0.2,1.2];
hbaBinCtr = mean([hbaBinEdg(1:end-1);hbaBinEdg(2:end)]);

hbaBin.edges = [-1.2,-0.2,0.2,1.2];
hbaBin.centers = mean([hbaBin.edges(1:end-1);hbaBin.edges(2:end)]);
hbaBin.count = numel(hbaBin.centers);

phzBin.edges = linspace(0.5,2*pi-0.5,4);
phzBin.centers = (binPhzs(1:end-1)+binPhzs(2:end))./2;
phzBin.count = numel(phzBin.centers);

havBin.edges = [-0.3,-0.15,0.15,0.3];

hvang = filter(copy(xyz),'ButFilter',4,2,'low');
xycoor = cat(2,...
             hvang(:,'spine_upper',[1,2])-hvang(:,'bcom',[1,2]),... 
             hvang(:,'nose',[1,2])-hvang(:,'hcom',[1,2]));
hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% Positive: CCW (Left)     Negative: CW (Right) 
dc.hvang.data = circ_dist(circshift(hvang.data(:,2),-10),...
                          circshift(hvang.data(:,2),+10));


% $$$ figure();
ucounter = 1;
egoMeanRmapPosHba = [];
egoMaxRmapPosHba = [];
egoSizeHba = [];
egoMeanRmapRateHba = [];
egoMaxRmapRateHba = [];

for t = 1:numel(Trials)
for u = 1:numel(unitsEgo{t})
unit = unitsEgo{t}(u);
for p = 1:3;
for a = 1:3,
binSubsetX = abs(pfs{t}{p,a}.adata.bins{1})<300;
binSubsetY = abs(pfs{t}{p,a}.adata.bins{2})<300;
mapPosition = cell([1,2]);
[mapPosition{:}] = ndgrid(pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
                          pfs{t}{p,a}.adata.bins{2}(binSubsetY));
mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
rmap = plot(pfs{t}{p,a},unit,[],[],[],false);
rmap = rmap(binSubsetX,binSubsetY);
nanmap = double(~isnan(rmap));
nanmap(nanmap==0) = nan;
rmap = rmap.*fliplr(nanmap);
rmap(rmap<2) = 0;
% $$$ subplot2(3,3,4-p,a);
% $$$ hold('on');
% $$$ imagescnan({pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
% $$$                    pfs{t}{p,a}.adata.bins{2}(binSubsetY),...
% $$$                    rmap'});
nrmap =rmap./sum(rmap(:),'omitnan');
rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
egoSizeHba(ucounter,p,a) = sum(nniz(nrmap(:)));
egoMeanRmapRateHba(ucounter,p,a,:) = mean(rmap(nniz(rmap(:))),'omitnan');

% $$$ plot(rmapCenter(1),rmapCenter(2),'*m')
% $$$ axis('tight')
% $$$ Lines([],0,'g');
% $$$ Lines(0,[],'g');
egoMeanRmapPosHba(ucounter,p,a,:) = rmapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);

[~,maxPos] = max(nrmap(:));
if ~isempty(maxPos)
[maxX,maxY] = ind2sub(size(nrmap),maxPos);
egoMaxRmapPosHba(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
else
egoMaxRmapPosHba(ucounter,p,a,:) = nan([1,1,1,2]);
end
end
end
ucounter = ucounter+1;
end
end




% $$$ figure();






% $$$ 
% $$$ ucounter = 1;
% $$$ egoPfsStatsHba.control.meanPos = []; % egoMeanRmapPosHba
% $$$ egoPfsStatsHba.control.peakPos = []; % egoMaxRmapPosHba
% $$$ egoPfsStatsHba.control.size = []; %egoSizeHba = []; 
% $$$ egoPfsStatsHba.control.maxRate = [] %egoMeanRmapRateHba = [];
% $$$ egoPfsStatsHba.control.meanRate = []; %egoMaxRmapRateHba = [];
% $$$ 
% $$$ for t = 1:numel(Trials)
% $$$     for u = 1:numel(unitsEgo{t})
% $$$         unit = unitsEgo{t}(u);
% $$$         for p = 1:3;
% $$$             for a = 1:3,
% $$$                 binSubsetX = abs(pfs{t}{p,a}.adata.bins{1})<300;
% $$$                 binSubsetY = abs(pfs{t}{p,a}.adata.bins{2})<300;
% $$$                 mapPosition = cell([1,2]);
% $$$                 [mapPosition{:}] = ndgrid(pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
% $$$                                           pfs{t}{p,a}.adata.bins{2}(binSubsetY));
% $$$                 mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
% $$$                 rmap = plot(pfs{t}{p,a},unit,[],[],[],false);
% $$$                 rmap = rmap(binSubsetX,binSubsetY);
% $$$                 nanmap = double(~isnan(rmap));
% $$$                 nanmap(nanmap==0) = nan;
% $$$                 rmap = rmap.*fliplr(nanmap);
% $$$                 rmap(rmap<2) = 0;
% $$$                 nrmap =rmap./sum(rmap(:),'omitnan');
% $$$                 rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
% $$$                 egoPfsStatsHba.control.size(ucounter,p,a) = sum(nniz(nrmap(:)));
% $$$                 egoPfsStatsHba.control.maxRate(ucounter,p,a,:) = mean(rmap(nniz(rmap(:))),'omitnan');
% $$$                 egoPfsStatsHba.control.meanPos(ucounter,p,a,:) = rmapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);
% $$$                 [~,maxPos] = max(nrmap(:));
% $$$                 if ~isempty(maxPos)
% $$$                     [maxX,maxY] = ind2sub(size(nrmap),maxPos);
% $$$                     egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
% $$$                 else
% $$$                     egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = nan([1,1,1,2]);
% $$$                 end
% $$$             end
% $$$         end
% $$$         ucounter = ucounter+1;
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ ucounter = 1;
% $$$ egoPfsStatsHba.shuffle.meanPos = []; % egoMeanRmapPosHba
% $$$ egoPfsStatsHba.shuffle.peakPos = []; % egoMaxRmapPosHba
% $$$ egoPfsStatsHba.shuffle.size = []; %egoSizeHba = []; 
% $$$ egoPfsStatsHba.shuffle.maxRate = [] %egoMeanRmapRateHba = [];
% $$$ egoPfsStatsHba.shuffle.meanRate = []; %egoMaxRmapRateHba = [];
% $$$ 
% $$$ for t = 1:numel(Trials)
% $$$     for u = 1:numel(unitsEgo{t})
% $$$         unit = unitsEgo{t}(u);
% $$$         for p = 1:3;
% $$$             for a = 1:3,
% $$$                 binSubsetX = abs(pfs{t}{p,a}.adata.bins{1})<300;
% $$$                 binSubsetY = abs(pfs{t}{p,a}.adata.bins{2})<300;
% $$$                 mapPosition = cell([1,2]);
% $$$                 [mapPosition{:}] = ndgrid(pfs{t}{p,a}.adata.bins{1}(binSubsetX),...
% $$$                                           pfs{t}{p,a}.adata.bins{2}(binSubsetY));
% $$$                 mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
% $$$                 rmap = plot(pfs{t}{p,a},unit,[],[],[],false);
% $$$                 rmap = rmap(binSubsetX,binSubsetY);
% $$$                 nanmap = double(~isnan(rmap));
% $$$                 nanmap(nanmap==0) = nan;
% $$$                 rmap = rmap.*fliplr(nanmap);
% $$$                 rmap(rmap<2) = 0;
% $$$                 nrmap =rmap./sum(rmap(:),'omitnan');
% $$$                 rmapCenter = sq(sum(sum(bsxfun(@times,nrmap,mapPosition),'omitnan'),'omitnan'))';
% $$$                 egoPfsStatsHba.control.size(ucounter,p,a) = sum(nniz(nrmap(:)));
% $$$                 egoPfsStatsHba.control.maxRate(ucounter,p,a,:) = mean(rmap(nniz(rmap(:))),'omitnan');
% $$$                 egoPfsStatsHba.control.meanPos(ucounter,p,a,:) = rmapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);
% $$$                 [~,maxPos] = max(nrmap(:));
% $$$                 if ~isempty(maxPos)
% $$$                     [maxX,maxY] = ind2sub(size(nrmap),maxPos);
% $$$                     egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
% $$$                 else
% $$$                     egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = nan([1,1,1,2]);
% $$$                 end
% $$$             end
% $$$         end
% $$$         ucounter = ucounter+1;
% $$$     end
% $$$ end
% $$$ 
% $$$ 

%%% new method



ucounter = 1;
egoPfsStatsHba.control.meanPos = []; % egoMeanRmapPosHba
egoPfsStatsHba.control.peakPos = []; % egoMaxRmapPosHba
egoPfsStatsHba.control.size = []; %egoSizeHba = []; 
egoPfsStatsHba.control.maxRate = [] %egoMeanRmapRateHba = [];
%egoPfsStatsHba.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:3;
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = rmap{trialInd}( abs(xbins)<=300,...
                                          abs(ybins)<=300,...
                                          unit,...
                                          phzInd,...
                                          hbaInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoPfsStatsHba.control.size   (ucounter,phzInd,hbaInd)   = sum(nniz(nratemap(:)));
                egoPfsStatsHba.control.maxRate(ucounter,phzInd,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoPfsStatsHba.control.meanPos(ucounter,phzInd,hbaInd,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = mapPosition(maxX,maxY,:);
                else
                    egoPfsStatsHba.control.peakPos(ucounter,p,a,:) = nan([1,1,1,2]);
                end
            end
        end
        ucounter = ucounter+1;
    end
end

ucounter = 1;
egoPfsStatsHba.shuffle.meanPos = []; % egoMeanRmapPosHba
egoPfsStatsHba.shuffle.peakPos = []; % egoMaxRmapPosHba
egoPfsStatsHba.shuffle.size = []; %egoSizeHba = []; 
egoPfsStatsHba.shuffle.maxRate = [] %egoMeanRmapRateHba = [];
%egoPfsStatsHba.shuffle.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:size(rmapShuff{trialInd},4)
            for hbaInd = 1:size(rmapShuff{trialInd},5)
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(rmapShuff{trialInd},6)
                    ratemap = sq(rmapShuff{trialInd}( abs(xbins)<=300,...
                                                      abs(ybins)<=300,...
                                                      unit,           ...
                                                      phzInd,         ...
                                                      hbaInd,         ...
                                                      iter));
                    nanmap = double(~isnan(ratemap));
                    nanmap(nanmap==0) = nan;
                    ratemap = ratemap.*fliplr(nanmap);
                    ratemap(ratemap<2) = 0;
                    nratemap =ratemap./sum(ratemap(:),'omitnan');
                    ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                    egoPfsStatsHba.shuffle.size   (ucounter,phzInd,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoPfsStatsHba.shuffle.maxRate(ucounter,phzInd,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoPfsStatsHba.shuffle.meanPos(ucounter,phzInd,hbaInd,iter,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(t,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoPfsStatsHba.shuffle.peakPos(ucounter,p,a,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoPfsStatsHba.shuffle.peakPos(ucounter,p,a,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        end% phzInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


u = 76
figure,histogram(sq(egoPfsStatsHba.shuffle.meanPos(u,1,3,:,2))- ...
                 sq(egoPfsStatsHba.shuffle.meanPos(u,1,2,:,2)),linspace([-10,10,30]));

figure,histogram((egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,3,2)-mean(egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,1,3,:,2),4))...
    ./std(egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,1,3,:,2),[], ...
          4),linspace([-30,60,60]))


egoPfsStatsHba.perm.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2])
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoPfsStatsHba.perm.ca1.zscore(:,phzInd,hbai,x) = ( (egoPfsStatsHba.control.meanPos(unitsEgoCA1,phzInd,hbaInd(1),x) - egoPfsStatsHba.control.meanPos(unitsEgoCA1,phzInd,hbaInd(2),x))-(mean(egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),4)))./std(egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoPfsStatsHba.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end


figure();
hold('on');
for phzInd = 1:phzBin.count
    subplot2(phzBin.count, 1, phzBin.count+1-phzInd, 1);
    hold('on')
    for hbaInd = 1:hbaBin.count
        cdfplot(egoPfsStatsHba.perm.ca1.zscore(:,phzInd,hbaInd,2));
    end
    Lines(3.35279,[],'r');
    Lines(-3.35279,[],'r');
    xlim([-20,20]);
end

% $$$ Lines(4.314457,[],'r');
% $$$ Lines(4.798339,[],'r');

% 1-(1-0.05)^(1/numel(unitsEgoCA1)) -> 0.000410262174549647 -> zscore > 3.352795
% 1-(1-0.001)^(1/numel(unitsEgoCA1)) -> 8.00397063671632e-06 -> zscore > 4.314457
% 1-(1-0.0001)^(1/numel(unitsEgoCA1)) -> 8.00397063671632e-06 -> zscore > 4.798339

figure();
hold('on');
plot(egoPfsStatsHba.perm.zscoreL, egoPfsStatsHba.perm.zscoreR,'.');
circle(0,0,3.35279,'r');
circle(0,0,4.314457,'r');
circle(0,0,4.798339,'r');
grid('on');


figure();
hold('on');
plot(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,2,2)-egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,1,2),...
     egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,3,2)-egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,2,2),...
     '.');
grid('on');


figure();
hold('on');
plot(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,2,2), ...
     egoPfsStatsHba.perm.zscoreC,'.');
plot(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,3,2), ...
     egoPfsStatsHba.perm.zscoreR,'.r');
plot(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,1,2), ...
     egoPfsStatsHba.perm.zscoreL,'.g');
%[RHO,PVAL]=corr(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,2,2),egoPfsStatsHba.perm.zscoreC)
[RHO,PVAL]=corr(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,1,2), ...
     egoPfsStatsHba.perm.zscoreL)
[RHO,PVAL]=corr(egoPfsStatsHba.control.meanPos(unitsEgoCA1,1,3,2), ...
     egoPfsStatsHba.perm.zscoreR)
xlim([-30,30]);