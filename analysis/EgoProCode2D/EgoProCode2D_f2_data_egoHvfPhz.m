
%egoHvfPhzRmaps = load(fullfile(Trials{1}.path.project,'analysis','EgoProCode2D_compute_egoratemaps_conditioned_on_hvf_and_phz_DATA.mat'));
egoHvfPhzRmaps = load(fullfile(Trials{1}.path.project,'analysis','EgoProCode2D_compute_egoratemaps_conditioned_on_hvf_and_phz_DATA.mat'));


xbins = egoHvfPhzRmaps.xbins;
ybins = egoHvfPhzRmaps.ybins;

ucounter = 1;
egoHvfPhz.control.meanPos = [];   %egoMeanRmapPosHvf
egoHvfPhz.control.peakPos = [];   %egoMaxRmapPosHvf
egoHvfPhz.control.size = [];      %egoSizeHvf = []; 
egoHvfPhz.control.maxRate = [];    %egoMeanRmapRateHvf = [];
%egoHvfPhz.control.meanRate = []; %egoMaxRmapRateHvf = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHvfPhzRmaps.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgoHvf{trialInd})
        for phzInd = 1:phzBin.count;
            for hvfInd = 1:hvfBin.count,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHvfPhzRmaps.rmap{trialInd}( abs(xbins)<=300,...
                                                         abs(ybins)<=300,...
                                                         unit,...
                                                         phzInd,...
                                                         hvfInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoHvfPhz.control.size   (ucounter,phzInd,hvfInd)   = sum(nniz(nratemap(:)));
                egoHvfPhz.control.maxRate(ucounter,phzInd,hvfInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHvfPhz.control.meanPos(ucounter,phzInd,hvfInd,:) = ratemapCenter./10 + offset * ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHvfPhz.control.peakPos(ucounter,phzInd,hvfInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHvfPhz.control.peakPos(ucounter,phzInd,hvfInd,:) = nan([1,1,1,2]);
                end
            end
        end
        ucounter = ucounter+1;
    end
end

ucounter = 1;
egoHvfPhz.shuffle.meanPos = []; % egoMeanRmapPosHvf
egoHvfPhz.shuffle.peakPos = []; % egoMaxRmapPosHvf
egoHvfPhz.shuffle.size = []; %egoSizeHvf = []; 
egoHvfPhz.shuffle.maxRate = []; %egoMeanRmapRateHvf = [];
%egoHvfPhz.shuffle.meanRate = []; %egoMaxRmapRateHvf = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHvfPhzRmaps.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgoHvf{trialInd})
        for phzInd = 1:phzBin.count
            for hvfInd = 1:hvfBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHvfPhzRmaps.rmapShuff{trialInd},6)
                    ratemap = sq(egoHvfPhzRmaps.rmapShuff{trialInd}( ...
                        abs(xbins)<=300,...
                        abs(ybins)<=300,...
                        unit,           ...
                        phzInd,         ...
                        hvfInd,         ...
                        iter));
                    nanmap = double(~isnan(ratemap));
                    nanmap(nanmap==0) = nan;
                    ratemap = ratemap.*fliplr(nanmap);
                    ratemap(ratemap<2) = 0;
                    nratemap =ratemap./sum(ratemap(:),'omitnan');
                    ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                    egoHvfPhz.shuffle.size   (ucounter,phzInd,hvfInd,iter)   = sum(nniz(nratemap(:)));
                    egoHvfPhz.shuffle.maxRate(ucounter,phzInd,hvfInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHvfPhz.shuffle.meanPos(ucounter,phzInd,hvfInd,iter,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHvfPhz.shuffle.peakPos(ucounter,phzInd,hvfInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHvfPhz.shuffle.peakPos(ucounter,phzInd,hvfInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hvfInd
        end% phzInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


% $$$ u = 76
% $$$ figure,histogram(sq(egoHvfPhz.shuffle.meanPos(u,1,3,:,2))- ...
% $$$                  sq(egoHvfPhz.shuffle.meanPos(u,1,2,:,2)),linspace([-10,10,30]));
% $$$ 
% $$$ figure,histogram((egoHvfPhz.control.meanPos(unitsEgoCA1,1,3,2)-mean(egoHvfPhz.shuffle.meanPos(unitsEgoCA1,1,3,:,2),4))...
% $$$     ./std(egoHvfPhz.shuffle.meanPos(unitsEgoCA1,1,3,:,2),[], ...
% $$$           4),linspace([-30,60,60]))





egoHvfPhz.perm.ca1.zscore = zeros([numel(unitsEgoHvfCA1),phzBin.count,hvfBin.count,2]);
for phzInd = 1:phzBin.count
    hvfi = 1;
    for hvfInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
            egoHvfPhz.perm.ca1.zscore(:,phzInd,hvfi,x) = ( (egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(1),x) ...
                                                            - egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(2),x))-(mean(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(1),:,x)-egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(2),:,x),4)))./std(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(1),:,x)-egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd(2),:,x),[],4);
        end
        hvfi = hvfi+1;
    end
end
egoHvfPhz.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA1)) ));


egoHvfPhz.boot.ca1.zscore = zeros([numel(unitsEgoHvfCA1),phzBin.count,hvfBin.count,2]);
for phzInd = 1:phzBin.count
    for hvfInd = 1:hvfBin.count
        for x = 1:2
            egoHvfPhz.boot.ca1.zscore(:,phzInd,hvfInd,x) = ...
                ( egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,hvfInd,x) ...
                  - mean(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd,:,x),4))./std(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA1,phzInd,hvfInd,:,x),[],4);
        end
    end
end
egoHvfPhz.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA1)) ));



egoHvfPhz.perm.ca3.zscore = zeros([numel(unitsEgoHvfCA3),phzBin.count,hvfBin.count,2]);
for phzInd = 1:phzBin.count
    hvfi = 1;
    for hvfInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHvfPhz.perm.ca3.zscore(:,phzInd,hvfi,x) = ( (egoHvfPhz.control.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(1),x) - egoHvfPhz.control.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(2),x))-(mean(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(1),:,x)-egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(2),:,x),4)))./std(egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(1),:,x)-egoHvfPhz.shuffle.meanPos(unitsEgoHvfCA3,phzInd,hvfInd(2),:,x),[],4);
        end
        hvfi = hvfi+1;
    end
end
egoHvfPhz.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA3)) ));

egoHvfPhz.xpos = egoHvfPhzRmaps.xpos;
egoHvfPhz.ypos = egoHvfPhzRmaps.ypos;
egoHvfPhz.xbins = egoHvfPhzRmaps.xbins;
egoHvfPhz.ybins = egoHvfPhzRmaps.ybins;
egoHvfPhz.rmap =  egoHvfPhzRmaps.rmap;


% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$ plot(1:4,sq(mean(egoHvfPhz.control.meanPos(unitsEgoHvfCA3,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$     plot(1:4,sq(mean(egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end

% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$     plot(1:4,sq(std(egoHvfPhz.control.meanPos(unitsEgoHvfCA1,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ for phzInd = 1:phzBin.count
% $$$     plot(1:4,sq(mean(egoHvfPhz.control.meanPos(unitsEgoHvfCA3,phzInd,:,1))),'-+','Color',phzBin.color(phzInd,:));
% $$$ end
% $$$ 
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ %histogram(egoHvfPhz.boot.ca1.zscore(:,1,4,1),linspace(-20,20,40))
% $$$ for phzInd = 1:phzBin.count,
% $$$     h = cdfplot(egoHvfPhz.boot.ca1.zscore(:,phzInd,4,1));
% $$$     h.Color = phzBin.color(phzInd,:);
% $$$ end
% $$$ Lines( egoHvfPhz.boot.ca1.sig,[],'r');
% $$$ Lines(-egoHvfPhz.boot.ca1.sig,[],'r');
% $$$ 
% $$$ 
% $$$ phzInd = 1;
% $$$ [h,p] = ttest2(egoHvfPhz.control.meanPos(unitsEgoHvfCA1(1:end),phzInd,2,1),...
% $$$               egoHvfPhz.control.meanPos(unitsEgoHvfCA1(1:end),phzInd,4,1))
% $$$ 
% $$$ 
% $$$ phzInd = 3;
% $$$ [h,p] = ttest2(egoHvfPhz.control.meanPos(unitsEgoHvfCA3(1:end),phzInd,2,1),...
% $$$                egoHvfPhz.control.meanPos(unitsEgoHvfCA3(1:end),phzInd,4,1))
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for phzInd = 1:phzBin.count
% $$$     subplot(3,1,phzInd)
% $$$     hold('on');
% $$$     for hvfInd = 2:hvfBin.count
% $$$         h = cdfplot(egoHvfPhz.control.meanPos(unitsEgoHvfCA1, ...
% $$$                                               phzInd,hvfInd,1));
% $$$         h.Color = hbaBin.color(hvfInd-1,:);
% $$$     end
% $$$ end
