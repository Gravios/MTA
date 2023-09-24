
egoHbaPhzRmaps = load(fullfile(Trials{1}.path.project,'analysis','EgoProCode2D_compute_egoratemaps_conditioned_on_hba_and_phz_DATA.mat'));


xbins = egoHbaPhzRmaps.xbins;
ybins = egoHbaPhzRmaps.ybins;

ucounter = 1;
egoHbaPhz.control.meanPos = [];   %egoMeanRmapPosHba
egoHbaPhz.control.peakPos = [];   %egoMaxRmapPosHba
egoHbaPhz.control.size = [];      %egoSizeHba = []; 
egoHbaPhz.control.maxRate = []    %egoMeanRmapRateHba = [];
%egoHbaPhz.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:3;
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaPhzRmaps.rmap{trialInd}( abs(xbins)<=300,...
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
                egoHbaPhz.control.size   (ucounter,phzInd,hbaInd)   = sum(nniz(nratemap(:)));
                egoHbaPhz.control.maxRate(ucounter,phzInd,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHbaPhz.control.meanPos(ucounter,phzInd,hbaInd,:) = ratemapCenter./10 + offset * ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHbaPhz.control.peakPos(ucounter,phzInd,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHbaPhz.control.peakPos(ucounter,phzInd,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        end
        ucounter = ucounter+1;
    end
end

ucounter = 1;
egoHbaPhz.shuffle.meanPos = []; % egoMeanRmapPosHba
egoHbaPhz.shuffle.peakPos = []; % egoMaxRmapPosHba
egoHbaPhz.shuffle.size = []; %egoSizeHba = []; 
egoHbaPhz.shuffle.maxRate = [] %egoMeanRmapRateHba = [];
%egoHbaPhz.shuffle.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:phzBin.count
            for hbaInd = 1:hbaBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHbaPhzRmaps.rmapShuff{trialInd},6)
                    ratemap = sq(egoHbaPhzRmaps.rmapShuff{trialInd}( ...
                        abs(xbins)<=300,...
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
                    egoHbaPhz.shuffle.size   (ucounter,phzInd,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHbaPhz.shuffle.maxRate(ucounter,phzInd,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHbaPhz.shuffle.meanPos(ucounter,phzInd,hbaInd,iter,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHbaPhz.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHbaPhz.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        end% phzInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


% $$$ u = 76
% $$$ figure,histogram(sq(egoHbaPhz.shuffle.meanPos(u,1,3,:,2))- ...
% $$$                  sq(egoHbaPhz.shuffle.meanPos(u,1,2,:,2)),linspace([-10,10,30]));
% $$$ 
% $$$ figure,histogram((egoHbaPhz.control.meanPos(unitsEgoCA1,1,3,2)-mean(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,1,3,:,2),4))...
% $$$     ./std(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,1,3,:,2),[], ...
% $$$           4),linspace([-30,60,60]))


egoHbaPhz.perm.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
            egoHbaPhz.perm.ca1.zscore(:,phzInd,hbai,x) = ( (egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd(1),x) ...
                                                            - egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd(2),x))-(mean(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
egoHbaPhz.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));


egoHbaPhz.boot.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        for x = 1:2
            egoHbaPhz.boot.ca1.zscore(:,phzInd,hbaInd,x) = ...
                ( egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,x) ...
                  - mean(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),4))./std(egoHbaPhz.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),[],4);
        end
    end
end
egoHbaPhz.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));



egoHbaPhz.perm.ca3.zscore = zeros([numel(unitsEgoCA3),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHbaPhz.perm.ca3.zscore(:,phzInd,hbai,x) = ( (egoHbaPhz.control.meanPos(unitsEgoCA3,phzInd,hbaInd(1),x) - egoHbaPhz.control.meanPos(unitsEgoCA3,phzInd,hbaInd(2),x))-(mean(egoHbaPhz.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhz.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhz.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhz.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
egoHbaPhz.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));


