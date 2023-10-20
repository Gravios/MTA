
egoHbaPhzRmaps_loc   = load(fullfile(Trials{1}.path.project,...
                            'analysis',...
                            'EgoProCode2D_compute_egoratemaps_conditioned_on_hba_and_phz_DATA_loc.mat'));



xbins = egoHbaPhzRmaps_loc.xbins;
ybins = egoHbaPhzRmaps_loc.ybins;

ucounter = 1;
egoHbaPhzLoc.control.meanPos = [];   %egoMeanRmapPosHba
egoHbaPhzLoc.control.peakPos = [];   %egoMaxRmapPosHba
egoHbaPhzLoc.control.size = [];      %egoSizeHba = []; 
egoHbaPhzLoc.control.maxRate = []    %egoMeanRmapRateHba = [];
%egoHbaPhzLoc.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps_loc.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:3;
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaPhzRmaps_loc.rmap{trialInd}( abs(xbins)<=300,...
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
                egoHbaPhzLoc.control.size   (ucounter,phzInd,hbaInd)   = sum(nniz(nratemap(:)));
                egoHbaPhzLoc.control.maxRate(ucounter,phzInd,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHbaPhzLoc.control.meanPos(ucounter,phzInd,hbaInd,:) = ratemapCenter./10 + offset*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHbaPhzLoc.control.peakPos(ucounter,phzInd,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHbaPhzLoc.control.peakPos(ucounter,phzInd,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        end
        ucounter = ucounter+1;
    end
end

ucounter = 1;
egoHbaPhzLoc.shuffle.meanPos = []; % egoMeanRmapPosHba
egoHbaPhzLoc.shuffle.peakPos = []; % egoMaxRmapPosHba
egoHbaPhzLoc.shuffle.size = []; %egoSizeHba = []; 
egoHbaPhzLoc.shuffle.maxRate = [] %egoMeanRmapRateHba = [];
%egoHbaPhzLoc.shuffle.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps_loc.rmapShuff{trialInd})
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
                for iter = 1:size(egoHbaPhzRmaps_loc.rmapShuff{trialInd},6)
                    ratemap = sq(egoHbaPhzRmaps_loc.rmapShuff{trialInd}( ...
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
                    egoHbaPhzLoc.shuffle.size   (ucounter,phzInd,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHbaPhzLoc.shuffle.maxRate(ucounter,phzInd,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHbaPhzLoc.shuffle.meanPos(ucounter,phzInd,hbaInd,iter,:) = ratemapCenter./10 + offset*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHbaPhzLoc.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHbaPhzLoc.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        end% phzInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


% $$$ u = 76
% $$$ figure,histogram(sq(egoHbaPhzLoc.shuffle.meanPos(u,1,3,:,2))- ...
% $$$                  sq(egoHbaPhzLoc.shuffle.meanPos(u,1,2,:,2)),linspace([-10,10,30]));
% $$$ 
% $$$ figure,histogram((egoHbaPhzLoc.control.meanPos(unitsEgoCA1,1,3,2)-mean(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,1,3,:,2),4))...
% $$$     ./std(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,1,3,:,2),[], ...
% $$$           4),linspace([-30,60,60]))


egoHbaPhzLoc.perm.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHbaPhzLoc.perm.ca1.zscore(:,phzInd,hbai,x) = ( (egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd(1),x) - egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd(2),x))-(mean(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
egoHbaPhzLoc.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));

egoHbaPhzLoc.perm.ca3.zscore = zeros([numel(unitsEgoCA3),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHbaPhzLoc.perm.ca3.zscore(:,phzInd,hbai,x) = ( (egoHbaPhzLoc.control.meanPos(unitsEgoCA3,phzInd,hbaInd(1),x) - egoHbaPhzLoc.control.meanPos(unitsEgoCA3,phzInd,hbaInd(2),x))-(mean(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
 egoHbaPhzLoc.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));



 
egoHbaPhzLoc.boot.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        for x = 1:2
            egoHbaPhzLoc.boot.ca1.zscore(:,phzInd,hbaInd,x) = ...
                ( egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd,x) ...
                  - mean(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),4))./std(egoHbaPhzLoc.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),[],4);
        end
    end
end
egoHbaPhzLoc.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));
