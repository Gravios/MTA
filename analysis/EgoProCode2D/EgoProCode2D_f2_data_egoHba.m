egoHbaRmaps = load(fullfile(Trials{1}.path.project,...
              'analysis',...
              'EgoProCode2D_compute_egoratemaps_conditioned_on_hba_DATA.mat'));



xbins = egoHbaRmaps.xbins;
ybins = egoHbaRmaps.ybins;
ucounter = 1;
egoHba.control.meanPos = []; % egoMeanRmapPosHba
egoHba.control.peakPos = []; % egoMaxRmapPosHba
egoHba.control.size = []; %egoSizeHba = []; 
egoHba.control.maxRate = [] %egoMeanRmapRateHba = [];
%egoHba.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaRmaps.rmap{trialInd}( abs(xbins)<=300,...
                                          abs(ybins)<=300,...
                                          unit,...
                                          hbaInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoHba.control.size   (ucounter,hbaInd)   = sum(nniz(nratemap(:)));
                egoHba.control.maxRate(ucounter,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHba.control.meanPos(ucounter,hbaInd,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHba.control.peakPos(ucounter,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHba.control.peakPos(ucounter,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        ucounter = ucounter+1;
    end
end



ucounter = 1;
egoHba.shuffle.meanPos = [];   % egoMeanRmapPosHba
egoHba.shuffle.peakPos = [];   % egoMaxRmapPosHba
egoHba.shuffle.size = [];      % egoSizeHba = []; 
egoHba.shuffle.maxRate = []    % egoMeanRmapRateHba = [];
%egoHba.shuffle.meanRate = []; % egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:hbaBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHbaRmaps.rmapShuff{trialInd},5)
                    ratemap = sq(egoHbaRmaps.rmapShuff{trialInd}( ...
                        abs(xbins)<=300,...
                        abs(ybins)<=300,...
                        unit,           ...
                        hbaInd,         ...
                        iter));
                    nanmap = double(~isnan(ratemap));
                    nanmap(nanmap==0) = nan;
                    ratemap = ratemap.*fliplr(nanmap);
                    ratemap(ratemap<2) = 0;
                    nratemap =ratemap./sum(ratemap(:),'omitnan');
                    ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                    egoHba.shuffle.size   (ucounter,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHba.shuffle.maxRate(ucounter,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHba.shuffle.meanPos(ucounter,hbaInd,iter,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHba.shuffle.peakPos(ucounter,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHba.shuffle.peakPos(ucounter,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


egoHba.perm.ca1.zscore = zeros([numel(unitsEgoCA1),hbaBin.count,2]);

hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHba.perm.ca1.zscore(:,hbai,x) = ( (egoHba.control.meanPos(unitsEgoCA1,hbaInd(1),x) - egoHba.control.meanPos(unitsEgoCA1,hbaInd(2),x))-(mean(egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),3)))./std(egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHba.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));


egoHba.perm.ca3.zscore = zeros([numel(unitsEgoCA3),hbaBin.count,2]);
hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHba.perm.ca3.zscore(:,hbai,x) = ( (egoHba.control.meanPos(unitsEgoCA3,hbaInd(1),x) - egoHba.control.meanPos(unitsEgoCA3,hbaInd(2),x))-(mean(egoHba.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHba.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),3)))./std(egoHba.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHba.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHba.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));



egoHba.boot.ca1.zscore = zeros([numel(unitsEgoCA1),hbaBin.count,2]);

for hbaInd = 1:hbaBin.count
    for x = 1:2
        egoHba.boot.ca1.zscore(:,hbaInd,x) = ...
            ( egoHba.control.meanPos(unitsEgoCA1,hbaInd,x) ...
              - mean(egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd,:,x),3) ...
              )./std(egoHba.shuffle.meanPos(unitsEgoCA1,hbaInd,:,x),[],3);
    end
end

egoHba.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));
