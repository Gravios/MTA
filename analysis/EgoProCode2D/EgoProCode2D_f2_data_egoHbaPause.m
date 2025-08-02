

%egoHbaRmaps_pause

xbins = egoHbaRmaps_pause.xbins;
ybins = egoHbaRmaps_pause.ybins;
ucounter = 1;
egoHbaPause.control.meanPos = []; % egoMeanRmapPosHba
egoHbaPause.control.peakPos = []; % egoMaxRmapPosHba
egoHbaPause.control.size = []; %egoSizeHba = []; 
egoHbaPause.control.maxRate = []; %egoMeanRmapRateHba = [];
%egoHbaPause.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps_pause.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaRmaps_pause.rmap{trialInd}( abs(xbins)<=300,...
                                          abs(ybins)<=300,...
                                          unit,...
                                          hbaInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoHbaPause.control.size   (ucounter,hbaInd)   = sum(nniz(nratemap(:)));
                egoHbaPause.control.maxRate(ucounter,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHbaPause.control.meanPos(ucounter,hbaInd,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHbaPause.control.peakPos(ucounter,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHbaPause.control.peakPos(ucounter,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        ucounter = ucounter+1;
    end
end



ucounter = 1;
egoHbaPause.shuffle.meanPos = [];   % egoMeanRmapPosHba
egoHbaPause.shuffle.peakPos = [];   % egoMaxRmapPosHba
egoHbaPause.shuffle.size = [];      % egoSizeHba = []; 
egoHbaPause.shuffle.maxRate = []    % egoMeanRmapRateHba = [];
%egoHbaPause.shuffle.meanRate = []; % egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps_pause.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:hbaBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHbaRmaps_pause.rmapShuff{trialInd},5)
                    ratemap = sq(egoHbaRmaps_pause.rmapShuff{trialInd}( ...
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
                    egoHbaPause.shuffle.size   (ucounter,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHbaPause.shuffle.maxRate(ucounter,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHbaPause.shuffle.meanPos(ucounter,hbaInd,iter,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHbaPause.shuffle.peakPos(ucounter,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHbaPause.shuffle.peakPos(ucounter,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


egoHbaPause.perm.ca1.zscore = zeros([numel(unitsEgoCA1),hbaBin.count,2]);

hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHbaPause.perm.ca1.zscore(:,hbai,x) = ( (egoHbaPause.control.meanPos(unitsEgoCA1,hbaInd(1),x) - egoHbaPause.control.meanPos(unitsEgoCA1,hbaInd(2),x))-(mean(egoHbaPause.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHbaPause.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),3)))./std(egoHbaPause.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHbaPause.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHbaPause.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));


egoHbaPause.perm.ca3.zscore = zeros([numel(unitsEgoCA3),hbaBin.count,2]);
hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHbaPause.perm.ca3.zscore(:,hbai,x) = ( (egoHbaPause.control.meanPos(unitsEgoCA3,hbaInd(1),x) - egoHbaPause.control.meanPos(unitsEgoCA3,hbaInd(2),x))-(mean(egoHbaPause.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHbaPause.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),3)))./std(egoHbaPause.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHbaPause.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHbaPause.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));



egoHbaLoc.boot.ca1.zscore = zeros([numel(unitsEgoCA1),hbaBin.count,2]);
    for hbaInd = 1:hbaBin.count
        for x = 1:2
            egoHbaLoc.boot.ca1.zscore(:,hbaInd,x) = ...
                ( egoHbaLoc.control.meanPos(unitsEgoCA1,hbaInd,x) ...
                  - mean(egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd,:,x),3))./std(egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd,:,x),[],3);
        end
    end
egoHbaLoc.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));
% $$$ 
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ subplot2(1, 1, 1, 1);
% $$$ hold('on')
% $$$ for hbaInd = 1:hbaBin.count
% $$$     if hbaInd == 2
% $$$         cdfplot(-egoHbaPause.perm.ca1.zscore(:,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaPause.perm.ca1.zscore(:,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaPause.perm.ca1.sig,[],'r');
% $$$ Lines(-egoHbaPause.perm.ca1.sig,[],'r');
% $$$ xlim([-20,30]);
% $$$ legend({'R-C','L-C','R-L'},'Location','southeast');
% $$$ xlabel('z-score');
% $$$ 
% $$$ 
% $$$ figure();
% $$$ hold('on');
% $$$ subplot2(1, 1, 1, 1);
% $$$ hold('on')
% $$$ for hbaInd = 1:hbaBin.count
% $$$     if hbaInd == 2
% $$$         cdfplot(-egoHbaPause.perm.ca3.zscore(:,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaPause.perm.ca3.zscore(:,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaPause.perm.ca3.sig,[],'r');
% $$$ Lines(-egoHbaPause.perm.ca3.sig,[],'r');
% $$$ xlim([-20,30]);
% $$$ legend({'R-C','L-C','R-L'},'Location','southeast');
% $$$ xlabel('z-score');
