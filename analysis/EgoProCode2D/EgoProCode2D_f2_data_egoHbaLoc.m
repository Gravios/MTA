

%egoHbaRmaps_loc

xbins = egoHbaRmaps_loc.xbins;
ybins = egoHbaRmaps_loc.ybins;
ucounter = 1;
egoHbaLoc.control.meanPos = []; % egoMeanRmapPosHba
egoHbaLoc.control.peakPos = []; % egoMaxRmapPosHba
egoHbaLoc.control.size = []; %egoSizeHba = []; 
egoHbaLoc.control.maxRate = []; %egoMeanRmapRateHba = [];
%egoHbaLoc.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps_loc.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaRmaps_loc.rmap{trialInd}( abs(xbins)<=300,...
                                          abs(ybins)<=300,...
                                          unit,...
                                          hbaInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoHbaLoc.control.size   (ucounter,hbaInd)   = sum(nniz(nratemap(:)));
                egoHbaLoc.control.maxRate(ucounter,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHbaLoc.control.meanPos(ucounter,hbaInd,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHbaLoc.control.peakPos(ucounter,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHbaLoc.control.peakPos(ucounter,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        ucounter = ucounter+1;
    end
end



ucounter = 1;
egoHbaLoc.shuffle.meanPos = [];   % egoMeanRmapPosHba
egoHbaLoc.shuffle.peakPos = [];   % egoMaxRmapPosHba
egoHbaLoc.shuffle.size = [];      % egoSizeHba = []; 
egoHbaLoc.shuffle.maxRate = []    % egoMeanRmapRateHba = [];
%egoHbaLoc.shuffle.meanRate = []; % egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaRmaps_loc.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
            for hbaInd = 1:hbaBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHbaRmaps_loc.rmapShuff{trialInd},5)
                    ratemap = sq(egoHbaRmaps_loc.rmapShuff{trialInd}( ...
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
                    egoHbaLoc.shuffle.size   (ucounter,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHbaLoc.shuffle.maxRate(ucounter,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHbaLoc.shuffle.meanPos(ucounter,hbaInd,iter,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHbaLoc.shuffle.peakPos(ucounter,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHbaLoc.shuffle.peakPos(ucounter,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


egoHbaLoc.perm.ca1.zscore = zeros([numel(unitsEgoCA1),hbaBin.count,2]);

hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHbaLoc.perm.ca1.zscore(:,hbai,x) = ( (egoHbaLoc.control.meanPos(unitsEgoCA1,hbaInd(1),x) - egoHbaLoc.control.meanPos(unitsEgoCA1,hbaInd(2),x))-(mean(egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),3)))./std(egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd(1),:,x)-egoHbaLoc.shuffle.meanPos(unitsEgoCA1,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHbaLoc.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));


egoHbaLoc.perm.ca3.zscore = zeros([numel(unitsEgoCA3),hbaBin.count,2]);
hbai = 1;
for hbaInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHbaLoc.perm.ca3.zscore(:,hbai,x) = ( (egoHbaLoc.control.meanPos(unitsEgoCA3,hbaInd(1),x) - egoHbaLoc.control.meanPos(unitsEgoCA3,hbaInd(2),x))-(mean(egoHbaLoc.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHbaLoc.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),3)))./std(egoHbaLoc.shuffle.meanPos(unitsEgoCA3,hbaInd(1),:,x)-egoHbaLoc.shuffle.meanPos(unitsEgoCA3,hbaInd(2),:,x),[],3);
    end
    hbai = hbai+1;
end
egoHbaLoc.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));





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
% $$$ figure();
% $$$ hold('on');
% $$$ subplot2(1, 1, 1, 1);
% $$$ hold('on')
% $$$ for hbaInd = 1:hbaBin.count
% $$$     if hbaInd == 2
% $$$         cdfplot(-egoHbaLoc.perm.ca1.zscore(:,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaLoc.perm.ca1.zscore(:,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaLoc.perm.ca1.sig,[],'r');
% $$$ Lines(-egoHbaLoc.perm.ca1.sig,[],'r');
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
% $$$         cdfplot(-egoHbaLoc.perm.ca3.zscore(:,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaLoc.perm.ca3.zscore(:,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaLoc.perm.ca3.sig,[],'r');
% $$$ Lines(-egoHbaLoc.perm.ca3.sig,[],'r');
% $$$ xlim([-20,30]);
% $$$ legend({'R-C','L-C','R-L'},'Location','southeast');
% $$$ xlabel('z-score');
