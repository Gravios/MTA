%egoHbaPhzRmaps_pause = load(fullfile(Trials{1}.path.project,'analysis','EgoProCode2D_compute_egoratemaps_conditioned_on_hba_and_phz_DATA_pause.mat'));

egoHbaRmaps_pause    = load(fullfile(Trials{1}.path.project,...
                            'analysis',...
                            'EgoProCode2D_compute_egoratemaps_conditioned_on_hba_DATA_pause.mat'));


xbins = egoHbaPhzRmaps_pause.xbins;
ybins = egoHbaPhzRmaps_pause.ybins;

ucounter = 1;
egoHbaPhzPause.control.meanPos = [];   %egoMeanRmapPosHba
egoHbaPhzPause.control.peakPos = [];   %egoMaxRmapPosHba
egoHbaPhzPause.control.size = [];      %egoSizeHba = []; 
egoHbaPhzPause.control.maxRate = []    %egoMeanRmapRateHba = [];
%egoHbaPhzPause.control.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps_pause.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgo{trialInd})
        for phzInd = 1:3;
            for hbaInd = 1:3,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHbaPhzRmaps_pause.rmap{trialInd}( abs(xbins)<=300,...
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
                egoHbaPhzPause.control.size   (ucounter,phzInd,hbaInd)   = sum(nniz(nratemap(:)));
                egoHbaPhzPause.control.maxRate(ucounter,phzInd,hbaInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHbaPhzPause.control.meanPos(ucounter,phzInd,hbaInd,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHbaPhzPause.control.peakPos(ucounter,phzInd,hbaInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHbaPhzPause.control.peakPos(ucounter,phzInd,hbaInd,:) = nan([1,1,1,2]);
                end
            end
        end
        ucounter = ucounter+1;
    end
end

ucounter = 1;
egoHbaPhzPause.shuffle.meanPos = []; % egoMeanRmapPosHba
egoHbaPhzPause.shuffle.peakPos = []; % egoMaxRmapPosHba
egoHbaPhzPause.shuffle.size = []; %egoSizeHba = []; 
egoHbaPhzPause.shuffle.maxRate = [] %egoMeanRmapRateHba = [];
%egoHbaPhzPause.shuffle.meanRate = []; %egoMaxRmapRateHba = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHbaPhzRmaps_pause.rmapShuff{trialInd})
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
                for iter = 1:size(egoHbaPhzRmaps_pause.rmapShuff{trialInd},6)
                    ratemap = sq(egoHbaPhzRmaps_pause.rmapShuff{trialInd}( ...
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
                    egoHbaPhzPause.shuffle.size   (ucounter,phzInd,hbaInd,iter)   = sum(nniz(nratemap(:)));
                    egoHbaPhzPause.shuffle.maxRate(ucounter,phzInd,hbaInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHbaPhzPause.shuffle.meanPos(ucounter,phzInd,hbaInd,iter,:) = ratemapCenter./10;%+ [-2,1.6]*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHbaPhzPause.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHbaPhzPause.shuffle.peakPos(ucounter,phzInd,hbaInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hbaInd
        end% phzInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


% $$$ u = 76
% $$$ figure,histogram(sq(egoHbaPhzPause.shuffle.meanPos(u,1,3,:,2))- ...
% $$$                  sq(egoHbaPhzPause.shuffle.meanPos(u,1,2,:,2)),linspace([-10,10,30]));
% $$$ 
% $$$ figure,histogram((egoHbaPhzPause.control.meanPos(unitsEgoCA1,1,3,2)-mean(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,1,3,:,2),4))...
% $$$     ./std(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,1,3,:,2),[], ...
% $$$           4),linspace([-30,60,60]))


egoHbaPhzPause.perm.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHbaPhzPause.perm.ca1.zscore(:,phzInd,hbai,x) = ( (egoHbaPhzPause.control.meanPos(unitsEgoCA1,phzInd,hbaInd(1),x) - egoHbaPhzPause.control.meanPos(unitsEgoCA1,phzInd,hbaInd(2),x))-(mean(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(1),:,x)-egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
egoHbaPhzPause.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));

egoHbaPhzPause.perm.ca3.zscore = zeros([numel(unitsEgoCA3),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    hbai = 1;
    for hbaInd = [3,1,3;...
                  2,2,1]
        for x = 1:2
egoHbaPhzPause.perm.ca3.zscore(:,phzInd,hbai,x) = ( (egoHbaPhzPause.control.meanPos(unitsEgoCA3,phzInd,hbaInd(1),x) - egoHbaPhzPause.control.meanPos(unitsEgoCA3,phzInd,hbaInd(2),x))-(mean(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhzPause.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),4)))./std(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(1),:,x)-egoHbaPhzPause.shuffle.meanPos(unitsEgoCA3,phzInd,hbaInd(2),:,x),[],4);
        end
        hbai = hbai+1;
    end
end
egoHbaPhzPause.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA3)) ));




 
egoHbaPhzPause.boot.ca1.zscore = zeros([numel(unitsEgoCA1),phzBin.count,hbaBin.count,2]);
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        for x = 1:2
            meanPos      = egoHbaPhzPause.control.meanPos(unitsEgoCA1,phzInd,hbaInd,x);
            meanPosShuff = mean(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),4);
            stdPosShuff  = std(egoHbaPhzPause.shuffle.meanPos(unitsEgoCA1,phzInd,hbaInd,:,x),[],4);
            
            zscores      = ( meanPos - meanPosShuff ) ./ stdPosShuff;
            
            egoHbaPhzPause.boot.ca1.zscore(:,phzInd,hbaInd,x) = zscores
        end
    end
end
egoHbaPhzPause.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoCA1)) ));


% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:phzBin.count
% $$$     subplot2(3, 1, 4-phzInd, 1);
% $$$     hold('on');
% $$$ for hbaInd = 1:hbaBin.count
% $$$     if hbaInd == 2
% $$$         cdfplot(-egoHbaPhzLoc.perm.ca1.zscore(:,phzInd,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaPhzLoc.perm.ca1.zscore(:,phzInd,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaPhzLoc.perm.ca1.sig,[],'r');
% $$$ Lines(-egoHbaPhzLoc.perm.ca1.sig,[],'r');
% $$$ xlim([-20,30]);
% $$$ legend({'R-C','L-C','R-L'},'Location','southeast');
% $$$ xlabel('z-score');
% $$$ end
% $$$ 
% $$$ 
% $$$ figure();
% $$$ for phzInd = 1:phzBin.count
% $$$     subplot2(3, 1, 4-phzInd, 1);
% $$$     hold('on');
% $$$ for hbaInd = 1:hbaBin.count
% $$$     if hbaInd == 2
% $$$         cdfplot(-egoHbaPhzPause.perm.ca1.zscore(:,phzInd,hbaInd,2));
% $$$     else            
% $$$         cdfplot(egoHbaPhzPause.perm.ca1.zscore(:,phzInd,hbaInd,2));
% $$$     end
% $$$ end
% $$$ Lines(egoHbaPhzPause.perm.ca1.sig,[],'r');
% $$$ Lines(-egoHbaPhzPause.perm.ca1.sig,[],'r');
% $$$ xlim([-20,30]);
% $$$ legend({'R-C','L-C','R-L'},'Location','southeast');
% $$$ xlabel('z-score');
% $$$ end