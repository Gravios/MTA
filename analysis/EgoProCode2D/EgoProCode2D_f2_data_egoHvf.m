% EgoProCode2D_f2_data_egoHvf.m
egoHvfRmaps = load(fullfile(Trials{1}.path.project,...
              'analysis',...
              'EgoProCode2D_compute_egoratemaps_conditioned_on_hvf_DATA.mat'));



xbins = egoHvfRmaps.xbins;
ybins = egoHvfRmaps.ybins;
ucounter = 1;
egoHvf.control.meanPos = []; % egoMeanRmapPosHvf
egoHvf.control.peakPos = []; % egoMaxRmapPosHvf
egoHvf.control.size = []; %egoSizeHvf = []; 
egoHvf.control.maxRate = []; %egoMeanRmapRateHvf = [];
%egoHvf.control.meanRate = []; %egoMaxRmapRateHvf = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHvfRmaps.rmap{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgoHvf{trialInd})
            for hvfInd = 1:hvfBin.count,
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1,mapPosition{:});
                ratemap = egoHvfRmaps.rmap{trialInd}( abs(xbins)<=300,...
                                          abs(ybins)<=300,...
                                          unit,...
                                          hvfInd);
                nanmap = double(~isnan(ratemap));
                nanmap(nanmap==0) = nan;
                ratemap = ratemap.*fliplr(nanmap);
                ratemap(ratemap<2) = 0;
                nratemap =ratemap./sum(ratemap(:),'omitnan');
                ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                egoHvf.control.size   (ucounter,hvfInd)   = sum(nniz(nratemap(:)));
                egoHvf.control.maxRate(ucounter,hvfInd,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                egoHvf.control.meanPos(ucounter,hvfInd,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                [~,maxPos] = max(nratemap(:));
                if ~isempty(maxPos)
                    [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                    egoHvf.control.peakPos(ucounter,hvfInd,:) = mapPosition(maxX,maxY,:);
                else
                    egoHvf.control.peakPos(ucounter,hvfInd,:) = nan([1,1,1,2]);
                end
            end
        ucounter = ucounter+1;
    end
end



ucounter = 1;
egoHvf.shuffle.meanPos = [];   % egoMeanRmapPosHvf
egoHvf.shuffle.peakPos = [];   % egoMaxRmapPosHvf
egoHvf.shuffle.size = [];      % egoSizeHvf = []; 
egoHvf.shuffle.maxRate = [];    % egoMeanRmapRateHvf = [];
%egoHvf.shuffle.meanRate = []; % egoMaxRmapRateHvf = [];

for trialInd = 1:numel(Trials)
    if isempty(egoHvfRmaps.rmapShuff{trialInd})
        continue;
    end
    for unit = 1:numel(unitsEgoHvf{trialInd})
            for hvfInd = 1:hvfBin.count
                mapPosition = cell([1,2]);
                [mapPosition{:}] = ndgrid(xbins(abs(xbins)<=300),...
                                          ybins(abs(ybins)<=300));
                mapPosition = cat(numel(mapPosition)+1, ...
                                  mapPosition{:});
                for iter = 1:size(egoHvfRmaps.rmapShuff{trialInd},5)
                    ratemap = sq(egoHvfRmaps.rmapShuff{trialInd}( ...
                        abs(xbins)<=300,...
                        abs(ybins)<=300,...
                        unit,           ...
                        hvfInd,         ...
                        iter));
                    nanmap = double(~isnan(ratemap));
                    nanmap(nanmap==0) = nan;
                    ratemap = ratemap.*fliplr(nanmap);
                    ratemap(ratemap<2) = 0;
                    nratemap =ratemap./sum(ratemap(:),'omitnan');
                    ratemapCenter = sq(sum(sum(bsxfun(@times,nratemap,mapPosition),'omitnan'),'omitnan'))';
                    egoHvf.shuffle.size   (ucounter,hvfInd,iter)   = sum(nniz(nratemap(:)));
                    egoHvf.shuffle.maxRate(ucounter,hvfInd,iter,:) = mean(ratemap(nniz(ratemap(:))),'omitnan');
                    egoHvf.shuffle.meanPos(ucounter,hvfInd,iter,:) = ratemapCenter./10+ offset*ismember(trialInd,[3:5,18:26,30]);
                    [~,maxPos] = max(nratemap(:));
                    if ~isempty(maxPos)
                        [maxX,maxY] = ind2sub(size(nratemap),maxPos);
                        egoHvf.shuffle.peakPos(ucounter,hvfInd,iter,:) = mapPosition(maxX,maxY,:);
                    else
                        egoHvf.shuffle.peakPos(ucounter,hvfInd,iter,:) = nan([1,1,1,1,2]);
                    end
                end% iter
            end% hvfInd
        ucounter = ucounter+1;
    end% unit
end% trialInd


egoHvf.perm.ca1.zscore = zeros([numel(unitsEgoHvfCA1),hvfBin.count,2]);

hvfi = 1;
for hvfInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHvf.perm.ca1.zscore(:,hvfi,x) = ( (egoHvf.control.meanPos(unitsEgoHvfCA1,hvfInd(1),x) - egoHvf.control.meanPos(unitsEgoHvfCA1,hvfInd(2),x))-(mean(egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd(1),:,x)-egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd(2),:,x),3)))./std(egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd(1),:,x)-egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd(2),:,x),[],3);
    end
    hvfi = hvfi+1;
end
egoHvf.perm.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA1)) ));


egoHvf.perm.ca3.zscore = zeros([numel(unitsEgoHvfCA3),hvfBin.count,2]);
hvfi = 1;
for hvfInd = [3,1,3;...
              2,2,1]
    for x = 1:2
        egoHvf.perm.ca3.zscore(:,hvfi,x) = ( (egoHvf.control.meanPos(unitsEgoHvfCA3,hvfInd(1),x) - egoHvf.control.meanPos(unitsEgoHvfCA3,hvfInd(2),x))-(mean(egoHvf.shuffle.meanPos(unitsEgoHvfCA3,hvfInd(1),:,x)-egoHvf.shuffle.meanPos(unitsEgoHvfCA3,hvfInd(2),:,x),3)))./std(egoHvf.shuffle.meanPos(unitsEgoHvfCA3,hvfInd(1),:,x)-egoHvf.shuffle.meanPos(unitsEgoHvfCA3,hvfInd(2),:,x),[],3);
    end
    hvfi = hvfi+1;
end
egoHvf.perm.ca3.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA3)) ));



egoHvf.boot.ca1.zscore = zeros([numel(unitsEgoHvfCA1),hvfBin.count,2]);

for hvfInd = 1:hvfBin.count
    for x = 1:2
        egoHvf.boot.ca1.zscore(:,hvfInd,x) = ...
            ( egoHvf.control.meanPos(unitsEgoHvfCA1,hvfInd,x) ...
              - mean(egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd,:,x),3) ...
              )./std(egoHvf.shuffle.meanPos(unitsEgoHvfCA1,hvfInd,:,x),[],3);
    end
end

egoHvf.boot.ca1.sig = abs(norminv(1-(1-0.05)^(1/numel(unitsEgoHvfCA1)) ));


egoHvf.xpos = egoHvfRmaps.xpos;
egoHvf.ypos = egoHvfRmaps.ypos;
egoHvf.xbins = egoHvfRmaps.xbins;
egoHvf.ybins = egoHvfRmaps.ybins;
egoHvf.rmap =  egoHvfRmaps.rmap;

% $$$ 
% $$$ figure,
% $$$ hold('on');
% $$$ cdfplot(egoHvf.boot.ca1.zscore(:,4,1))
% $$$ 
% $$$ 
% $$$ figure
% $$$ plot(egoHvf.control.meanPos(unitsEgoHvfCA1,2,1),egoHvf.control.meanPos(unitsEgoHvfCA1,3,1),'.');
% $$$ daspect([1,1,1]);
% $$$ ylim([-10,10]);
% $$$ xlim([-10,10]);
% $$$ hold('on');
% $$$ grid('on');
% $$$ line([-10,10],[-10,10]);
% $$$ 
% $$$ figure
% $$$ plot(egoHvf.control.meanPos(unitsEgoHvfCA3,2,1),egoHvf.control.meanPos(unitsEgoHvfCA3,4,1),'.');
% $$$ daspect([1,1,1]);
% $$$ ylim([-5,15]);
% $$$ xlim([-5,15]);
% $$$ hold('on');
% $$$ grid('on');
% $$$ line([-5,15],[-5,15]);


