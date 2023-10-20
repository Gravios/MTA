

configure_default_args();
EgoProCode2D_load_data(); 
EgoProCode2D_f2_data_egoHba()
EgoProCode2D_f2_data_egoHbaPhz()







figure();
subplot(131);
hold('on');
plot(egoHbaPause.perm.ca1.zscore(:,3,2),egoHbaLoc.perm.ca1.zscore(:,3,2),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Left Vs Right'});
xlabel('z-score (pause)');
ylabel('z-score (loc)');
daspect([1,1,1]);

subplot(132);
hold('on');
plot(egoHbaPause.perm.ca1.zscore(:,1,2),egoHbaLoc.perm.ca1.zscore(:,1,2),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Right Vs Center'})
xlabel('z-score (pause)');
ylabel('z-score (loc)');
daspect([1,1,1]);

subplot(133);
hold('on');
plot(-egoHbaPause.perm.ca1.zscore(:,2,2),-egoHbaLoc.perm.ca1.zscore(:,2,2),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
circle(0,0,egoHba.perm.ca1.sig,'r-');
linkax('xy');
title({'Permutation Test','HBA - Left Vs Center'});
xlabel('z-score (pause)');
ylabel('z-score (loc)');
daspect([1,1,1]);



figpath = '/storage/share/Projects/EgoProCode2D/'
figure();
comp = 2;
dim = 1; dimlabel = 'forward'
dim = 2; dimlabel = 'lateral'
subplot(121);
hold('on');
plot(egoHbaPause.perm.ca1.zscore(:,comp,dim),egoHbaLoc.perm.ca1.zscore(:,comp,dim),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Left Vs Right'});
xlabel([dimlabel,' z-score (pause)']);
ylabel([dimlabel,' z-score (loc)']);
daspect([1,1,1]);

subplot(122);
hold('on');
plot(egoHbaPause.control.meanPos(unitsEgoCA1,comp,dim),egoHbaLoc.control.meanPos(unitsEgoCA1,comp,dim),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
%plot(egoHbaPause.perm.ca1.zscore(:,3,2),egoHbaLoc.perm.ca1.zscore(:,3,2),'.');
%circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Left Vs Right'});
xlabel([dimlabel,' displacement cm (pause)']);
ylabel([dimlabel,' displacement cm (loc)']);
daspect([1,1,1]);
saveas(gcf,fullfile(figpath,['hba_permutation_',dimlabel,'_LvR.pdf']),'pdf')


figure,
complabels = {'left','center','right'};
for comp = 1:3
subplot(1,3,comp)
    plot(egoHbaLoc.control.meanPos(unitsEgoCA1,comp,2),egoHbaLoc.control.meanPos(unitsEgoCA1,comp,1),'.');
    daspect(gca(),[1,1,1]);grid(gca(),'on');
    title(complabels{comp});
end
dimlabel = 'fwdVlat_loc';
saveas(gcf,fullfile(figpath,['hba_permutation_',dimlabel,'_LvR.pdf']),'pdf')



figure();
subplot(121);
hold('on');
plot(egoHbaPause.perm.ca1.zscore(:,3,2),egoHbaLoc.perm.ca1.zscore(:,3,2),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Left Vs Right'});
xlabel('z-score (pause)');
ylabel('z-score (loc)');
daspect([1,1,1]);

subplot(122);
hold('on');
plot(egoHbaPause.control.meanPos(unitsEgoCA1,3,2),egoHbaLoc.control.meanPos(unitsEgoCA1,3,2),'.');daspect(gca(),[1,1,1]);grid(gca(),'on');
%plot(egoHbaPause.perm.ca1.zscore(:,3,2),egoHbaLoc.perm.ca1.zscore(:,3,2),'.');
%circle(0,0,egoHba.perm.ca1.sig,'r-');
title({'Permutation Test','HBA - Left Vs Right'});
xlabel('lateral displacement (pause)');
ylabel('lateral displacement (loc)');
daspect([1,1,1]);


% LOC 
figure,
for comp = 1:3
    for u =1:10
        subplot2(10,3,u,comp);
            hold(gca(),'on');
            set(pcolor(egoHbaRmaps_loc.xpos,...
                       egoHbaRmaps_loc.ypos,...
                       fliplr(rot90(egoHbaRmaps_loc.rmap{20}(:,:,u,comp)',-1))...
                       ),...
                'edgecolor','none');
            xlim(gca(),[-200,200])
            ylim(gca(),[-150,250])
            plot(egoHbaLoc.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2)*10,...
                 egoHbaLoc.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1)*10,...
                 '*m');
            plot(egoHbaLoc.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2),...
                 egoHbaLoc.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1),...
                 '*g');
end
end

% PAUSE 
figure,
for comp = 1:3
    for u =1:10
        subplot2(10,3,u,comp);
            hold(gca(),'on');
            set(pcolor(egoHbaRmaps_pause.xpos,...
                       egoHbaRmaps_pause.ypos,...
                       fliplr(rot90(egoHbaRmaps_pause.rmap{20}(:,:,u,comp)',-1))...
                       ),...
                'edgecolor','none');
            xlim(gca(),[-200,200])
            ylim(gca(),[-150,250])
            plot(egoHbaPause.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2)*10,...
                 egoHbaPause.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1)*10,...
                 '*m');
            plot(egoHbaPause.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2),...
                 egoHbaPause.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1),...
                 '*g');
end
end




figure,
for comp = 1:3
    for u =1:10
        subplot2(10,3,u,comp);
            hold(gca(),'on');
            set(pcolor(egoHbaRmaps.xpos,...
                       egoHbaRmaps.ypos,...
                       fliplr(rot90(egoHbaRmaps.rmap{20}(:,:,u,comp)',-1))...
                       ),...
                'edgecolor','none');
            xlim(gca(),[-200,200])
            ylim(gca(),[-150,250])
            plot(egoHba.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2)*10,...
                 egoHba.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1)*10,...
                 '*m');
            plot(egoHba.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,2),...
                 egoHba.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),comp,1),...
                 '*g');
end
end


% Ctrl Phz
figure,
for comp = 1:3
    for u =1:10
        subplot2(10,3,u,comp);
            hold(gca(),'on');
            set(pcolor(egoHbaPhzRmaps.xpos,...
                       egoHbaPhzRmaps.ypos,...
                       fliplr(rot90(egoHbaPhzRmaps.rmap{20}(:,:,u,3,comp)',-1))...
                       ),...
                'edgecolor','none');
            xlim(gca(),[-200,200])
            ylim(gca(),[-150,250])
            plot(egoHbaPhz.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),3,comp,2)*10,...
                 egoHbaPhz.control.meanPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),3,comp,1)*10,...
                 '*m');
            plot(egoHbaPhz.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),3,comp,2),...
                 egoHbaPhz.control.peakPos(ismember(egoCluSessionMap,[20,unitsEgo{20}(u)],'rows'),3,comp,1),...
                 '*g');
end
end



figure,
for comp = 1:3
for phzInd = 1:3
subplot2(3,3,4-phzInd,comp);
plot(egoHbaPhzPause.perm.ca1.zscore(:,phzInd,comp,2),egoHbaPhzLoc.perm.ca1.zscore(:,phzInd,comp,2),'.');hold(gca(),'on');circle(0,0,egoHbaPhzLoc.perm.ca1.sig,'r-');
end;end;



figure,
hbaClr = 'bgr'
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        subplot2(3,3,4-phzInd,hbaInd);
        hold('on');
        plot(egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,2), ...
             egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,1),['.',hbaClr(hbaInd)]);
        xlim(gca(),[-15,15]);
        ylim(gca(),[-15,20]);
        grid(gca(),'on');
    end;
end



figure,
for hbaInd = 1:3
    subplot(3,1,hbaInd)
histcirc(atan2(egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,2), ...
               egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,1))-0.24);
end

figure,
for hbaInd = 1:3
    subplot2(3,3,1,hbaInd);
histogram(atan2(egoHbaPhz.control.meanPos(:,phzInd,hbaInd,2), ...
                egoHbaPhz.control.meanPos(:,phzInd,hbaInd,1)),...
          linspace(-pi,pi,36));
    subplot2(3,3,2,hbaInd)
histogram(atan2(egoHbaPhzPause.control.meanPos(:,phzInd,hbaInd,2), ...
                egoHbaPhzPause.control.meanPos(:,phzInd,hbaInd,1)),...
          linspace(-pi,pi,36));
    subplot2(3,3,3,hbaInd)
histogram(atan2(egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,2), ...
                egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,1)),...
          linspace(-pi,pi,36));
end


figure,
for hbaInd = 1:3
    subplot(3,1,hbaInd);
plot(egoHbaPhzPause.control.meanPos(:,phzInd,hbaInd,2), ...
     egoHbaPhzPause.control.meanPos(:,phzInd,hbaInd,1),'.')
        xlim(gca(),[-15,15]);
        ylim(gca(),[-15,20]);
end


hbaInd = 1;
circ_wwtest(atan2(egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,2), ...
                  egoHbaPhzLoc.control.meanPos(:,phzInd,hbaInd,1))-0.24,...
            atan2(egoHbaPhzLoc.control.meanPos(:,phzInd,2,2), ...
               egoHbaPhzLoc.control.meanPos(:,phzInd,2,1))-0.24)


figure,
hbaClr = 'bgr'
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        subplot2(3,3,4-phzInd,hbaInd);
        hold('on');
        plot(egoHbaPhz.control.meanPos(:,phzInd,hbaInd,2), ...
             egoHbaPhz.control.meanPos(:,phzInd,hbaInd,1),['.',hbaClr(hbaInd)]);
        xlim(gca(),[-15,15]);
        ylim(gca(),[-15,20]);
        grid(gca(),'on');
    end;
end

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

figure();
hold('on');
for phzInd = 1:phzBin.count
    subplot2(phzBin.count, 1, phzBin.count+1-phzInd, 1);
    hold('on')
    for hbaInd = 1:hbaBin.count
        cdfplot((1-double(hbaInd==1)*2)*egoHbaPhz.boot.ca1.zscore(:,phzInd,hbaInd,2));
    end
    Lines(egoHbaPhz.boot.ca1.sig,[],'r');
    Lines(-egoHbaPhz.boot.ca1.sig,[],'r');
    xlim([-20,20]);
end

figure,
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    hold(gca(),'on');
    for hbaInd = 1:hbaBin.count
        set([histogram(egoHbaPhz.boot.ca1.zscore(:,phzInd,hbaInd,2), ...
                      linspace(-20,20,16)),'FaceAlpha',0.3);
    end
end


figure,
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    hold(gca(),'on');
    for hbaInd = 1:hbaBin.count
% $$$         set(histogram(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2), ...
% $$$                       linspace(-15,15,32)),'FaceAlpha',0.3);
        [ehpcmpKDE,dxi] = ksdensity(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2));
        plot(dxi,ehpcmpKDE,'-','color',hbaBin.color(hbaInd,:))
        med = median(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd));
    end
end
linkax('xy')
xlim([-10,10])
ylim([0,0.16])




hbaBin.color = 'gbr';
figure,
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    hold(gca(),'on');
    for hbaInd = 1:hbaBin.count
% $$$         set(histogram(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2), ...
% $$$                       linspace(-15,15,32)),'FaceAlpha',0.3);
        [ehpcmpKDE,dxi] = ksdensity(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1));
        plot(dxi,ehpcmpKDE,['-',hbaBin.color(hbaInd)])
        med = median(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1))
        [~,xi] = NearestNeighbour(dxi,med);
        line(dxi(xi)*[1,1],[0,ehpcmpKDE(xi)],'color',hbaBin.color(hbaInd));
    end
end
linkax('xy')
xlim([-15,20])
ylim([0,0.16])


% UNIT HBA lateral 
figure,
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
% $$$     boxplot(reshape(sq(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,:,2)),[],1),...
% $$$             reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
% $$$             'orientation','horizontal');
    boxplot(reshape(sq(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,:,2)),[],1),...
            reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
            'orientation','horizontal');

end
linkax('x')
xlim([-15,15])



% UNIT HBA anteroposterior 
figure,
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    boxplot(reshape(sq(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,:,1)),[],1),...
            reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
            'symbol',     'r.',...
            'plotstyle',  'traditional',...
            'orientation','vertical');    
end
linkax('x')
ylim([-15,25])

figure
for hbaInd = 1:hbaBin.count
    subplot(1,3,hbaInd);
    rose(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,2), ...
               egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,1)-2),64);
end
figure
tids = {3:5,18:25,29};
for tid = 1:3
for hbaInd = 1:hbaBin.count
    subplot2(3,3,tid,hbaInd);
    uinds = unitsEgoCA1(ismember(egoCluSessionMap(unitsEgoCA1,1),tids{tid}));
    rose(atan2(egoHbaPhz.control.meanPos(uinds,3,hbaInd,2), ...
               egoHbaPhz.control.meanPos(uinds,3,hbaInd,1)-2),16);
end
end

figure
tids = {3:5,[6,7,27],18:25,29};
for tid = 1:4
for hbaInd = 1:hbaBin.count
    subplot2(4,3,tid,hbaInd);
    uinds = ismember(egoCluSessionMap(:,1),tids{tid});
    rose(atan2(egoHbaPhz.control.meanPos(uinds,3,hbaInd,2), ...
               egoHbaPhz.control.meanPos(uinds,3,hbaInd,1)-2),16);
end
end



figure,
for hbaInd = 1:hbaBin.count
    subplot(1,3,hbaInd);
plot(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,2), ...
           egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,1)-2),...
     sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,2).^2 ...
                      +(egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,1)-2).^2),...
                 '.');
end

uang = sq(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,2), ...
                     egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,1)-2));
udist =sq(sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,2).^2 ...
                      +(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,1)-2).^2));
uang = reshape(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,2), ...
                     egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,1)-2),[],1);
udist =reshape(sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,2).^2 ...
                      +(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,1)-2).^2),[],1);
figure
hold('on');
plot(uang(:,1),udist(:,1),'.g');
plot(uang(:,2),udist(:,2),'.b');
plot(uang(:,3),udist(:,3),'.r');

figure,
boxplot(reshape(udist,[],1),...
        reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
            'symbol',     'r.',...
            'plotstyle',  'traditional',...
            'orientation','vertical');    
figure,
boxplot(reshape(egoHbaPhz.control.meanPos(unitsEgoCA1,3,:,1)-2,[],1),...
        reshape(ones([numel(unitsEgoCA1),1])*[1,2,3],[],1),...
            'symbol',     'r.',...
            'plotstyle',  'traditional',...
            'orientation','vertical');    




figure,
rose(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,2), ...
           egoHbaPhz.control.meanPos(unitsEgoCA1,3,hbaInd,1)),33);




figure,
hold('on');
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2),egoHbaPhz.boot.ca1.zscore(:,phzInd,3,2),'.');
Lines([],egoHbaPhz.boot.ca1.sig,'r');
grid(gca(),'on')
daspect(gca(),[1,1,1]);


figure,
hold('on');
plot(egoHbaPhz.boot.ca1.zscore(:,phzInd,3,2),-egoHbaPhz.boot.ca1.zscore(:,phzInd,1,2),'.');
circle(0,0,egoHbaPhz.boot.ca1.sig,'r-');
grid(gca(),'on')
daspect(gca(),[1,1,1]);

for phzInd = 1:phzBin.count
    subplot2(phzBin.count, 1, phzBin.count+1-phzInd, 1);
    hold('on')
    for hbaInd = 1:hbaBin.count+1
        if hbaInd == 2
            cdfplot(-egoHbaPhz.perm.ca1.zscore(:,phzInd,hbaInd,2));
        else            
            cdfplot(egoHbaPhz.perm.ca1.zscore(:,phzInd,hbaInd,2));
        end
    end
    Lines(egoHbaPhz.perm.ca1.sig,[],'r');
    Lines(-egoHbaPhz.perm.ca1.sig,[],'r');
    xlim([-20,20]);
end

figure();
hold('on');
for phzInd = 1:phzBin.count
    subplot2(phzBin.count, 1, phzBin.count+1-phzInd, 1);
    hold('on')
    for hbaInd = 1:hbaBin.count
        cdfplot(-egoHbaPhz.perm.ca1.zscore(:,phzInd,hbaInd,1));
    end
    Lines(egoHbaPhz.perm.ca1.sig,[],'r');
    Lines(-egoHbaPhz.perm.ca1.sig,[],'r');
    xlim([-20,20]);
end



figure();
hold('on');
for phzInd = 1:phzBin.count
    subplot2(phzBin.count, 1, phzBin.count+1-phzInd, 1);
    hold('on')
    for hbaInd = 1:hbaBin.count
        if hbaInd == 2
            cdfplot(-egoHbaPhz.perm.ca3.zscore(:,phzInd,hbaInd,2));
        else            
            cdfplot(egoHbaPhz.perm.ca3.zscore(:,phzInd,hbaInd,2));
        end
    end
    Lines(egoHbaPhz.perm.ca3.sig,[],'r');
    Lines(-egoHbaPhz.perm.ca3.sig,[],'r');
    xlim([-20,20]);
end

% 1-(1-0.05)^(1/numel(unitsEgoCA1)) -> 0.000410262174549647 -> zscore > 3.352795
% 1-(1-0.001)^(1/numel(unitsEgoCA1)) -> 8.00397063671632e-06 -> zscore > 4.314457
% 1-(1-0.0001)^(1/numel(unitsEgoCA1)) -> 8.00397063671632e-06 -> zscore > 4.798339

figure();
hold('on');
plot(egoHbaPhz.perm.ca1.zscore(:,phzInd,1,2), -egoHbaPhz.perm.ca1.zscore(:,phzInd,2,2),'.');
circle(0,0,3.35279,'r');
circle(0,0,4.314457,'r');
circle(0,0,4.798339,'r');
grid('on');

sum(sum([egoHbaPhz.perm.ca1.zscore(:,phzInd,1,2), -egoHbaPhz.perm.ca1.zscore(:,phzInd,2,2)]>3.35,2)==2)./numel(unitsEgoCA1)
% 20.8% of CA1 place cell shift thier prospective spatial field
% laterally past the head towards the side of 
sum(sum([egoHbaPhz.perm.ca1.zscore(:,phzInd,1,2), -egoHbaPhz.perm.ca1.zscore(:,phzInd,2,2)]>3.35,2)>=1)./numel(unitsEgoCA1)


% Center head-body-ang ego-field lateral difference to left and right
% head-body-ang egofields
figure();
phzInd = 3;
hold('on');
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2),...
     '.');
grid('on');
xlim([-15,15]);
ylim([-15,15]);


% Center head-body-ang ego-field foward difference to left and right
figure();
phzInd = 3;
hold('on');
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,1),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,1),...
     '.');
grid('on');
xlim([-15,15]);
ylim([-15,15]);


figure,
phzInd = 3;
% Right
hold('on');                                 
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,1,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,1,2,1),...
     '.r');
% Left
hold('on');                                 
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,1,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,1,2,1),...
     '.g');
grid('on');
xlim([-15,15]);
ylim([-15,20]);

figure,
subplot(131);
rose(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)),...
     32);
subplot(132);
rose(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)),...
     32);
subplot(133);
rose(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)),...
     32);

circ_mean(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)))

circ_mean(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)))

circ_mean(atan2(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,2),...
     egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,1)-egoHbaPhz.control.meanPos(unitsEgoCA1,2,2,1)))

figure();
hold('on');
% $$$ plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2), ...
% $$$      egoHbaPhz.perm.ca1.zscore(:,phzInd,3,2),'.');
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,1,3,2), ...
     -egoHbaPhz.perm.ca1.zscore(:,phzInd,2,2),'.r');
plot(egoHbaPhz.control.meanPos(unitsEgoCA1,1,1,2), ...
     egoHbaPhz.perm.ca1.zscore(:,phzInd,1,2),'.g');

%[RHO,PVAL]=corr(egoHbaPhz.control.meanPos(unitsEgoCA1,1,2,2),egoHbaPhz.perm.zscoreC)
[RHO,PVAL]=corr(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2), ...
     egoHbaPhz.perm.ca1.zscore(:,phzInd,1,2))
[RHO,PVAL]=corr(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2), ...
     egoHbaPhz.perm.ca1.zscore(:,phzInd,2,2))
xlim([-30,30]);

figure
dat = load('sunspot.dat');
rose(dat(:,2));
texthandles = findall(gca,'Type','Text');
set(texthandles,'Rotation',-90)
view([0,0]);



figure,
for a = 1:3
    subplot(3,1,a);
    hold('on')
        xlim([-15,15]);
    ylim([-15,15]);
    Lines(0,[],'k');
    plot(egoHba.control.meanPos(unitsEgoCA1,a,2), ...
         egoHba.control.meanPos(unitsEgoCA1,a,1),'.')
end

figure
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2), ...
         egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,3,2),...
         '.');
    daspect([1,1,1]);
    grid(gca(),'on');
end
linkax('xy');    
xlim([-15,15]);
ylim([-20,20]);


figure
for phzInd = 1:phzBin.count
    subplot(3,1,phzBin.count+1-phzInd);
    plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,1,2)-egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2), ...
         egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,2,2),...
         '.');
    daspect([1,1,1]);
    grid(gca(),'on');
end
linkax('xy');    
xlim([-15,15]);
ylim([-20,20]);


figure();
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        subplot2(3,3,phzBin.count+1-phzInd,hbaInd);    
    plot(egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2),...
         egoHbaPhz.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1), '.');
    daspect(gca(),[1,1,1]);
    grid(gca(),'on');
    end
end
linkax('xy');
xlim([-15,15]);
ylim([-20,20]);


figure();
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        subplot2(3,3,phzBin.count+1-phzInd,hbaInd);    
    plot(egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2),...
         egoHbaPhzLoc.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1), '.');
    daspect(gca(),[1,1,1]);
    grid(gca(),'on');
    end
end
linkax('xy');
xlim([-15,15]);
ylim([-20,20]);


figure();
for phzInd = 1:phzBin.count
    for hbaInd = 1:hbaBin.count
        subplot2(3,3,phzBin.count+1-phzInd,hbaInd);    
    plot(egoHbaPhzPause.control.meanPos(unitsEgoCA1,phzInd,hbaInd,2),...
         egoHbaPhzPause.control.meanPos(unitsEgoCA1,phzInd,hbaInd,1), '.');
    daspect(gca(),[1,1,1]);
    grid(gca(),'on');
    end
end
linkax('xy');
xlim([-15,15]);
ylim([-20,20]);




figure,
for a = 1:3
    hold('on')
    xlim([-15,15]);
    Lines(0,[],'k');
    cdfplot(egoHba.control.meanPos(unitsEgoCA1,a,2));
end

figure
hold('on');
cdfplot(egoHba.control.meanPos(unitsEgoCA1,1,2)-egoHba.control.meanPos(unitsEgoCA1,2,2));
cdfplot(egoHba.control.meanPos(unitsEgoCA1,3,2)-egoHba.control.meanPos(unitsEgoCA1,2,2));
xlim([-15,15]);

figure,
for a = 1:3
    hold('on')
    xlim([-15,15]);
    Lines(0,[],'k');
    cdfplot(egoHba.control.meanPos(unitsEgoCA1,a,1));
end


figure,
for a = 1:3
    subplot(3,1,a)
    if a==1,Lines(0,[],'k');end    
    plot(egoHba.control.meanPos(unitsEgoCA3,a,2), ...
         egoHba.control.meanPos(unitsEgoCA3,a,1),'.')
    xlim([-15,15]);
    ylim([-15,15]);
end


hdist = sqrt(egoHbaPhz.control.meanPos(unitsEgoCA1, phzInd, :, lat).^2 ...
                +(egoHbaPhz.control.meanPos(unitsEgoCA1, phzInd, :, fwd)-2).^2);


[h,p] = ttest(hdist(:,1),hdist(:,2))
[h,p] = ttest(hdist(:,1),hdist(:,3))

[h,p] = ttest(hdist(:,2),hdist(:,3))

1-(1-0.05)^(1/numel(unitsEgoCA1))
