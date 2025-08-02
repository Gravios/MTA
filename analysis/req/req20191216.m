
% REQUIRES
% req20191104


hpchBinEdges = linspace(-1.4,0,16);
hpchBinCenters = mean([hpchBinEdges(2:end); hpchBinEdges(1:end-1)]);
edy = numel(hpchBinCenters);
hpchBinInd = discretize(dfet,hpchBinEdges);


eds = 4;
hbangBinEdges = linspace(-1.2,1.2,eds);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);


hvelBinEdges = linspace(0,1.5,10);
hvelBinCenters = mean([hvelBinEdges(2:end); hvelBinEdges(1:end-1)]);
hvelBinInd = discretize(log10(dxyvel),hvelBinEdges);



numIter = 100;
qtls = [0.25,0.50,0.75];

%%%<<< Ego centric phase precession bootstrap (hbang,hvang) 
qntlavHBP = nan(  [numel(qtls), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
                2,                        ...
                numIter,                  ...
                2                         ...               
             ]                            ...
           );

% pxav( x, y, hba, hav, sts, fet, quad, itr)
shifts = 0:8:2^8;
for sts = 1:numel(stsGrps),
    ind =   logical(dstcm(:,1))                             ...
            & any(dstcm(:, stsGrps{sts}),2)                  ...
            & dpostI                                        ...
            & duincI                                        ...
            & dhdist<380;
    ind = ind & ~ismember(dtind,[3,4,5]);
    ind = ind & dxyvel>1 & dxyvel<32;
    disp(['[Info] State: ',stsLbls{sts}]);
    
    tdphz = dphz(ind);
    thpchBinInd = hpchBinInd(ind);
    thbangBinInd = hbangBinInd(ind);
    tsmMask = smMask(ind);
    tvel = hvelBinInd(ind);

    for e = 1:2
        tferr = ferr{e}(ind);
        disp(['[Info] e: ',num2str(e)]);
        for f = 1:edy,
            for b = 1:edx,
                tic
                for i = 1:numIter,
                    indb = thpchBinInd == f & thbangBinInd == b & circshift(tsmMask, randsample(shifts,1));
                    %indb = thpchBinInd == f & ismember(thbangBinInd,[3,4,5]) & circshift(tsmMask, randsample(shifts,1));
                    indv = [];
                    for v = 1:numel(hvelBinCenters),
                        if sum(tvel(indb)==v)<400
                            indv = [indv;find(tvel(indb)==v)];
                        else
                            indv = [indv;randsample(find(tvel(indb)==v),400,false)];
                        end
                    end
                    if numel(indv) < 1000, 
                        continue; 
                    else
                        indv = randsample(indv,round(numel(indv)./2));
                    end;
                    vferr = tferr(indb);
                    vferr = vferr(indv);
                    vdphz = tdphz(indb);
                    vdphz = vdphz(indv);
                    
                    
                    pxavnHBP =               ...
                        histcounts2(vferr,          ...
                                    vdphz,             ...
                                    ferrorBinEdges{e},       ...
                                    phzBins);
                    pxavnHBP = cumsum(bsxfun(@rdivide,pxavnHBP,sum(pxavnHBP)));
                    for p = 1:numel(phzBinCenters),
                        try,
                            [x, index] = unique(pxavnHBP(:,p),'first');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavHBP(:,p,b,f,sts,e,i,1) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                        try,
                            [x, index] = unique(pxavnHBP(:,p),'last');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavHBP(:,p,b,f,sts,e,i,2) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                    end
                end
                toc
            end
        end
    end
end
temp_qntlavHBP = qntlavHBP;
temp_qntlavHBP2 = qntlavHBP;
qntlavHBP = mean(qntlavHBP,8,'omitnan');


distHBPa = sqrt(sq(sum(diff(qntlavHBP(2,[7,3],:,:,1,:,:),1,2).^2,6)));

distHBP = sqrt(sq(sum(diff(mean(qntlavHBP(2,[7,3],:,:,1,:,:),7,'omitnan'),1,2).^2,6,'omitnan')));
lonHBP =  sq(diff(median(qntlavHBP(2,[7,3],:,:,1,1,:),7,'omitnan')));
latHBP =  sq(diff(median(qntlavHBP(2,[7,3],:,:,1,2,:),7,'omitnan')));

figure
subplot(291);
imagesc(round(hbangBinCenters,1),hpchBinCenters,lonHBP');
axis('xy');
xlabel('Yaw (rad)')
ylabel('Pitch (rad)');
title('Lon');
subplot(292);
imagesc(round(hbangBinCenters,1),hpchBinCenters,latHBP');
xlabel('Yaw (rad)')
caxis([-70,70]);
axis('xy');
title('Lat');
subplot(293);
imagesc(round(hbangBinCenters,1),hpchBinCenters,distHBP');
title('Range');
xlabel('Yaw (rad)')
caxis([-20,120]);
colormap('jet');
axis('xy');
subplot(2,9,[4:9]);
boxplot(sq(distHBPa(2,1:end,:))','labels',round(hpchBinCenters,2),'labelorientation','inline')
title({'Bootstrapped Median Phase Precession Distance',...
       'Vertical Head Pitch VS PPDist '});
xlabel('Head Pitch (rad)')
ylabel('PP Range (mm)');
subplot(2,9,[13:18]);
plot(hpchBinCenters,-acos(distHBP(2,1:end)./105),'.');
line([-1.4,-0.2],[-1.4,-0.2]+0.1)
legend({'data','y=x+0.1'},'Location','southeast');
xlabel('Head Pitch (rad)')
ylabel('PP Pitch (rad)');
axes('Position',[0,0,1,1],'Visible','off')
text(0.1,0.4,{'Subject: jg05','','Each bin resampled to have ','uniform head speed distribution','between 1-30 cm/s'});




scl = 100;
figure
subplot(121);
hold('on');
cmap = cool(3);
for p = 1:3,
    plot(cos(hpchBinCenters(1:end)).*scl,...
         distHBP(p,1:end),...
         '.',...
         'MarkerSize',10,...
         'Color', ...
         cmap(p,:));
end
plot(cos(hpchBinCenters(2:end-1)).*scl,median(distHBP(2:end-1,2:end-1)),'-+');
line([40,120],[40,120]);
subplot(122);
hold('on');
cmap = cool(11);
for p = 2:12,
    plot(cos(hpchBinCenters(2:end-1)).*scl,...
         distHBP(2:end-1,:),...
         '.',...
         'MarkerSize',10,...
         'Color', ...
         cmap(p-1,:));
end
plot(cos(hpchBinCenters(2:end-1)).*scl,median(distHBP(2:end-1,2:end-1),2),'-+');

figure();
subplot(121);
plot(hpchBinCenters,(median(distHBP,1)),'-+');
plot(hpchBinCenters,-acos(median(distHBP,1)./),'-+');
ylim([-1.5,0]);xlim([-1.5,0]);

subplot(122);
plot(hpchBinCenters,median(distHBP,2),'-+');

figure,    plot(hpchBinCenters(1:end),distHBP(7,:),'.')
hold('on');
line([40,120],[40,120]);