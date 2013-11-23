function rear_stats(Session,ccg_name)


Session = MTASession('jg05-20120310',{'CluRes'});

numClu = size(Session.map,1);

ccg_name = 'rdpf_nrhp_sur_par1';

load([Session.spath.analysis Session.filebase '.ccg.' ccg_name '.mat'])

fccg = Bccg.filter(gausswin(9));


%% Compute main CCG statistics
fconf = 0.95;

upbd = sq(max(fccg,[],1));
supbd = sort(upbd,3);
upperConfidenceBoundary = sq(supbd(:,3,round(size(fccg,5)*fconf)));

lwbd = sq(min(fccg,[],1));
slwbd = sort(lwbd,3,'descend');
lowerConfidenceBoundary = sq(slwbd(:,3,round(size(fccg,5)*fconf)));


unit = 92;
clf
bar(Bccg.tbin,fccg(:,unit,2,1,1))
Lines([],upperConfidenceBoundary(unit),'r')
Lines([],lowerConfidenceBoundary(unit),'r')
hold on
plot(Bccg.tbin,mean(fccg(:,unit,3,1,:),5),'g')




% $$$ srf = sq(sort(fccg(:,:,1,:),1));
% $$$ p95 = sq(srf(round(size(fccg,1)*fconf),:,:));
% $$$ dsrf = sq(sort(fccg(:,:,1,:),1,'descend'));
% $$$ p05 = sq(dsrf(round(size(fccg,1)*fconf),:,:));
% $$$ 
%p_ron  = sq(sum(fccg(:,:,1,:,:)>fccg(:,:,3,:,:),1)./size(fccg,5));
%p_roff = sq(sum(fccg(:,:,2,:,:)>fccg(:,:,3,:,:),1)./size(fccg,5));


%% Compute spiking significance bins

RearSpkSig = zeros(length(Bccg.tbin),numClu,2,2);
for b = 1:2,
    for unit = 1:numClu,
        for t = 2:length(Bccg.tbin),
            RearSpkSig(t,unit,b,1) = RearSpkSig(t-1,unit,b,1) + sum(fccg(t,unit,b,1,1)<lowerConfidenceBoundary(unit));
            RearSpkSig(t,unit,b,2) = RearSpkSig(t-1,unit,b,2) + sum(fccg(t,unit,b,1,1)>upperConfidenceBoundary(unit));
        end
    end
end




%% Compute max/min spiking sig time lag
rssi = cat(1,zeros(1,numClu,2,2),diff(RearSpkSig));
rssi(rssi>0) = 1;
mmstl = zeros(numClu,2,2,2);
for b = 1:2,
    for unit = 1:numClu,
        if ~isempty(min(fccg(rssi(:,unit,b,1)==1,unit,b,1,1),[],1)),
            [mmstl(unit,1,b,1),mmstl(unit,2,b,1)] = min(fccg(rssi(:,unit,b,1)==1,unit,b,1,1),[],1);
        end
        if ~isempty(max(fccg(rssi(:,unit,b,2)==1,unit,b,1,1),[],1)),
            [mmstl(unit,1,b,2),mmstl(unit,2,b,2)] = max(fccg(rssi(:,unit,b,2)==1,unit,b,1,1),[],1);
        end
    end
end


[mmstl(:,1,:,1),mmstl(:,2,:,1)] = min(fccg(rssi(:,:,:,1)==1,:,:,1,1));
[mmstl(:,1,:,2),mmstl(:,2,:,2)] = max(fccg(rssi(:,:,:,2)==1,:,:,1,1));




%% Load Place Fields
Session.Pfs = [];
Session = Session.load_Pfs;

pf_search.mazeName = 'cof';
pf_search.trialName = Session.trialName;
pf_search.trackingMarker = Session.trackingMarker;
pf_search.stateLabel = 'head';
pf_search.spk_shuffle = 'n';
pf_search.pos_shuffle = 0;
pf_search.numBSiterations = 1;
pf_search.numZslices = 1;
pf_search.nbins = 50;
pf_search.smooth = 0.03;




%% Display Everything
figure(102)
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
    clf
    %% Rate Map
    subplot2(6,2,[1,2],1);
    pf_search.stateLabel = 'head';
    Pfs = Session.getPfs(pf_search);
    ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})
    hold on
    scatter(wpos(gwi,1),wpos(gwi,2),20,'m','Marker','*')
    title(num2str(unit))
    %
    subplot2(6,2,[1,2],2);
    pf_search.stateLabel = 'rear';
    Pfs = Session.getPfs(pf_search);
    ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})
    hold on
    scatter(rpos(gri,1),rpos(gri,2),20,'m','Marker','*')

    %% CCGs
    subplot2(6,2,3,1);
    bar(tbin/1000,wfccg(:,2,unit));axis tight; grid on
    %
    subplot2(6,2,4,1);
    bar(tbin/1000,wfccg(:,1,unit));axis tight; grid on
    %
    subplot2(6,2,3,2);
    bar(tbin/1000,rfccg(:,1,unit));axis tight; grid on
    hold on
    Lines([],[upperConfidenceBoundary(unit),lowerConfidenceBoundary(unit)],'r')
    Lines([],mean(mean(fccg(:,:,1,unit),1),2),'g')
    if max(sq(rfccg(:,1,unit)))<upperConfidenceBoundary(unit), ylim([0,1+upperConfidenceBoundary(unit)]);,end
    plot(tbin,p95(:,unit),'m')
    plot(tbin,p05(:,unit),'m')
    subplot2(6,2,4,2);
    bar(tbin/1000,rfccg(:,2,unit));axis tight; grid on
    hold on
    Lines([],[upperConfidenceBoundary(unit),lowerConfidenceBoundary(unit)],'r')
    Lines([],mean(mean(bfccg(:,:,1,unit),1),2),'g')
    if max(sq(rfccg(:,2,unit)))<upperConfidenceBoundary(unit), ylim([0,1+upperConfidenceBoundary(unit)]);,end
    plot(tbin,p95(:,unit),'m')
    plot(tbin,p05(:,unit),'m')

    %% Rasters
    subplot2(6,2,5,2);
    plot(tron(conr==unit)/1000,tonr(conr==unit),'.','markersize',5);axis tight; grid on
    subplot2(6,2,6,2);
    plot(troff(coffr==unit)/1000,toffr(coffr==unit),'.','markersize',5);axis tight; grid on
    subplot2(6,2,5,1);
    plot(twoff(coffw==unit)/1000,toffw(coffw==unit),'.','markersize',5);axis tight; grid on
    subplot2(6,2,6,1);
    plot(twon(conw==unit)/1000,tonw(conw==unit),'.','markersize',5);axis tight; grid on

    %% ReportFig
    %reportfig(102,[Session.filebase '.rearing'],[],['unit: ' num2str(unit)],[],0);
    
    %% Figure controls
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #: ');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
    end
end
