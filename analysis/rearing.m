function rearing(Session)
%,trialName)
% $$$ 
% $$$ trialName = 'crt1';
% $$$ 
% $$$ if ~isa(Session,'MTASession'),
% $$$     Session = MTASession(Session);
% $$$ else
% $$$     Session = MTATrial(Session,trialName);
% $$$ end

% $$$    wind = 64;
% $$$    nbins = round(sum(diff(stsp,1,2)+1)/wind);
% $$$    
% $$$    winres = round(myRes/wind);
% $$$    subs = [winres(winres~=0),myClu(winres~=0)];
% $$$    vals = ones(size(subs,1),1);
% $$$    spkrate = accumarray(subs,vals,[nbins,size(Map,1)],@(x) sum(x)/(wind/Session.lfpSampleRate));

%Session = MTASession('jg05-20120309');
Session = MTASession('jg05-20120310',{'CluRes','ang'});
%Session = MTASession('jg04-20120128');
%Session = MTASession('jg04-20120130');
%Session = MTASession(Session);
%Session = MTATrial(Session,'all');



%% Calculate rearing and non-rearing periods
rper = Session.Bhv.getState('rear').state;
% $$$ wper = Session.Bhv.getState('walk').state;
% $$$ hper = Session.Bhv.getState('head').state;

[~,rfet] = rear(Session,'com');
% $$$ [~,wfet] = walk(Session,'com');
% $$$ [~,hfet] = head(Session,'com');

rpos = sq(Session.xyz(rper(:,1),Session.Model.gmi(Session.trackingMarker),[1,2]));
% $$$ wpos = sq(Session.xyz(wper(:,1),Session.Model.gmi(Session.trackingMarker),[1,2]));

rdur = diff(rper,1,2);
% $$$ wdur = diff(wper,1,2);
% $$$ hdur = diff(hper,1,2);

% $$$ 
% $$$ %% Select Walking onset and offset without previous/future walk/rear
% $$$ rbcom = zeros(size(Session.xyz,1),2,3);
% $$$ rb = Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
% $$$ rbcom(:,1,:) = sq(Session.com(rb));
% $$$ rb = Session.Model.rb({'head_back','head_left','head_front','head_right'});
% $$$ rbcom(:,2,:) = sq(Session.com(rb));
% $$$ 
% $$$ %% Selected Walking inds
% $$$ fet = GetSegs(wfet,round(wper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ fet = reshape(fet,round(6*Session.xyzSampleRate),size(wper,1),2);
% $$$ hrfet = rbcom(:,2,3).*Session.ang(:,Session.Model.gmi('spine_middle'),Session.Model.gmi('spine_upper'),2);
% $$$ hrfet = GetSegs(hrfet,round(wper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ hrfet = reshape(hrfet,round(6*Session.xyzSampleRate),size(wper,1),2);
% $$$ gwi = find(max(fet(1:300,:,1))<1&wdur'>1.5&max(hrfet(:,:,1))<0);
% $$$ gwiu = find(max(fet(450:end,:,2))<1&wdur'>1.5&max(hrfet(:,:,1))<0);
% $$$ 
% $$$ %% Selected Head inds
% $$$ fet = GetSegs(hfet,round(hper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ fet = reshape(fet,round(6*Session.xyzSampleRate),size(hper,1),2);
% $$$ hrfet = rbcom(:,2,3).*Session.ang(:,Session.Model.gmi('spine_middle'),Session.Model.gmi('spine_upper'),2);
% $$$ hrfet = GetSegs(hrfet,round(hper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ hrfet = reshape(hrfet,round(6*Session.xyzSampleRate),size(hper,1),2);
% $$$ ghi = find(max(fet(100:300,:,1))<1&hdur'>1.5&max(hrfet(:,:,1))<0);
% $$$ ghiu = find(max(fet(450:650,:,2))<1&hdur'>1.5&max(hrfet(:,:,1))<0);
% $$$ 
% $$$ %% Selected Rear inds
fet = GetSegs(Session.ang(:,3,4,2),round(rper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
fet = reshape(fet,round(6*Session.xyzSampleRate),size(rper,1),2);
gri = find(max(fet(1:300,:,1))<0&rdur'>1.5);
griu = find(max(fet(450:end,:,2))<0&rdur'>1.5);



%% Calculate walk CCG and TrigRaster

binSize = 100;%ms
halfBins = 50;
normalization = 'hz';

% $$$ per_on  = wper(gwi,1);
% $$$ per_off = wper(gwiu,2);
% $$$ poni  = round((per_on-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
% $$$ poffi = round((per_off-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
% $$$ clear('fccg')
% $$$ [fccg tbin ] = Trains2CCG({poni,poffi,Res},{1,2,Clu},binSize,halfBins,Session.lfpSampleRate,normalization);
% $$$ uClu = unique(Clu);
% $$$ sccg = zeros(size(fccg,1),2,size(Map,1));
% $$$ sccg(:,:,uClu) = sq(fccg(:,1:2,3:end));
% $$$ ffccg = zeros(size(fccg,1),2,size(Map,1));
% $$$ for k=1:size(Map,1)
% $$$     for l=1:2
% $$$         ffccg(:,l,k) = Filter0(gausswin(3),sccg(:,l,k));
% $$$     end
% $$$ end
% $$$ wfccg = ffccg;
% $$$ [twon,tonw,conw] =  TrigRasters(poni,10000,Res,Clu,Session.lfpSampleRate,1,[]);
% $$$ [twoff,toffw,coffw] = TrigRasters(poffi,10000,Res,Clu,Session.lfpSampleRate,1,[]);

%% Calculate rear CCG and TrigRaster
per_on  = rper(gri,1);
per_off = rper(griu,2);
poni  = round((per_on-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
poffi = round((per_off-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
clear('fccg')
clear('ffccg')
[fccg tbin ] = Trains2CCG({poni,poffi,Session.res},{1,2,Session.clu},binSize,halfBins,Session.lfpSampleRate,normalization);
uClu = unique(Session.clu);
sccg = zeros(size(fccg,1),2,size(Session.map,1));
sccg(:,:,uClu) = sq(fccg(:,1:2,3:end));
ffccg = zeros(size(fccg,1),2,size(Session.map,1));
for k=1:size(Session.map,1)
    for l=1:2
        ffccg(:,l,k) = Filter0(gausswin(3),sccg(:,l,k));
    end
end

rfccg = ffccg;

% $$$ [tron,tonr,conr] =  TrigRasters(poni,10000,Res,Clu,Session.lfpSampleRate,1,[]);
% $$$ [troff,toffr,coffr] = TrigRasters(poffi,10000,Res,Clu,Session.lfpSampleRate,1,[]);


%% Calculate Non rearing movment periods
trimLength = 120*5;%samples
temp_non_rearing_periods = [rper(:,1)-trimLength,rper(:,2)+trimLength];
tnrp = [[1;temp_non_rearing_periods(:,2)],[temp_non_rearing_periods(:,1);size(Session.xyz,1)]];
nrp = tnrp(diff(tnrp,1,2)>0,:);
nrwp = IntersectRanges(nrp,hper);
fhp = [];
for i = 1:length(nrwp)
    fhp = cat(1,fhp,[nrwp(i,1):nrwp(i,2)]');
end


% $$$ %Start-bootstrap-ccg
% $$$ numIterations = 1;
% $$$ clear('fccg')
% $$$ bfccg = [];
% $$$ for i=1:numIterations,
% $$$     %sprintf('Iter #: %d',i)
% $$$     randInd = zeros(size(rper,1),1);
% $$$     for j = 1:size(rper,1)
% $$$         randInd(j)=round((fhp(randi(size(fhp,1)))-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
% $$$     end
% $$$     %posind(i,:) = randInd;
% $$$     [fccg tbin ] = Trains2CCG({randInd,Res},{1,Clu},binSize,halfBins,Session.lfpSampleRate,normalization);
% $$$     uClu = unique(Clu);
% $$$     sgccg = zeros(size(fccg,1),1,size(Map,1));
% $$$     sgccg(:,:,uClu) = sq(fccg(:,1,2:end));
% $$$     for k=1:size(Map,1)
% $$$         bfccg(i,:,1,k) = Filter0(gausswin(3)/sum(gausswin(3)),sgccg(:,1,k));
% $$$     end
% $$$ end

%End-bootstrap-ccg


%% Compute main CCG statistics
% $$$ fconf = 0.95;
% $$$ 
% $$$ upbd = sq(max(bfccg,[],2));
% $$$ supbd = sort(upbd,1);
% $$$ upperConfidenceBoundary = sq(supbd(round(size(bfccg,1)*fconf),:,:));
% $$$ 
% $$$ lwbd = sq(min(bfccg,[],2));
% $$$ slwbd = sort(lwbd,1,'descend');
% $$$ lowerConfidenceBoundary = sq(slwbd(round(size(bfccg,1)*fconf),:,:));
% $$$ 
% $$$ srf = sq(sort(bfccg(:,:,1,:),1));
% $$$ p95 = sq(srf(round(size(bfccg,1)*fconf),:,:));
% $$$ dsrf = sq(sort(bfccg(:,:,1,:),1,'descend'));
% $$$ p05 = sq(dsrf(round(size(bfccg,1)*fconf),:,:));

% $$$ p_ron  = sq(sum(repmat(rfccg(:,1,:),size(bfccg,1),1)>bfccg(:,:,1,:),1)./size(bfccg,1));
% $$$ p_roff = sq(sum(repmat(rfccg(:,1,:),size(bfccg,1),1)>bfccg(:,:,2,:),1)./size(bfccg,1));

% $$$ 
% $$$ %% Compute spiking significance bins
% $$$ 
% $$$ RearSpkSig = zeros(size(Map,1),length(tbin),2,2);
% $$$ for b = 1:2,
% $$$     for unit = 1:size(Map,1),
% $$$         for t = 2:length(tbin),
% $$$             RearSpkSig(unit,t,b,1) = RearSpkSig(unit,t-1,b,1) + sum(rfccg(t,b,unit)<lowerConfidenceBoundary(unit));
% $$$             RearSpkSig(unit,t,b,2) = RearSpkSig(unit,t-1,b,2) + sum(rfccg(t,b,unit)>upperConfidenceBoundary(unit));
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ %% Compute max/min spiking sig time lag
% $$$ rssi = cat(2,zeros(size(Map,1),1,2,2),diff(RearSpkSig,[],2));
% $$$ rssi(rssi>0) = 1;
% $$$ mmstl = zeros(size(Map,1),2,2,2);
% $$$ for b = 1:2,
% $$$     for unit = 1:size(Map,1),
% $$$         if ~isempty(min(rfccg(rssi(unit,:,b,1)==1,b,unit))),
% $$$             [mmstl(unit,1,b,1),mmstl(unit,2,b,1)] = min(rfccg(rssi(unit,:,b,1)==1,b,unit));
% $$$         end
% $$$         if ~isempty(max(rfccg(rssi(unit,:,b,2)==1,b,unit))),
% $$$             [mmstl(unit,1,b,2),mmstl(unit,2,b,2)] = max(rfccg(rssi(unit,:,b,2)==1,b,unit));
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ [mmstl(:,1,:,1),mmstl(:,2,:,1)] = min(rfccg(rssi(:,:,:,1)==1,:,:));
% $$$ [mmstl(:,1,:,2),mmstl(:,2,:,2)] = max(rfccg(rssi(:,:,:,2)==1,:,:));



%figure,plot(tbin,sq(RearSpkSig(4,:,2,:)))
% $$$ 
% $$$ 
% $$$ 
% $$$ %% Compute Unit CCGs
% $$$ if strcmp(Session.trialName,'all'),
% $$$ [accg atbin] = autoccg(Session);
% $$$ else    
% $$$ fs = MTASession(Session.name,Session.Maze.name);
% $$$ [accg atbin] = autoccg(fs);
% $$$ clear('fs')
% $$$ end




%% Compute Spectrograms around rear & walk - onset/offset
% $$$ %spec = load([Session.spath.analysis Session.name '.spectrums.mat']);
% $$$ channels = [65:96];
% $$$ lfp = Session.loadlfp(channels);
% $$$ wlfp = WhitenSignal(lfp,[],1);
% $$$ 
% $$$ frange = [1,30];
% $$$ nffts = 2^12;
% $$$ windows =2^11;
% $$$ yr = [];
% $$$ yw = [];
% $$$ t = [];
% $$$ f = [];
% $$$ for b = 1:2,
% $$$     rwlfp = GetSegs(wlfp,...
% $$$                     round(((rper(:,b)./Session.xyzSampleRate)-3).*Session.lfpSampleRate),...
% $$$                     8*Session.lfpSampleRate,0);
% $$$     for j = 1:32;
% $$$         for i = 1:size(rwlfp,2),
% $$$             [yr(:,:,i,b,j),f,t] = mtchglong(rwlfp(:,i,j),...
% $$$                                             nffts,...
% $$$                                             Session.lfpSampleRate,...
% $$$                                             windows,...
% $$$                                             windows*0.875,[],[],[],frange);
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ for b = 1:2,
% $$$     wwlfp = GetSegs(wlfp,...
% $$$                     round(((wper(:,b)./Session.xyzSampleRate)-3).*Session.lfpSampleRate),...
% $$$                     8*Session.lfpSampleRate,0);
% $$$     for j = 1:32;
% $$$         for i = 1:size(wwlfp,2),
% $$$             [yw(:,:,i,b,j),f,t] = mtchglong(wwlfp(:,i,j),...
% $$$                                             nffts,...
% $$$                                             Session.lfpSampleRate,...
% $$$                                             windows,...
% $$$                                             windows*0.875,[],[],[],frange);
% $$$         end        
% $$$     end
% $$$ end
% $$$ t = t+diff(t(1:2))/2;
% $$$ 
% $$$ thpow_rear = sq(median(yr(:,f<12&f>6,:,:,:),2));
% $$$ thpow_walk = sq(median(yw(:,f<12&f>6,:,:,:),2));
% $$$ 
% $$$ ThetaWalkVsRearSig = zeros(length(t),length(channels));
% $$$ for tind = 1:length(t),
% $$$     for chan = 1:length(channels),
% $$$         ThetaWalkVsRearSig(tind,chan) = ranksum(sq(thpow_rear(tind,gri,1,chan)),sq(thpow_walk(tind,gwi,1,chan)));
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ frbcom_gkern = gausswin(5)/sum(gausswin(5));
% $$$ frbcom = zeros(size(rbcom));
% $$$ for mar = 1:size(rbcom,2),
% $$$     for dim = 1:size(rbcom,3),
% $$$         frbcom(:,mar,dim) = Filter0(frbcom_gkern,rbcom(:,mar,dim));
% $$$     end
% $$$ end
% $$$ 
% $$$ rbc_off_segs = GetSegs(frbcom,round(rper(griu,2)-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ rbc_off_segs = reshape(rbc_off_segs,round(6*Session.xyzSampleRate),size(rper(griu,2),1),size(frbcom,2),3);
% $$$ rbcv_off_segs = abs(diff(sq(sqrt(sum(rbc_off_segs.^2,4))),1));
% $$$ rbc_on_segs = GetSegs(frbcom,round(rper(gri,1)-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ rbc_on_segs = reshape(rbc_on_segs,round(6*Session.xyzSampleRate),size(rper(gri,1),1),size(frbcom,2),3);
% $$$ rbcv_on_segs = abs(diff(sq(sqrt(sum(rbc_on_segs.^2,4))),1));
% $$$ 
% $$$ wbc_off_segs = GetSegs(frbcom,round(wper(gwiu,2)-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ wbc_off_segs = reshape(wbc_off_segs,round(6*Session.xyzSampleRate),size(wper(gwiu,2),1),size(frbcom,2),3);
% $$$ wbcv_off_segs = abs(diff(sq(sqrt(sum(wbc_off_segs.^2,4))),1));
% $$$ wbc_on_segs = GetSegs(frbcom,round(wper(gwi,1)-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ wbc_on_segs = reshape(wbc_on_segs,round(6*Session.xyzSampleRate),size(wper(gwi,1),1),size(frbcom,2),3);
% $$$ wbcv_on_segs = abs(diff(sq(sqrt(sum(wbc_on_segs.^2,4))),1));
% $$$ 
% $$$ hbc_off_segs = GetSegs(frbcom,round(hper(ghiu,2)-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
% $$$ hbc_off_segs = reshape(hbc_off_segs,round(6*Session.xyzSampleRate),size(hper(ghiu,2),1),size(frbcom,2),3);
% $$$ hbcv_off_segs = abs(diff(sq(sqrt(sum(hbc_off_segs.^2,4))),1));
% $$$ 
% $$$ 
% $$$ figure,plot(rbcv_on_segs(:,:,2).*Session.xyzSampleRate/10,'r')
% $$$ hold on, plot(hbcv_off_segs(:,:,2).*Session.xyzSampleRate/10,'g')
% $$$ hold on, plot(wbcv_off_segs(:,:,2).*Session.xyzSampleRate/10,'b')
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,plot(rbcv_off_segs(:,:,2).*Session.xyzSampleRate/10,'r')
% $$$ hold on, plot(wbcv_on_segs(:,:,2).*Session.xyzSampleRate/10,'b')
% $$$ 
% $$$ 
% $$$ gkern = gausswin(31);
% $$$ gkern = gkern/sum(gkern);
% $$$ 
% $$$ frbcv_on_segs = zeros(size(rbcv_on_segs));
% $$$ for i=1:size(rbcv_on_segs,2)
% $$$     for j=1:2,
% $$$         frbcv_on_segs(:,i,j) = Filter0(gkern,rbcv_on_segs(:,i,j));    
% $$$     end
% $$$ end
% $$$ 
% $$$ fwbcv_off_segs = zeros(size(wbcv_off_segs));
% $$$ for i=1:size(wbcv_off_segs,2)
% $$$     for j = 1:2,
% $$$         fwbcv_off_segs(:,i,j) = Filter0(gkern,wbcv_off_segs(:,i,j));
% $$$     end
% $$$ end
% $$$ 
% $$$ fhbcv_off_segs = zeros(size(hbcv_off_segs));
% $$$ for i=1:size(hbcv_off_segs,2)
% $$$     for j = 1:2,
% $$$         fhbcv_off_segs(:,i,j) = Filter0(gkern,hbcv_off_segs(:,i,j));
% $$$     end
% $$$ end
% $$$ 
% $$$ figure,plot(frbcv_on_segs(:,:,1).*Session.xyzSampleRate/10,'r')
% $$$ hold on, plot(fhbcv_off_segs(:,:,1).*Session.xyzSampleRate/10,'g')
% $$$ hold on, plot(fwbcv_off_segs(:,:,1).*Session.xyzSampleRate/10,'b')
% $$$ 
% $$$ 
% $$$ figure,plot(frbcv_on_segs(:,:,2).*Session.xyzSampleRate/10,'r')
% $$$ hold on, plot(fhbcv_off_segs(:,:,2).*Session.xyzSampleRate/10,'g')
% $$$ hold on, plot(fwbcv_off_segs(:,:,2).*Session.xyzSampleRate/10,'b')
% $$$ 
% $$$ 
% $$$ SpeedWalkVsRearSig = zeros(size(frbcv_on_segs,1),size(frbcom,2));
% $$$ for tind = 1:size(rvel_on_segs,1),
% $$$     for j = 1:size(frbcom,2),
% $$$         SpeedWalkVsRearSig(tind,j) = ranksum(frbcv_on_segs(tind,:,j),fwbcv_off_segs(tind,:,j));
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,plot(SpeedWalkVsRearSig,'.')
% $$$ 


save([Session.spath.analysis Session.name '.' Session.Maze.name '.' Session.trialName '.REARING.' num2str(randi(10000,1)) '.mat'],'bfccg','-v7.3');

% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %% Load Place Fields
% $$$ Session.Pfs = [];
% $$$ Session = Session.loadPfs;
% $$$ 
% $$$ pf_search.mazeName = 'cof';
% $$$ pf_search.trialName = Session.trialName;
% $$$ pf_search.trackingMarker = Session.trackingMarker;
% $$$ pf_search.stateLabel = 'head';
% $$$ pf_search.spk_shuffle = 'n';
% $$$ pf_search.pos_shuffle = 0;
% $$$ pf_search.numBSiterations = 1;
% $$$ pf_search.numZslices = 1;
% $$$ pf_search.nbins = 50;
% $$$ pf_search.smooth = 0.03;
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %% Display Everything
% $$$ figure(102)
% $$$ set(gcf,'CurrentCharacter','l');
% $$$ unit = 1;
% $$$ while 1,
% $$$     clf
% $$$     %% Rate Map
% $$$     subplot2(6,2,[1,2],1);
% $$$     pf_search.stateLabel = 'head';
% $$$     Pfs = Session.getPfs(pf_search);
% $$$     ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})
% $$$     hold on
% $$$     scatter(wpos(gwi,1),wpos(gwi,2),20,'m','Marker','*')
% $$$     title(num2str(unit))
% $$$     %
% $$$     subplot2(6,2,[1,2],2);
% $$$     pf_search.stateLabel = 'rear';
% $$$     Pfs = Session.getPfs(pf_search);
% $$$     ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})
% $$$     hold on
% $$$     scatter(rpos(gri,1),rpos(gri,2),20,'m','Marker','*')
% $$$ 
% $$$     %% CCGs
% $$$     subplot2(6,2,3,1);
% $$$     bar(tbin/1000,wfccg(:,2,unit));axis tight; grid on
% $$$     %
% $$$     subplot2(6,2,4,1);
% $$$     bar(tbin/1000,wfccg(:,1,unit));axis tight; grid on
% $$$     %
% $$$     subplot2(6,2,3,2);
% $$$     bar(tbin/1000,rfccg(:,1,unit));axis tight; grid on
% $$$     hold on
% $$$     Lines([],[upperConfidenceBoundary(unit),lowerConfidenceBoundary(unit)],'r')
% $$$     Lines([],mean(mean(bfccg(:,:,1,unit),1),2),'g')
% $$$     if max(sq(rfccg(:,1,unit)))<upperConfidenceBoundary(unit), ylim([0,1+upperConfidenceBoundary(unit)]);,end
% $$$     plot(tbin,p95(:,unit),'m')
% $$$     plot(tbin,p05(:,unit),'m')
% $$$     subplot2(6,2,4,2);
% $$$     bar(tbin/1000,rfccg(:,2,unit));axis tight; grid on
% $$$     hold on
% $$$     Lines([],[upperConfidenceBoundary(unit),lowerConfidenceBoundary(unit)],'r')
% $$$     Lines([],mean(mean(bfccg(:,:,1,unit),1),2),'g')
% $$$     if max(sq(rfccg(:,2,unit)))<upperConfidenceBoundary(unit), ylim([0,1+upperConfidenceBoundary(unit)]);,end
% $$$     plot(tbin,p95(:,unit),'m')
% $$$     plot(tbin,p05(:,unit),'m')
% $$$ 
% $$$     %% Rasters
% $$$     subplot2(6,2,5,2);
% $$$     plot(tron(conr==unit)/1000,tonr(conr==unit),'.','markersize',5);axis tight; grid on
% $$$     subplot2(6,2,6,2);
% $$$     plot(troff(coffr==unit)/1000,toffr(coffr==unit),'.','markersize',5);axis tight; grid on
% $$$     subplot2(6,2,5,1);
% $$$     plot(twoff(coffw==unit)/1000,toffw(coffw==unit),'.','markersize',5);axis tight; grid on
% $$$     subplot2(6,2,6,1);
% $$$     plot(twon(conw==unit)/1000,tonw(conw==unit),'.','markersize',5);axis tight; grid on
% $$$ 
% $$$     %% ReportFig
% $$$     %reportfig(102,[Session.filebase '.rearing'],[],['unit: ' num2str(unit)],[],0);
% $$$     
% $$$     %% Figure controls
% $$$     waitforbuttonpress
% $$$     whatkey = get(gcf,'CurrentCharacter');
% $$$     switch double(whatkey)
% $$$       case double('i')
% $$$         unit = input('Enter unit #: ');
% $$$       case double('n')
% $$$         unit = unit+1;
% $$$       case double('p')
% $$$         unit=unit-1;
% $$$       case double('q')
% $$$         return
% $$$     end
% $$$ end
 
% $$$ 
% $$$ figure,
% $$$ subplot(2,2,1),bar(tbin,sum(rssi(:,:,1,2))./size(rssi,1).*100),axis tight,title(gca,'
% $$$ subplot(2,2,3),bar(tbin,sum(rssi(:,:,1,1))./size(rssi,1).*100),axis tight
% $$$ subplot(2,2,2),bar(tbin,sum(rssi(:,:,2,2))./size(rssi,1).*100),axis tight
% $$$ subplot(2,2,4),bar(tbin,sum(rssi(:,:,2,1))./size(rssi,1).*100),axis tight

% $$$ figure,bar(atbin,accg(:,6))


% $$$ %Fig - RCCG
% $$$ x =1;
% $$$ y =1;
% $$$ gsp = 26;
% $$$ for unit = 1:size(Map,1),
% $$$     subplot2(2*round(size(Map,1)/gsp),gsp,y,x);
% $$$     bar(tbin./1000,ffccg(:,1,unit));axis tight; grid on
% $$$     title(num2str(unit));
% $$$     subplot2(2*round(size(Map,1)/gsp),gsp,y+1,x);
% $$$     bar(tbin./1000,ffccg(:,2,unit));axis tight; grid on
% $$$     title(num2str(unit));
% $$$     if x==gsp,x=1;y=y+2;else,x=x+1;end
% $$$ end
% $$$ 
% $$$ %% Part 2 - The TrigRaster Monster
% $$$ inter_rear_dur = [rp(1,1);abs(rp(1:end-1,2)-rp(2:end,1));size(Session.xyz,1)-rp(end)];
% $$$ rear_dur = diff(rp,1,2);
% $$$ rd_thresh = 120; % samples @ xyzSampleRate ~ 1sec
% $$$ ird_thresh = 600;% samples @ xyzSampleRate ~ 5sec
% $$$ %% Selection of Rears
% $$$ % Only onsets and offset whith out another following within 5 seconds
% $$$ 
% $$$ rt_on  = rp(inter_rear_dur(1:end-1)>ird_thresh&rear_dur>rd_thresh,1);
% $$$ rt_off = rp(inter_rear_dur(2:end)>ird_thresh&rear_dur>rd_thresh,2);
% $$$ rpon_Res  = round((rt_on-1.05)/Session.xyzSampleRate*Session.lfpSampleRate);
% $$$ rpoff_Res = round((rt_off-1.05)/Session.xyzSampleRate*Session.lfpSampleRate);
% $$$ 
% $$$ total_rear_count = size(rp,1);
% $$$ unit_rear_count = zeros(size(Map,1),2);
% $$$ unit_rear_firingrate = cell(size(Map,1),2);
% $$$ unit_rear_position = cell(size(Map,1),2);
% $$$ unit_rear_raster = cell(size(Map,1),2); 
% $$$ unit_rear_rtind = cell(size(Map,1),2); 
% $$$ 
% $$$ frWindow1 = 2; % sec
% $$$ frWindow2 = 2; % sec
% $$$ fr_on = zeros(size(rt_on,1),size(Map,1));
% $$$ fr_off = zeros(size(rt_off,1),size(Map,1));
% $$$ for i = 1:size(rt_on,1),
% $$$     [sRes sind] = SelectPeriods(Res,round([rpon_Res(i),rpon_Res(i)]+[-Session.lfpSampleRate*frWindow1,Session.lfpSampleRate*frWindow2]),'d',1,1);
% $$$     sClu = Clu(sind);
% $$$     fr_on(i,unique(sClu)) = FiringRate(sRes,sClu,[],Session.lfpSampleRate);
% $$$ end
% $$$ for i = 1:size(rt_off,1),
% $$$     [sRes sind] = SelectPeriods(Res,round([rpoff_Res(i),rpoff_Res(i)]+[-Session.lfpSampleRate*frWindow1,Session.lfpSampleRate*frWindow2]),'d',1,1);
% $$$     sClu = Clu(sind);
% $$$     fr_off(i,unique(sClu)) = FiringRate(sRes,sClu,[],Session.lfpSampleRate);
% $$$ end
% $$$ fr_on(isnan(fr_on))=0;
% $$$ fr_off(isnan(fr_off))=0;
% $$$ 
% $$$ for unit = 1:size(Map,1),
% $$$     tic
% $$$     %% Selection of rears
% $$$     [frs_on,ifrs_on]=sort(fr_on(:,unit),'descend');
% $$$     hspkrind_on = ifrs_on(find(frs_on~=0));
% $$$     bins_on = prctile(fr_on(hspkrind_on,unit),[0:30:100]);
% $$$     [count_on,ind_on] = histcI(fr_on(hspkrind_on,unit),bins_on);
% $$$     [frs_off,ifrs_off]=sort(fr_off(:,unit),'descend');
% $$$     hspkrind_off = ifrs_off(find(frs_off~=0));
% $$$     bins_off = prctile(fr_off(hspkrind_off,unit),[0:30:100]);
% $$$     [count_off,ind_off] = histcI(fr_off(hspkrind_off,unit),bins_off);
% $$$     unit_rear_count(unit,1) = length(hspkrind_on);
% $$$     unit_rear_count(unit,2) = length(hspkrind_off);
% $$$     %% Raster Plot
% $$$     if length(Res(Clu==unit))>1,
% $$$         if size(hspkrind_on,1)>1,
% $$$             [unit_rear_raster{unit,1},unit_rear_rtind{unit,1},~] = TrigRasters(rpon_Res(hspkrind_on)  ,10000,Res(Clu==unit),Clu(Clu==unit),Session.lfpSampleRate,0,ind_on);
% $$$         end
% $$$         if size(hspkrind_off,1)>1,
% $$$             [unit_rear_raster{unit,2},unit_rear_rtind{unit,2},~] = TrigRasters(rpoff_Res(hspkrind_off),10000,Res(Clu==unit),Clu(Clu==unit),Session.lfpSampleRate,0,ind_off);
% $$$         end
% $$$     end
% $$$     toc
% $$$     unit
% $$$ end
% $$$ 
% $$$ rrwindow = 125;%samples
% $$$ overlap = 62;
% $$$ 
% $$$ wind = -2500:overlap:2500;
% $$$ 
% $$$ %% Part 4 Rearing and Theta
% $$$ 
% $$$ % trim to fit window for matrix reshaping
% $$$ window = 60;%samples
% $$$ numSamplesToTrim = mod(size(rb_com_xyz,1),window);
% $$$ 
% $$$ xyzsegs = reshape(sq(rb_com_xyz(1:end-numSamplesToTrim,1,:)),window,[],size(rb_com_xyz,2));
% $$$ distsegs = sqrt(sum(xyzsegs.^2,3));
% $$$ velsegs = abs(diff(distsegs));
% $$$ accsegs = abs(diff(velsegs));
% $$$ accsegs = diff(velsegs);
% $$$ 
% $$$ average_speed_segs = mean(velsegs,1)';
% $$$ average_accel_segs = mean(accsegs,1)';
% $$$ 
% $$$ %tsegs = d2t(average_speed_segs,window/Session.xyzSampleRate,0);
% $$$ 
% $$$ tsegs = [1:size(average_speed_segs,1)].*window/Session.xyzSampleRate;
% $$$ tsegs = tsegs-diff(tsegs(1:2))/2;
% $$$ 
% $$$ % Compute lfp and theta power
% $$$ %% missing spectral calculations
% $$$ 
% $$$ timeToTrim = numSamplesToTrim/Session.xyzSampleRate;
% $$$ thpow = sq(mean(y(:,f<12&f>6,:),2));
% $$$ thpow = thpow(1:(find(t<t(end)-timeToTrim,1,'last')-1),:);
% $$$ 
% $$$ 
% $$$ tt = t(1:(find(t<t(end)-timeToTrim,1,'last')-1));
% $$$ 
% $$$ 
% $$$ thp = zeros(size(average_speed_segs,1),32);
% $$$ tind = zeros(2,1);
% $$$ for i = 1:size(average_speed_segs,1)-1,
% $$$     tind = zeros(2,1);
% $$$     tind(1) = find(tt<tsegs(i),1,'last');
% $$$     tind(2) = find(tt>tsegs(i),1,'first');
% $$$     thp(i,:) = mean(thpow(tind(tind~=0),:),1);
% $$$ end
% $$$  
% $$$ 
% $$$ % speed vs theta power
% $$$ plot(tsegs,unity(average_speed_segs),tsegs,unity(thp(:,18)),'.');
% $$$ plot(unity(average_speed_segs(average_speed_segs~=0)),unity(thp(average_speed_segs~=0,18)))
% $$$ cov(unity(average_speed_segs(average_speed_segs~=0)),unity(thp(average_speed_segs~=0,18)))
% $$$ 
% $$$ rears_per = Session.Bhv.getState('rear').state/Session.xyzSampleRate;
% $$$ rthp = [];
% $$$ ravs = [];
% $$$ for i = 1:size(rears_per,1),
% $$$     rthp = cat(1,rthp,thp(rears_per(i,1)<tsegs&rears_per(i,2)>tsegs,:));
% $$$     ravs = cat(1,ravs,average_speed_segs(rears_per(i,1)<tsegs&rears_per(i,2)>tsegs));
% $$$ end
% $$$ 
% $$$ walks_per = Session.Bhv.getState('walk').state/Session.xyzSampleRate; 
% $$$ 
% $$$ walks_per = walks_per+repmat([1,-1],size(walks_per,1),1);
% $$$ 
% $$$ wthp = [];
% $$$ wavs = [];
% $$$ for i = 1:size(walks_per,1),
% $$$     wthp = cat(1,wthp,thp(walks_per(i,1)<tsegs&walks_per(i,2)>tsegs,:));
% $$$     wavs = cat(1,wavs,average_speed_segs(walks_per(i,1)<tsegs&walks_per(i,2)>tsegs));
% $$$ end
% $$$ 
% $$$ 
% $$$ plot(unity(ravs(ravs~=0)),unity(rthp(ravs~=0,18)),'.',unity(wavs(wavs~=0)),unity(wthp(wavs~=0,18)),'.')
% $$$ 
% $$$ 
% $$$  
% $$$ %% Part 5 Theta rear onset offset
% $$$ 
% $$$ 
% $$$ rthp = sq(median(y{1}(:,fr<12&fr>5,:,:,:),2));
% $$$ 
% $$$ 
% $$$ rxyz = GetSegs(Session.xyz(:,7,3),rper-30,60,[]);
% $$$ rxyz = reshape(rxyz,[],size(rper,1),2);
% $$$ rvel = sq((rxyz(end,:,:)-rxyz(1,:,:))/60*Session.xyzSampleRate);
% $$$ 
% $$$ 
% $$$ rxyz = GetSegs(Session.xyz,round(rper-3*Session.xyzSampleRate),round(6*Session.xyzSampleRate),[]);
% $$$ rxyz = reshape(rxyz,round(6*Session.xyzSampleRate),size(rper,1),2,size(Session.xyz,2),size(Session.xyz,3));
% $$$ 
% $$$ rang = GetSegs(Session.ang(:,:,:,1:2),round(rper-3*Session.xyzSampleRate),round(6*Session.xyzSampleRate),[]);
% $$$ rang = reshape(rang,round(6*Session.xyzSampleRate),size(rper,1),2,size(Session.xyz,2),size(Session.xyz,2),2);
% $$$ rang(isnan(rang))=0;
% $$$ 
% $$$ rfet = reshape(cat(4,rang(:,:,:,4,5,:),rang(:,:,:,5,7,:)),round(6*Session.xyzSampleRate),size(rper,1),2,2,2);
% $$$ %rfet = rxyz;
% $$$ 
% $$$ rswin = 12;
% $$$ rspeed = zeros(size(yr,1),size(rfet,2),size(rfet,3),size(rfet,4));
% $$$ racc = zeros(size(yr,1),size(rfet,2),size(rfet,3),size(rfet,4));
% $$$ trs = [1:size(rfet,1)]/Session.xyzSampleRate;
% $$$ tind = zeros(rswin,1);
% $$$ for i = 1:size(tr,1),
% $$$     tind = zeros(rswin,1);
% $$$     if ~isempty(find(trs<tr(i),rswin/2,'last')),
% $$$         tind(1:rswin/2) = find(trs<tr(i),rswin/2,'last');
% $$$     end
% $$$     if ~isempty(find(trs>tr(i),rswin/2,'first')),
% $$$         tind(rswin/2+1:rswin) = find(trs>tr(i),rswin/2,'first');
% $$$     end
% $$$     rspeed(i,:,:,:) = median(abs(sq(sqrt(sum(diff(rfet(tind(tind~=0),:,:,:,:),1).^2,5)))),1);
% $$$     racc(i,:,:,:) = median(diff(abs(sq(diff(sqrt(sum(rfet(tind(tind~=0),:,:,:,:).^2,5)),1))),1));
% $$$ end
% $$$ rdur = diff(rper,1,2)/Session.xyzSampleRate;
% $$$ 
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ figure
% $$$ % $$$ hold on
% $$$ % $$$ for i=1:size(yr,1),plot(rspeed(i,rdur>1,1,7),rthp(i,rdur>1,1,18),'.',rspeed(i,rdur>1,2,7),rthp(i,rdur>1,2,18),'.'),pause(.4),end
% $$$ 
% $$$ figure
% $$$ clf
% $$$ for i=20:31,
% $$$     plot(rspeed(i,:,1,2),rthp(i,:,1,30),'.',rspeed(i,:,2,2),rthp(i,:,2,30),'.'),
% $$$     title(num2str(i))
% $$$     xlim([0,.14])
% $$$     ylim([0,2000])    
% $$$     waitforbuttonpress
% $$$ end
% $$$ 
% $$$ % $$$ figure
% $$$ % $$$ 
% $$$ % $$$ clf
% $$$ % $$$ for i=1:size(yr,1),
% $$$ % $$$     plot(racc(i,:,1,7),rthp(i,:,1,18),'.',racc(i,:,2,7),rthp(i,:,2,18),'.')
% $$$ % $$$     title(num2str(i))
% $$$ % $$$     xlim([-0.5,0.5])
% $$$ % $$$     ylim([0,6000])
% $$$ % $$$     waitforbuttonpress
% $$$ % $$$ end
% $$$ 
% $$$ inter_rear_dur = [rper(1,1);abs(rper(1:end-1,2)-rper(2:end,1));size(Session.xyz,1)-rper(end)];
% $$$ rd_thresh = 1; % samples @ xyzSampleRate ~ 1sec
% $$$ ird_thresh = 360;% samples @ xyzSampleRate ~ 3sec
% $$$ rt_on  = rper(inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,1);
% $$$ rt_off = rper(inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,2);
% $$$ 
% $$$ rfitcof = zeros(size(yr,1),2,2,size(rspeed,4),32);
% $$$ rrankcof = zeros(size(yr,1),2,2,size(rspeed,4),32);
% $$$ rcorcof = zeros(size(yr,1),2,2,size(rspeed,4),32);
% $$$ rfdr = zeros(size(yr,1),2,size(rspeed,4),32);
% $$$ for chan = 1:32
% $$$     for mar = 1:size(rspeed,4),        
% $$$         for i=1:size(yr,1),
% $$$             rfitcof(i,1,:,mar,chan) = polyfit(rspeed(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,1,mar),rthp(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,1,chan),1);
% $$$             rfitcof(i,2,:,mar,chan) = polyfit(rspeed(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,2,mar),rthp(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,2,chan),1);
% $$$             rrankcof(i,1,:,mar,chan) = RankCorrelation(rspeed(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,1,mar),rthp(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,1,chan));
% $$$             rrankcof(i,2,:,mar,chan) = RankCorrelation(rspeed(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,2,mar),rthp(i,inter_rear_dur(1:end-1)>ird_thresh&rdur>rd_thresh,2,chan));
% $$$             [tcor,pval] = corrcoef(rspeed(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,1,mar),rthp(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,1,chan));            
% $$$             rcorcof(i,1,1,mar,chan) = tcor(2);
% $$$             rcorcof(i,1,2,mar,chan) = pval(2);
% $$$             [tcor,pval] = corrcoef(rspeed(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,2,mar),rthp(i,inter_rear_dur(2:end)>ird_thresh&rdur>rd_thresh,2,chan));
% $$$             rcorcof(i,2,1,mar,chan) = tcor(2);
% $$$             rcorcof(i,2,2,mar,chan) = pval(2);
% $$$         end        
% $$$         rfdr(:,1,mar,chan) = fdr_bh(rcorcof(:,1,2,mar,chan));
% $$$         rfdr(:,2,mar,chan) = fdr_bh(rcorcof(:,2,2,mar,chan));
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,imagesc(sq(rcorcof(:,2,2,2,:))'.*sq(rfdr(:,2,2,:))')
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ ch = 1;
% $$$ ma = 1;
% $$$ for chan = 1:2:32
% $$$     for mar = 1:2:size(rspeed,4),
% $$$         subplot2(5,16,ma,ch);
% $$$         %plot(sq(rfitcof(:,1,:,mar,chan)),'.');
% $$$         plot(sq(rcorcof(:,2,mar,chan)),'.');
% $$$         hold on,
% $$$         Lines(29,[],'r');
% $$$         ma = ma+1;
% $$$     end
% $$$     ma = 1;
% $$$     ch = ch+1;
% $$$ end
% $$$ %ForAllSubplots('xlim([0,55])')
% $$$ %ForAllSubplots('ylim([-1000,2500])')
% $$$ 
% $$$ % For angles
% $$$ figure
% $$$ hold on
% $$$ ch = 1;
% $$$ ma = 1;
% $$$ for chan = 1:2:32,
% $$$     for mar = 1:2,
% $$$         subplot2(2,16,ma,ch);
% $$$         %plot(sq(rfitcof(:,1,:,mar,chan)),'.');
% $$$         plot(sq(rcorcof(:,1,mar,chan)),'.');
% $$$         hold on,
% $$$         Lines(29,[],'r');
% $$$         ma = ma+1;
% $$$     end
% $$$     ma = 1;
% $$$     ch = ch+1;
% $$$ end
% $$$ %ForAllSubplots('xlim([0,55])')
% $$$ %ForAllSubplots('ylim([-1000,2500])')
% $$$ 
% $$$ 
% $$$ %% Part 6 Theta rear vs theta walk
% $$$ 
% $$$ Session = MTASession('jg05-20120310');
% $$$ Session = MTATrial(Session,'all');
% $$$ 
% $$$ lfp = Session.loadlfp(65:96);
% $$$ wlfp = WhitenSignal(lfp,[],1);
% $$$ 
% $$$ rper = Session.Bhv.getState('rear').state;
% $$$ wper = Session.Bhv.getState('walk').state;
% $$$ frange = [1,90];
% $$$ nffts = 2^10;
% $$$ windows =2^9;
% $$$ yr = [];
% $$$ yw = [];
% $$$ t = [];
% $$$ f = [];
% $$$ for b = 1:2,
% $$$     rwlfp = GetSegs(wlfp,round((rper(:,b)./Session.xyzSampleRate).*1250),10*1250,[]);
% $$$     wwlfp = GetSegs(wlfp,round((wper(:,b)./Session.xyzSampleRate).*1250),10*1250,[]);
% $$$     for j = 1:32;
% $$$         for i = 1:size(rwlfp,2),
% $$$             [yr(:,:,i,b,j),f,t] = mtchglong(rwlfp(:,i,j),nffts,Session.lfpSampleRate,windows,0.75*windows,[],[],[],frange);
% $$$         end
% $$$         for i = 1:size(wwlfp,2),
% $$$             [yw(:,:,i,b,j),f,t] = mtchglong(wwlfp(:,i,j),nffts,Session.lfpSampleRate,windows,0.75*windows,[],[],[],frange);
% $$$         end        
% $$$     end
% $$$ end
% $$$ t = t+diff(t(1:2))/2;
% $$$ 
% $$$ rb_com_xyz = zeros(size(Session.xyz,1),2,3);
% $$$ rb = Session.Model.rb({'head_back','head_left','head_front','head_right'});
% $$$ rb_com_xyz(:,1,:) = sq(Session.com(rb));
% $$$ rb = Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
% $$$ rb_com_xyz(:,2,:) = sq(Session.com(rb));
% $$$ 
% $$$ %% Needs rear walk selection 
% $$$ 
% $$$ rper = rper(goodRearInd==1,:);
% $$$ rdur = rdur(goodRearInd==1);
% $$$ 
% $$$ wper = wper(goodWalkInd==1,:);
% $$$ wdur = wdur(goodWalkInd==1);
% $$$ 
% $$$ rthp = sq(median(yr(:,f<12&f>5,:,:,:),2));
% $$$ wthp = sq(median(yw(:,f<12&f>5,:,:,:),2));
% $$$ 
% $$$ wthp = sq(wthp(:,goodWalkInd==1,1,:));
% $$$ yw = sq(yw(:,:,goodWalkInd==1,1,:));
% $$$ wgmp = sq(median(yw(:,f<80&f>60,:,:),2));
% $$$ wspeed = sq(wspeed(:,goodWalkInd==1,1,:));
% $$$ wacc = sq(wacc(:,goodWalkInd==1,1,:));
% $$$ 
% $$$ rthp = sq(rthp(:,goodRearInd==1,1,:));
% $$$ yr = sq(yr(:,:,goodRearInd==1,1,:));
% $$$ rgmp = sq(median(yr(:,f<80&f>60,:,:),2));
% $$$ rspeed = sq(rspeed(:,goodRearInd==1,1,:));
% $$$ racc = sq(racc(:,goodRearInd==1,1,:));
% $$$ 
% $$$ %% start here fix indexing error
% $$$ 
% $$$ for i = 1:size(wper,1),
% $$$     wthp(t>=wdur(i),i,:) = 0;
% $$$     yw(t>=wdur(i),:,i,:) = 0;
% $$$     wgmp(t>=wdur(i),i,:) = 0;
% $$$     wspeed(t>=wdur(i),i,:) = 0;
% $$$     wacc(t>=wdur(i),i,:) = 0;
% $$$ end     
% $$$ 
% $$$ for i = 1:size(rper,1),
% $$$     rthp(t>=rdur(i),i,:) = 0;
% $$$     yr(t>=rdur(i),:,i,:) = 0;
% $$$     rgmp(t>=rdur(i),i,:) = 0;
% $$$     rspeed(t>=rdur(i),i,:) = 0;
% $$$     racc(t>=rdur(i),i,:) = 0;
% $$$ end     
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ for i = 1:size(rspeed,2)
% $$$     if size(wspeed(wspeed(:,i,1)~=0,i,1),1)~=size(wthp(wthp(:,i,18)~=0,i,18),1),continue,end
% $$$     plot(rspeed(rspeed(:,i,1)~=0,i,1),rthp(rthp(:,i,18)~=0,i,18),'.',wspeed(wspeed(:,i,1)~=0,i,1),wthp(wthp(:,i,18)~=0,i,18),'.')
% $$$ end
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ for i = 1:size(rspeed,2)
% $$$     if size(wspeed(wspeed(:,i,2)~=0,i,2),1)~=size(wthp(wthp(:,i,18)~=0,i,18),1),continue,end
% $$$     plot(rspeed(rspeed(:,i,2)~=0,i,2),rthp(rthp(:,i,18)~=0,i,18),'.',wspeed(wspeed(:,i,2)~=0,i,2),wthp(wthp(:,i,18)~=0,i,18),'.')
% $$$ end
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ for i = 1:size(rspeed,2)
% $$$     if size(wspeed(wspeed(:,i,3)~=0,i,3),1)~=size(wthp(wthp(:,i,18)~=0,i,18),1),continue,end
% $$$     plot(rspeed(rspeed(:,i,3)~=0,i,3),rthp(rthp(:,i,18)~=0,i,18),'.',wspeed(wspeed(:,i,3)~=0,i,3),wthp(wthp(:,i,18)~=0,i,18),'.')
% $$$ end
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ for i = 1:size(rspeed,2)
% $$$     if size(wspeed(wspeed(:,i,1)~=0,i,1),1)~=size(wgmp(wgmp(:,i,18)~=0,i,18),1),continue,end
% $$$     plot(rspeed(rspeed(:,i,1)~=0,i,1),rgmp(rgmp(:,i,18)~=0,i,18),'.',wspeed(wspeed(:,i,1)~=0,i,1),wgmp(wgmp(:,i,18)~=0,i,18),'.')
% $$$ end
% $$$ 
% $$$ figure
% $$$ hold on
% $$$ for i = 1:size(rspeed,2)
% $$$     if size(wspeed(wspeed(:,i,3)~=0,i,3),1)~=size(wgmp(wgmp(:,i,18)~=0,i,18),1),continue,end
% $$$     plot(rspeed(rspeed(:,i,3)~=0,i,3),rgmp(rgmp(:,i,18)~=0,i,18),'.',wspeed(wspeed(:,i,3)~=0,i,3),wgmp(wgmp(:,i,18)~=0,i,18),'.')
% $$$ end
% $$$ 
% $$$ plot(rspeed(rspeed~=0),rthp(rthp~=0),'.',wspeed(wspeed~=0),wgmp(wgmp~=0),'.')
% $$$ 
% $$$ 
% $$$ rcorcof = zeros(size(yr,1),2,size(rspeed,3),32);
% $$$ rfdr = zeros(size(yr,1),size(rspeed,3),32);
% $$$ for chan = 1:32
% $$$     for mar = 1:size(rspeed,3),        
% $$$         for i=1:size(t,1),
% $$$             [tcor,pval] = corrcoef(rspeed(i,:,mar),rthp(i,:,chan));            
% $$$             rcorcof(i,1,mar,chan) = tcor(2);
% $$$             rcorcof(i,2,mar,chan) = pval(2);
% $$$         end        
% $$$         rfdr(:,mar,chan) = fdr_bh(rcorcof(:,2,mar,chan));
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ figure,imagesc(sq(rcorcof(:,2,3,:))')
% $$$ 
% $$$ 
% $$$ %rfitcof = zeros(size(yr,1),2,size(rspeed,3),32);
% $$$ %rfitcof(i,:,mar,chan) = polyfit(rspeed(i,:,mar),rthp(i,:,chan),1);
% $$$ 
% $$$ 
% $$$ 
% $$$ 

