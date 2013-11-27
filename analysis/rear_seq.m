function rear_seq(Session)

Session = MTASession('jg05-20120310',{'CluRes','ang','ufr'});

rper = Session.Bhv.getState('rear').state;
[~,rfet] = rear(Session,'com');
rpos = sq(Session.xyz(rper(:,1),Session.Model.gmi(Session.trackingMarker),[1,2]));
rdur = diff(rper,1,2);
fet = GetSegs(Session.ang(:,3,4,2),round(rper-3.*Session.xyzSampleRate),round(6*Session.xyzSampleRate),0);
fet = reshape(fet,round(6*Session.xyzSampleRate),size(rper,1),2);
gri = find(max(fet(1:300,:,1))<0&rdur'>1.5);
griu = find(max(fet(450:end,:,2))<0&rdur'>1.5);

binSize = 100;%ms
halfBins = 50;
normalization = 'hz';


Trains2CCG
binSize = 100;%ms
halfBins = 50;
normalization = 'hz';

per_on  = rper(gri,1);
per_off = rper(griu,2);
poni  = round((per_on-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
poffi = round((per_off-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
[fccg tbin ] = Trains2CCG({poni,poffi,Session.res},{1,2,Session.clu},binSize,halfBins,Session.lfpSampleRate,normalization);
uClu = unique(Session.clu);
sccg = zeros(size(fccg,1),2,size(Session.map,1));
sccg(:,:,uClu) = sq(fccg(:,1:2,3:end));
ffccg = zeros(size(fccg,1),2,size(Session.map,1));
for k=1:size(Session.map,1)
    for l=1:2
        ffccg(:,l,k) = Filter0(gausswin(7)/sum(gausswin(7)),sccg(:,l,k));
    end
end
rfccg = ffccg;

[u,s,v]=svd(sq(fccg(:,1,3:end)));
imagesc(v')
figure,imagesc(u')

[u,s,v]=svd(sq(rfccg(:,1,:)));
% $$$ 
% $$$ figure,imagesc(v')
% $$$ figure,imagesc(u')
% $$$ figure,imagesc(s')
% $$$ 
% $$$ dcm_obj = datacursormode(gcf);
% $$$ ci_obj = getCursorInfo(dcm_obj);
% $$$ units = ci_obj.Position;
% $$$ figure(21331)
% $$$ subplot(211)
% $$$ bar(tbin,sq(rfccg(:,1,units(1))));
% $$$ subplot(212)
% $$$ bar(tbin,sq(rfccg(:,1,units(2))));

urfccg = u*sq(rfccg(:,1,:));
vrfccg = sq(rfccg(:,1,:))*v;


% $$$ figure(21332)
% $$$ imagesc(vrfccg')
% $$$ 
% $$$ subplot(211)
% $$$ bar(tbin,sq(rfccg(:,1,units(1))));
% $$$ subplot(212)
% $$$ bar(tbin,sq(rfccg(:,1,units(2))));




fpc = u(2,:)*sq(rfccg(:,1,:));
spc = u(3,:)*sq(rfccg(:,1,:));

[~,rcells] = sort(spc,'descend');


uvrfccg = u*vrfccg(:,i);
[~,rcells] = sort(spc,'descend');

bar(tbin,sq(rfccg(:,1,2)))




unrcg = unity(sq(rfccg(:,1,:)));
[u,s,v]=svd(unrcg);
urcg = u(25:75,1)'*unrcg(25:75,:);
[~,rcellsupon] = sort(urcg);
urcg = u(25:75,2)'*unrcg(25:75,:);
[~,rcellsdownon] = sort(urcg);
rcellsdownon = flipdim(rcellsdownon,1);

rsenson = cat(2,rcellsupon(1:40),rcellsdownon(1:40));
rsenson = unique(rsens);


unrcg = unity(sq(rfccg(:,2,:)));
[u,s,v]=svd(unrcg);
urcg = u(25:75,1)'*unrcg(25:75,:);
[~,rcellsupon] = sort(urcg);
urcg = u(25:75,2)'*unrcg(25:75,:);
[~,rcellsdownon] = sort(urcg);
rcellsdownon = flipdim(rcellsdownon,1);

rsensoff = cat(2,rcellsupon(1:40),rcellsdownon(1:40));
rsensoff = unique(rsens);

rsens = cat(2,rsenson,rsensoff);
rsens = unique(rsens);


figure,
for i =1:length(rsens),
subplotfit(i,length(rsens))
%bar(tbin,vrcg(:,i))
%bar(tbin,unrcg(:,rcells(i)))
%bar(tbin,uvrfccg(:,i))
bar(tbin,sq(rfccg(:,1,rsens(i))));
end
ForAllSubplots('axis tight')


Session = Session.load_ufr();


ufr = Session.ufr(1:200000,rsens)+1;

nbin = 1200;
trim = mod(size(ufr,1),nbin);

tpvec = zeros(nbin,round((size(ufr,1)-nbin-1)/nbin-1),size(ufr,2));
pvec = [];
for i = 0:nbin-1,
tpvec(:,:,:) = sq(reshape(ufr(1+i:end-trim-nbin+i,:),nbin,[],size(ufr,2)));
abin = tpvec(1:nbin/3,:,:);
bbin = tpvec(nbin/3+1:2*nbin/3,:,:);
cbin = tpvec(2*nbin/3+1:nbin,:,:);
pvec(1+i,:,:) = (sum(bbin,1)./sum(abin,1)+sum(bbin,1)./sum(cbin,1))./(sum(abin,1)+sum(bbin,1));
end

gvec = reshape(pvec,[],64);
