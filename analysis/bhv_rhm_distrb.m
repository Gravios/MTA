function bhv_rhm_distrb(Trial,varargin)
[feature] = DefaultArgs(varargin,{'v'});

sname = 'jg05-20120310';
tname = 'all';
marker = 'spine_lower';
stc_mode = 'auto_wbhr';


% sname = 'er06-20130614';
% tname = 'all-cof';
% stc_mode = 'manual_tmknsrw';
% stc_mode = 'qda_filtp5';


Trial = MTATrial(sname,tname);
Trial.stc.updateMode(stc_mode);
Trial.stc.load;

% $$$ ang = Trial.ang.copy;
% $$$ ang.load(Trial);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.2,xyz.sampleRate));

switch feature
  case 'z'
    vel = MTADxyz('data',xyz(:,'head_front',3),'sampleRate',xyz.sampleRate);
    vnn = nniz(vel);
  case 'v'
    vel = xyz.vel(marker,[1,2]);
    vel.data = log10(vel.data);
    vnn = nniz(vel);
end


[rhm,fs] = fet_rhm(Trial,xyz.sampleRate,'wspectral');
rhm.data  = log10(rhm.data);
rhm.data(isinf(rhm.data))=nan;


vlim =prctile(vel(Trial.stc{'w'}),[2,98]);



s = 'w';
srhm = rhm(Trial.stc{s},:);
svel = vel(Trial.stc{s},:);

testrhm  = (srhm-repmat(nanmean(rhm(vnn,:)),[size(srhm,1),1]))./repmat(nanstd(rhm(vnn,:)),[size(srhm,1),1]);
%testrhm  = srhm;



vbins = 50;
vedgs = linspace(vlim(1),vlim(2),vbins);

[~,vbs] = histc(vel(Trial.stc{s}),vedgs);
%vel.data = repmat(vel.data,[1,rhm.size(2)]);

mrv = [];
srv = [];
for f =1:rhm.size(2),
    mrv(:,f) = accumarray(vbs(vbs~=0),testrhm(vbs~=0,f),[vbins,1],@nanmedian);
    srv(:,f) = accumarray(vbs(vbs~=0),testrhm(vbs~=0,f),[vbins,1],@nanstd);
end
%figure,imagesc(vedgs,fs,mrv'),axis xy
%figure,imagesc(vedgs,fs,srv'),axis xy
figure,figH = imagesc(vedgs,fs,mrv'./srv');axis xy,%caxis([0,4])
figure,figH = imagesc(vedgs,fs,mrv');axis xy,%caxis([0,4])


spowa = MTADlfp('data',log10(median(10.^rhm(:,fs>7&fs<11),2)./sum(10.^rhm(:,:),2)),'sampleRate',Trial.xyz.sampleRate);

figure,hist2([vel(vnn),clip(spowa(vnn),-5.8,.2)],100,100)



splim = prctile(spowa(Trial.stc{'t'},:),[2,98]);
vlim =prctile(vel(Trial.stc{'t'}),[2,98]);

s = 'w';
vbins = 25;
vedgs = linspace(vlim(1),vlim(2),vbins);
spbins = 25;
spedgs = linspace(splim(1),splim(2),spbins);
figure,hist2([vel(Trial.stc{s}),clip(spowa(Trial.stc{s}),-5.8,.2)],vedgs,spedgs)
title(Trial.stc{s}.label);
xlabel('height cm/s')
%xlabel('height cm/s')
ylabel('log10 rhm power 7-11 hz')

wvel = vel.data;
wvel = WhitenSignal(wvel,[],1,ARmodel);

[ys,fs,ts,phi] = mtchglong(vel.data,2^8,xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);

% [ys,fs,ts,phi,~,yn] = mtcsdlong(vel.data,2^8,xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);
% for sigs = 1:size(wvel,2),
%         ys(:,:,sigs,sigs) = mtcsglong(vel(:,sigs),2^8,xyz.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
% end

ts = ts+(2^6)/vel.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),2^7/xyz.sampleRate].*ssr);
szy = size(ys);
ysd = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),...
             'sampleRate',ssr);



wfet(:,1) = mean(log10(abs(yst(:,fs<10,1,1))),2)./mean(log10(abs(yst(:,fs>15,1,1))),2);
wfet(:,2) = mean(log10(abs(yst(:,fs<10,2,2))),2)./mean(log10(abs(yst(:,fs>15,2,2))),2);

wfet = mean(log10(abs(ysd(:,fs>25,2,2))),2);

wcoh = abs(yn(:,:,1,2))./sqrt((ys(:,:,1,1).*ys(:,:,2,2)));
% $$$ 
% $$$ sbins = 25;
% $$$ sedges = [-5,-2];
% $$$ %sedges = [50,160];
% $$$ 
% $$$ vbins = 25;
% $$$ vedges = [0,2];
% $$$ 
% $$$ wper = Trial.stc{'w'}.copy;
% $$$ wper.cast('TimeSeries');
% $$$ wper.resample(yad);
% $$$ 
% $$$ %spow = xyz(:,7,3);
% $$$ spow = spowa;
% $$$ spow = clip(spow,sedges(1),sedges(2));
% $$$ sind = spow>sedges(1)&spow<sedges(2);
% $$$ vlog = clip(log10(v.data),vedges(1),vedges(2));
% $$$ vind = vlog>vedges(1)&vlog<vedges(2)&~isnan(vlog);
% $$$ aind = sind&vind&wper.data;
% $$$ 
% $$$ [spow_count,shind] = histc(spow(aind),sedges(1):abs(diff(sedges))/sbins:sedges(2));
% $$$ [v_count,vhind] = histc(vlog(aind),vedges(1):abs(diff(vedges))/vbins:vedges(2));
% $$$ 
% $$$ %tpow = log10(mean(yld(aind,fl>=6&fl<=12,3),2));
% $$$ tpow = log10(mean(yld(aind,fl>=4&fl<=16,1,1),2));
% $$$ %tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
% $$$ tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$ 
% $$$ %A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
% $$$ %figure,
% $$$ %subplot(133),imagescnan({vedges,sedges,A'},[-5,-3],[],1,[0,0,0]),axis xy,
% $$$ %clf
% $$$ %imagescnan({vedges,sedges,A'},[1.2,1.9],[],1,[0,0,0]),axis xy,
% $$$ %subplot(122),imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ % CA1 LM 81;
% $$$ % DG  G  85;
% $$$ % CA3 ?  95;
% $$$ chan = find(chans == 71);
% $$$ numIter = 10000;
% $$$ %tpow = log10(mean(yld(aind,fl>6&fl<12,chan),2));
% $$$ %tpow = log10(mean(yld(aind,fl<=4,chan),2));
% $$$ tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
% $$$ %tpow = log10(mean(yld(aind,fh>50&fh<80,chan),2));
% $$$ 
% $$$ B=nan(vbins,sbins,numIter);
% $$$ A=nan(vbins,sbins,numIter);
% $$$ %S=nan(vbins,sbins,numIter);
% $$$ tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$ B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
% $$$ A(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
% $$$ %S(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
% $$$ for i = 2:numIter,
% $$$ tpow = tpow(randperm(numel(tpow)));
% $$$ tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$ A(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
% $$$ %S(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
% $$$ end
% $$$ 
% $$$ AS = sort(A,3);
% $$$ P = 1./sum(repmat(A(:,:,1),[1,1,numIter])>A,3);
% $$$ P(isinf(P)) = nan;
% $$$ 
% $$$ SIG = P<=0.0002;
% $$$ ASIG = A; 
% $$$ ASIG(~SIG)=nan;
% $$$ ASIG(B<10)=nan;
% $$$ Aclims = [prctile(ASIG(~isnan(ASIG)),5),prctile(ASIG(~isnan(ASIG)),95)];
% $$$ 
% $$$ figure
% $$$ 
% $$$ subplot(131),imagescnan({vedges,sppedges,B'./yad.sampleRate},[],[],1,[0,0,0]),axis xy,
% $$$ title('Occupancy in seconds')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ subplot(132),imagescnan({vedges,sedges,   A(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
% $$$ title('Mean Power 1-4Hz given 10Hz Osc. Power dist(SU,HB) VS Vel(HF)')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ subplot(133),imagescnan({vedges,sedges,ASIG(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
% $$$ title('P<0.0002 and Occupancy > 1.33 seconds')
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ ticks_lin2log(gca,'x')
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ flim=[1,4;6,12;20,27;30,40];
% $$$ %flim=[40,60;60,80;80,100;100,120];
% $$$ %mychans = [71,73,81,85,95];
% $$$ mychans = 1:2:8;
% $$$ B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
% $$$ figure
% $$$ for c = 1:numel(mychans),
% $$$     for i = flim',
% $$$         subplot2(numel(mychans),size(flim,1),c,find(i(1)==flim(:,1)));
% $$$         tpow = log10(mean(yld(aind,fl>i(1)&fl<i(2),find(chans==mychans(c))),2));
% $$$         %tpow = log10(mean(yld(aind,fh>i(1)&fh<i(2),find(chans==mychans(c))),2));
% $$$         tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
% $$$         AFB = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
% $$$         AFB(B<5)=nan;
% $$$         AFBclims = [prctile(AFB(~isnan(AFB)),5),prctile(AFB(~isnan(AFB)),95)];
% $$$         imagescnan({vedges,sedges,AFB'},AFBclims,[],1,[0,0,0]);
% $$$         axis xy,
% $$$         ticks_lin2log(gca,'x')
% $$$         title(['C: ' num2str(mychans(c)) 'Mean P(' num2str(i(1)) '-' num2str(i(2)) ')'])
% $$$     end
% $$$ end
% $$$ 
% $$$ text(.1,.1,['Mean P(' num2str(flim(1)) '-' num2str(flim(end)) ') given 10Hz Osc. Power dist(SU,HB) VS Vel(SL)'])
% $$$ ylabel('log10(P(10Hz)) of Spine to Head Distance')
% $$$ xlabel('Head Speed (cm/s)')
% $$$ 
% $$$ 
% $$$ A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@nanmean,nan);
% $$$ figure,imagescnan({vedges,sedges,clip(A,0,.015)'},[],[],1,[0,0,0]),axis xy,
% $$$ % $$$ 
% $$$ % $$$ %tbp_phase.resample(yad);
% $$$ % $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_median,nan);
% $$$ % $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_var,nan);
% $$$ % $$$ figure,imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ 
% $$$ % $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4|(fl>12&fl<18),chan),2);
% $$$ % $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4,chan),2);
