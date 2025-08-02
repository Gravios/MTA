% req20220217
%


% Data preparation
% COMPUTE the intantaneous frequency of the respiration signal
% SAMPLE random training segments to generate a uniform joint distribution in ifreq and head pitch domains
% 

% Checks
% Mean coherence of ncp with rhm and ncpSyn where { log10(ncp)>10 & goodPeriods-groom }
% Mean coherence of ncpSyn with ncp where { log10(ncpSyn)>10 & goodPeriods-groom }
%  

configure_default_args();

sessionList = get_session_list('hand_labeled');
Trials = af(@(s) MTATrial.validate(s), sessionList);


dataTrain = [];
dataTarget = [];

tinds = [3,6];

Trial = Trials{tinds(2)};
stc = copy(Trial.stc);

sampleRate = 250;

% CONDITIONAL VARS
pch = fet_HB_pitchB(Trial,sampleRate); 
vxy = vel(filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',4,2.4,'low'),{'hcom'},[1,2]);
rfq = fet_respiration_freq(Trial,sampleRate,'overwrite',true);

% TRAINING VARS
ang = create(MTADang,Trial,preproc_xyz(Trial,'trb',sampleRate));
rhm = fet_rhm(Trial,sampleRate);
bresp = copy(ang);
bresp.data = nunity(ang(:,'pelvis_root','spine_upper',3));
bresp.filter('ButFilter',4,[0.8,15],'bandpass');
lresp = copy(ang);
lresp.data = nunity(ang(:,'pelvis_root','head_back',3));
lresp.filter('ButFilter',4,[0.8,15],'bandpass');

% DEFVAR lowpass filtered nasal pressure sensor
ncpFilt = resample(filter(Trial.load('lfp',2),'ButFilter',4,16,'low'),rhm);
ncpFilt.data = nunity(ncpFilt);

% OUTPUT VARS
% ncpFilt

% head speed
% head pitch
% body pitch




% GET segments ---------------------------------------------------------------------
twin = 250;

gper = resample(cast([stc{'gper'}],'TimeSeries'),rhm);
gperSegs = gper.segs(1:twin/2:size(rhm,1),twin);
gperSegs = all(gperSegs)';

ncpSegs = ncpFilt.segs(1:twin/2:size(rhm,1),twin)';

pchSegs = pch.segs(1:twin/2:size(rhm,1),twin); pchSegs = pchSegs(:,:,1)';
pchSegsMean = pch.segs(1:twin/2:size(rhm,1),twin); pchSegsMean = sq(mean(pchSegsMean(:,:,1)))';
vxySegs = vxy.segs(1:twin/2:size(rhm,1),twin); vxySegs = sq(mean(vxySegs))';
rfqSegs = rfq.segs(1:twin/2:size(rhm,1),twin); rfqSegs = sq(mean(rfqSegs))';

rhmSegs   = clip(  rhm.segs(1:twin/2:size(rhm,1), twin)', -4, 4);
lrespSegs = clip(lresp.segs(1:twin/2:size(rhm,1), twin)', -1, 1);
brespSegs = clip(bresp.segs(1:twin/2:size(rhm,1), twin)', -1, 1);


% PRESELECT good periods -----------------------------------------------------------
% $$$ size(rhmSegs)   
% $$$ size(brespSegs) 
% $$$ size(lrespSegs)
% $$$ size(ncpSegs)
% $$$ size(pchSegs)
% $$$ size(vxySegs)
% $$$ size(gperSegs)

bind = nniz(rhmSegs)   & nniz(brespSegs) & nniz(lrespSegs) ...
       & nniz(ncpSegs) & nniz(pchSegs)   & nniz(vxySegs) & gperSegs;
ncpSegs(~bind,:) = [];
brespSegs(~bind,:) = [];
lrespSegs(~bind,:) = [];
rhmSegs(~bind,:) = [];
pchSegsMean(~bind,:) = [];
pchSegs(~bind,:) = [];
vxySegs(~bind,:) = [];
rfqSegs(~bind,:) = [];

clear('pchSegsInds','rfqSegsInds','vxySegsInds');
vxySegsInds = discretize(vxySegs,logspace(-3,2,6));
pchSegsInds(:,1) = discretize(pchSegsMean(:,1),linspace(-2,0.5,6));
rfqSegsInds(:,1) = discretize(rfqSegs(:,1),linspace(0,15,6));


dataTrain = [];
dataTarget = [];


sampleCount = 100;
for rf = 1:5
    for v = 1:5
        for ph = 1:5
            ind = find(   v == vxySegsInds ...
                          & ph == pchSegsInds(:,1) ...
                          & rf == rfqSegsInds(:,1));
            if numel(ind)>sampleCount
                ind = ind(randsample(numel(ind),sampleCount,true));
                dataTrain = cat(1,...
                                dataTrain,...
                                cat(2,...
                                    brespSegs(ind,:,:),...
                                    lrespSegs(ind,:,:),...
                                    rhmSegs(ind,:,:),...
                                    pchSegs(ind,:,:),...
                                    vxySegs(ind,:,:)) ...
                                );
                dataTarget = cat(1,dataTarget,ncpSegs(ind,:));
            end
        end
    end
end

figure
subplot(131);
imagesc(nunity(dataTarget(:,:)));caxis([-3,3])
subplot(132);
imagesc(dataTrain(:,101:150));caxis([-0.1,0.1])
subplot(133);
imagesc(dataTrain(:,1:50));caxis([-0.1,0.1])

%clear('net');
net2 = feedforwardnet([100,100]);
net2 = train(net2,dataTrain',dataTarget');


alrespSegs = clip(lresp.segs(1:size(rhm,1),twin)',-1,1);
abrespSegs = clip(bresp.segs(1:size(rhm,1),twin)',-1,1);
arhmSegs = clip(rhm.segs(1:size(rhm,1),twin)',-4,4);
apchSegs = pch.segs(1:size(rhm,1),twin);apchSegs = sq(mean(apchSegs));
avxySegs = vxy.segs(1:size(rhm,1),twin);avxySegs = sq(mean(avxySegs))';


nsegs = net(cat(2,abrespSegs,alrespSegs,arhmSegs,apchSegs,avxySegs)');
nnsegs = (1./twin).*nsegs(1,:)';
for s = 2:size(nsegs,1)
    nnsegs = nnsegs+(1./twin).*circshift(nsegs(s,:)',s);
end    
ncpSyn = copy(ncpFilt);
ncpSyn.data = nnsegs;



figure,hold('on');
% $$$ plot(nsegs(1,:)+8000,'b')
% $$$ plot(circshift(nsegs(20,:)',20)+8000,'c')
plot(nunity(nnsegs(:,1)),'r')
plot(nunity(ncpFilt(:,1)),'b');
plot(nunity(rhm.data)-7,'m')
plot(nunity(lresp.data)-10,'m')
plot(nunity(bresp.data)-14,'m')



sfet = copy(ncpFilt);
sfet.data = [sfet.data,nnsegs,rhm.data];
defspec=struct('nFFT',2^9,'Fs',sfet.sampleRate,...
               'WinLength',2^8,'nOverlap',2^8*.875,...
               'FreqRange',[1,20]);

[ys,fs,ts] = fet_spec(Trial,sfet,'mtcsdglong',false,[],defspec,[],true);

coh = abs(ys(:,:,1,2))./sqrt(ys(:,:,1,1).*ys(:,:,2,2));
cohr = abs(ys(:,:,1,3))./sqrt(ys(:,:,1,1).*ys(:,:,3,3));

figure,
subplot(411)
imagesc(ts,fs,coh')
axis('xy');
colormap(gca(),'jet');
subplot(412)
imagesc(ts,fs,angle(ys(:,:,1,2))');
axis('xy');
colormap(gca(),'hsv');
subplot(413)
imagesc(ts,fs,cohr')
axis('xy');
colormap(gca(),'jet');
subplot(414)
imagesc(ts,fs,angle(ys(:,:,1,3))');
axis('xy');
colormap(gca(),'hsv');
linkx()

% What do I do with this?
% 

ind = stc{'s'};
mcoh = sum(abs(ys(ind,:,1,2)))./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,2,2)));
mcohr = sum(abs(ys(ind,:,1,3)))./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,3,3)));

figure();
hold('on');
plot(mcoh);
plot(mcohr);

