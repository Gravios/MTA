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



% Ed01-20140707.cof.all
% Ed01-20140709.cof.all
% Ed01-20140717.cof.all
% Ed05-20140528.cof.all
% Ed05-20140529.ont.all
% Ed10-20140815.cof.all
% Ed10-20140816.cof.all
% Ed10-20140817.cof.gnd


% Find ncp troughs 
% get rhm or body 5-12 Hz phase at ncp troughs
% compute phase difference
% get frequencey

%Trial = MTATrial.validate('Ed01-20140717.cof.all');

configure_default_args();

sessionList = get_session_list_v2('ncp');
Trials = af(@(s) MTATrial.validate(s), sessionList);

dataTrain = [];
dataTarget = [];

tinds = [1,5];

Trial = Trials{tinds(2)};
stc = copy(Trial.stc);
stc = Trial.load('stc','msnn_ppsvd_raux');



sampleRate = 250;
xyz = filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',4,[30],'low');
hdo = transform_origin(Trial,xyz,'hcom','nose',{'head_left','head_right'});
dyaw = copy(xyz);
dyaw.data = circ_dist(circshift(hdo.direction,-1),circshift(hdo.direction,1));
dyaw.filter('ButFilter',4,[4,15],'bandpass');
dpch = copy(xyz);
dpch.data = circ_dist(circshift(hdo.pitch,-1),circshift(hdo.pitch,1));
dpch.filter('ButFilter',4,[4,15],'bandpass');
drll = copy(xyz);
drll.data = circ_dist(circshift(hdo.roll,-1),circshift(hdo.roll,1));
drll.filter('ButFilter',4,[4,15],'bandpass');

% CONDITIONAL VARS
pch = fet_HB_pitchB(Trial,sampleRate); 
vxy = vel(filter(preproc_xyz(Trial,'trb',sampleRate),'ButFilter',4,2.4,'low'),{'hcom'},[1,2]);
rfq = fet_respiration_freq(Trial,sampleRate,'overwrite',false);

% TRAINING VARS
ang = create(MTADang,Trial,xyz);
%rhm = fet_rhm(Trial,sampleRate);
bresp = copy(ang);
bresp.data = nunity(ang(:,'pelvis_root','spine_upper',3));
bresp.filter('ButFilter',4,[4,15],'bandpass');
lresp = copy(ang);
lresp.data = nunity(ang(:,'pelvis_root','head_back',3));
lresp.filter('ButFilter',4,[4,15],'bandpass');

% DEFVAR lowpass filtered nasal pressure sensor
ncpFilt = resample(filter(Trial.load('lfp',Trial.meta.channelGroup.respiration),'ButFilter',4,16,'low'),rhm);
ncpFilt.data = nunity(ncpFilt);

% OUTPUT VARS
% ncpFilt
% head speed
% head pitch
% body pitch




% GET segments ---------------------------------------------------------------------
twin = 64;

gper = resample(cast([stc{'gper-groom'}],'TimeSeries'),xyz);
gperSegs = gper.segs(1:twin/2:size(xyz,1),twin);
gperSegs = all(gperSegs)';

ncpSegs = ncpFilt.segs(1:twin/2:size(xyz,1),twin)';

dyawSegs = dyaw.segs(1:twin/2:size(xyz,1),twin)'; 
dpchSegs = dpch.segs(1:twin/2:size(xyz,1),twin)'; 
drllSegs = drll.segs(1:twin/2:size(xyz,1),twin)'; 

pchSegs = pch.segs(1:twin/2:size(xyz,1),twin); pchSegs = mean(pchSegs(:,:,1))';
vxySegs = vxy.segs(1:twin/2:size(xyz,1),twin); vxySegs = sq(mean(vxySegs))';
rfqSegs = rfq.segs(1:twin/2:size(xyz,1),twin); rfqSegs = sq(mean(rfqSegs))';

%rhmSegs   = clip(  rhm.segs(1:twin/2:size(xyz,1), twin)', -4, 4);
lrespSegs = clip(lresp.segs(1:twin/2:size(xyz,1), twin)', -1, 1);
brespSegs = clip(bresp.segs(1:twin/2:size(xyz,1), twin)', -1, 1);


% PRESELECT good periods -----------------------------------------------------------
% $$$ size(rhmSegs)   
% $$$ size(brespSegs) 
% $$$ size(lrespSegs)
% $$$ size(ncpSegs)
% $$$ size(pchSegs)
% $$$ size(vxySegs)
% $$$ size(gperSegs)

bind = nniz(brespSegs) & nniz(lrespSegs) ...
       & nniz(ncpSegs) & nniz(pchSegs)   & nniz(vxySegs) & gperSegs ...
       & nniz(drllSegs) & nniz(dpchSegs) & nniz(dyawSegs);
ncpSegs(~bind,:) = [];
dyawSegs(~bind,:) = [];
dpchSegs(~bind,:) = [];
drllSegs(~bind,:) = [];
brespSegs(~bind,:) = [];
lrespSegs(~bind,:) = [];
%rhmSegs(~bind,:) = [];
pchSegs(~bind,:) = [];
vxySegs(~bind,:) = [];
rfqSegs(~bind,:) = [];

clear('pchSegsInds','rfqSegsInds','vxySegsInds');
vxySegsInds = discretize(vxySegs,logspace(-3,2,6));
pchSegsInds(:,1) = discretize(pchSegs(:,1),linspace(-2,0.5,6));
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
                                    dyawSegs(ind,:,:),...
                                    dpchSegs(ind,:,:),...
                                    drllSegs(ind,:,:),...
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
net = feedforwardnet([32,32]);
net = train(net,dataTrain',dataTarget');

adyawSegs = dyaw.segs(1:size(xyz,1),twin)'; 
adpchSegs = dpch.segs(1:size(xyz,1),twin)'; 
adrllSegs = drll.segs(1:size(xyz,1),twin)'; 
alrespSegs = clip(lresp.segs(1:size(xyz,1),twin)',-1,1);
abrespSegs = clip(bresp.segs(1:size(xyz,1),twin)',-1,1);
%arhmSegs = clip(rhm.segs(1:size(rhm,1),twin)',-4,4);
apchSegs = pch.segs(1:size(xyz,1),twin);apchSegs = sq(mean(apchSegs(:,:,1)))';
avxySegs = vxy.segs(1:size(xyz,1),twin);avxySegs = sq(mean(avxySegs))';


nsegs = net(cat(2,abrespSegs,alrespSegs,adyawSegs,adpchSegs,adrllSegs,apchSegs,avxySegs)');
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

ncpFet = fet_respiration_reconstruction(Trial,250,2,'overwrite',false);

figure();
hold('on');
plot(ncpFilt.data,'c')
plot(ncpFet.data,'k')
plot(nnsegs(:,1),'m');


sfet = copy(ncpFilt);
sfet.data = [sfet.data,nnsegs,rhm.data];
defspec=struct('nFFT',2^10,'Fs',sfet.sampleRate,...
               'WinLength',2^9,'nOverlap',2^9*0.5,...
               'FreqRange',[1,15]);

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

ind = stc{'w'};
mcoh = sum(abs(ys(ind,:,1,2)))./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,2,2)));
mcohr = sum(abs(ys(ind,:,1,3)))./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,3,3)));

figure();
hold('on');
plot(mcoh);
plot(mcohr);

figure();
hold('on');
plot(mcoh-mcohr);

pchs = resample(copy(pch),ys);
% ncp, rsy, rhm
ncps.pow = sqrt(ys(:,:,1,1).^2);
rsys.pow = sqrt(ys(:,:,2,2).^2);
rhms.pow = sqrt(ys(:,:,3,3).^2);


rsys.edgs = linspace(-6,-1,31);
rsys.cntr = mean([rsys.edgs(1:end-1);rsys.edgs(2:end)]);
rsys.inds = discretize(mean(log10(rsys.pow(:,fs>2&fs<13)),2),rsys.edgs);
gperc = resample(cast([stc{'gper-groom'}],'TimeSeries'),ys);
gperc = resample(cast([stc{'lloc+lpause'}],'TimeSeries'),ys);
%gperc = resample(cast([stc{'hloc+hpause'}],'TimeSeries'),ys);
ccnt = [];
ccoh = [];
c = 2;
for bin = 1:numel(rhms.edgs)-1
    ind = rsys.inds==bin & gperc.data;
    if sum(ind)>10
        ccnt(bin) = sum(ind);
        ccoh(bin,:) = sum(abs(ys(ind,:,1,c)),'omitnan')./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,c,c)),'omitnan');
    else
        ccnt(bin) = 0;
        ccoh(bin,:) = zeros([1,numel(fs)]);;
    end
end


figure,
subplot(4,1,[1,2,3]);
imagesc(rsys.cntr,fs,ccoh');
axis('xy');
colormap('jet');
caxis([0.2,1]);
subplot(4,1,4);
bar(rsys.cntr,ccnt);




rhms.edgs = linspace(-9,-3.75,31);
rhms.cntr = mean([rhms.edgs(1:end-1);rhms.edgs(2:end)]);
rhms.inds = discretize(mean(log10(rhms.pow(:,fs>2&fs<13)),2),rhms.edgs);
%gperc = resample(cast([stc{'gper-groom'}],'TimeSeries'),ys);
gperc = resample(cast([stc{'lloc+lpause'}],'TimeSeries'),ys);
%gperc = resample(cast([stc{'hloc+hpause'}],'TimeSeries'),ys);
ccntr = [];
ccohr = [];
c = 3;
for bin = 1:numel(rhms.edgs)-1
    ind = rhms.inds==bin & gperc.data;
    if sum(ind)>10
        ccntr(bin) = sum(ind);
        ccohr(bin,:) = sum(abs(ys(ind,:,1,c)),'omitnan')./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,c,c)),'omitnan');
    else
        ccntr(bin) = 0;
        ccohr(bin,:) = zeros([1,numel(fs)]);;
    end
end


figure,
subplot(4,1,[1,2,3]);
imagesc(rhms.cntr,fs,ccohr');
axis('xy');
colormap('jet');
caxis([0.2,1]);
subplot(4,1,4);
bar(rhms.cntr,ccntr);



lrfq = resample(copy(rfq),ys);
rfqs.edgs = linspace(0.5,12.5,13);
rfqs.cntr = mean([rfqs.edgs(1:end-1);rfqs.edgs(2:end)]);
rfqs.inds = discretize(lrfq.data,rfqs.edgs);

gperc = resample(cast([stc{'gper-groom'}]+[-3,3],'TimeSeries'),ys);
ccntf = [];
ccohf = [];
c = 2;
for bin = 1:numel(rfqs.edgs)-1
    ind = rfqs.inds==bin&gperc.data;
    ccntf(bin) = sum(ind);
    ccohf(bin,:) = sum(abs(ys(ind,:,1,c)),'omitnan')./sum(sqrt(ys(ind,:,1,1).*ys(ind,:,c,c)),'omitnan');
end

figure,
subplot(4,1,[1,2,3]);
imagesc(rfqs.cntr,fs,ccohf');
axis('xy');
colormap('jet');
caxis([0.2,1]);
subplot(4,1,4);
bar(rfqs.cntr,ccntf);
xlim([0.5,12.5])



sresp = copy(rhm);
sresp.data = ang(:,'pelvis_root','spine_upper',3);
sresp.filter('ButFilter',4,[1,6],'bandpass');


figure();
plot(nunity(MedianFilter(ang(:,'pelvis_root','spine_upper',3)-mean(ang(:,'pelvis_root','spine_upper',3)),600)))
hold('on');
plot(nunity(sresp.data))
plot(nunity(ncpFilt.data))
plot(nunity(ncpSyn.data),'m')

figure,
hold('on');
plot(sqrt(conv(nunity(ncpSyn.data).^2,gausswin(512)./sum(gausswin(512)),'same')));
plot(nunity(sresp.data))
plot(nunity(ncpFilt.data)./10)
Lines([],0.1,'k');

ncpSynPow = copy(ncpSyn);
ncpSynPow.data = sqrt(conv(nunity(ncpSyn.data).^2,gausswin(512)./sum(gausswin(512)),'same'));

srespPeriods = ThreshCross(-ncpSynPow.data,-0.15,250);
ncpSynFinal = copy(ncpSyn);
for s = 1:size(srespPeriods,1)
    ncpSynFinal.data(srespPeriods(s,1):srespPeriods(s,2)) =  ...
        sresp.data(srespPeriods(s,1):srespPeriods(s,2));
end


figure,
hold('on');
plot(nunity(ncpFilt.data))
plot(nunity(ncpSyn.data),'k')
plot(nunity(ncpSynFinal.data),'m')



sfet = copy(ncpFilt);
sfet.data = [sfet.data,nunity(ncpSynFinal.data),rhm.data];
defspec=struct('nFFT',2^10,'Fs',sfet.sampleRate,...
               'WinLength',2^9,'nOverlap',2^9*.5,...
               'FreqRange',[1,15]);

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








sfet = copy(ncpFilt);
sfet.data = [sfet.data,nunity(bresp.data),nunity(rhm(:,1))];
defspec=struct('nFFT',2^10,'Fs',sfet.sampleRate,...
               'WinLength',2^9,'nOverlap',2^9*.5,...
               'FreqRange',[1,15]);
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

figure,
subplot(211);
imagesc(ts,fs,circ_dist(angle(ys(:,:,1,2))',angle(ys(:,:,1,3))'));
axis('xy');
colormap('hsv');
caxis([-pi,pi]);
subplot(212);
plot(ts,pchs(:,1));
hold('on');
plot([1:size(pch,1)]./pch.sampleRate,pch(:,1));
plot([1:size(rhm,1)]./rhm.sampleRate,rhm(:,1)*5);
plot([1:size(ncpFilt,1)]./ncpFilt.sampleRate,nunity(ncpFilt(:,1))./4);
linkx();



figure,
dpch = circ_dist(circshift(ang(:,'hcom','nose',2),-1),circshift(ang(:,'hcom','nose',2),1));
dyaw = circ_dist(circshift(ang(:,'hcom','nose',1),-1),circshift(ang(:,'hcom','nose',1),1));
drol = circ_dist(circshift(hobj.roll,-1),circshift(hobj.roll,1));

figure,plot(dyaw(1000:1100),dpch(1000:1100))

figure,plot(sqrt(sum([dyaw,dpch].^2,2)));


dphi = copy(ang);
dphi.data = [dyaw,dpch,drol];
dphi.filter('ButFilter',4,[5,12],'bandpass');

hobj = Trial.transform_origin(xyz,'hcom','nose',{'head_left','head_right'});


figure();
hold('on');
plot(nunity(dphi(:,1))+3,'b');
plot(nunity(dphi(:,2))+3,'r');
plot(nunity(dphi(:,3))+3,'c');
plot(nunity(ncpFilt.data),'k');
plot(hobj.roll*3-5);
plot(ang(:,'hcom','nose',2)-5);
Lines([],-5,'k');

figure();
plot(sqrt(sum(dphi.data.^2,2)));
hold('on');
plot(nunity(ncpFilt.data)/30);

