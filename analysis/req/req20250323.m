

MjgER2016_load_data();

UnitsInts = cf(@(T)   T.spk.get_unit_set( T, 'interneurons'), Trials); 

pyr = units{20};
ints = UnitsInts{20};
bhvstates = 'theta-groom-sit';
samplerate = 250;

stc = Trial.stc.copy();
xyz = preproc_xyz(Trial);
xyz.resample(samplerate);
spkp = Trial.load('spk',samplerate,bhvstates,pyr);
phz = load_theta_phase(Trial,samplerate);


sts = [stc{'loc',samplerate}];
figure();
plot( ...
    xyz(sts,'hcom',1),...
    xyz(sts,'hcom',2),....
    '.k');

res = spkp(21);
res = SelectPeriods(res, sts.data, 'd', 1,0);
hold('on');
plot( ...
    xyz(res,'hcom',1),...
    xyz(res,'hcom',2),....
    '.r','MarkerSize',10);

rpos = sq(xyz(sts,'hcom',[1:2]));
spos = xyz(res,'hcom',[1:2]);
nbin = 30;

ind = nniz(rpos);
o[rocc] = hist2(rpos(ind,:), linspace(-500,500,nbin), linspace(-500,500,nbin));
figure,
imagesc(rocc')
axis('xy');

ind = nniz(spos);
[socc] = hist2(spos(ind,:), linspace(-500,500,nbin), linspace(-500,500,nbin));
figure,
imagesc(socc')
axis('xy');

figure,
imagesc(socc'./rocc')
axis('xy');
colorbar()
colormap('jet')



ndims = 2;
sind = cell(1,ndims);
for i = 1:ndims,
    sind{i} = linspace(-500,500,nbin);
end
[sind{:}] = ndgrid(sind{:});
for i = 1:ndims,
    sind{i} = sind{i}.^2/50^2/2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));
figure,imagesc(Smoother)

srocc   = convn(rocc, Smoother,'same');
ssocc = convn(socc,Smoother,'same');


figure,
imagesc(ssocc'./(srocc/samplerate)')
axis('xy');
colorbar()
colormap('jet')






hpos = bsxfun(@minus,rpos,[-120,155]);

[hTH, hR] = cart2pol(hpos(:,1), hpos(:,2));

rhpos = sq(xyz(sts,'head_front',[1,2])-xyz(sts,'head_back',[1,2]));
[rTH, rR] = cart2pol(rhpos(:,1), rhpos(:,2));

sako = circ_dist( hTH,rTH);

hfdir = sako;
hfdir(WithinRanges(sako,[-pi/2,pi/2])) = 1;
hfdir(sako<-pi/2 | sako>pi/2) = -1;

lphz = phz(sts);

res = spkp(21);
res = SelectPeriods(res, sts.data, 'd', 1,1);

figure()
plot( ...
    hR(res).*hfdir(res),...
    lphz(res),...
    '.');







lfp = Trial.load('lfp',65:96);
figure();
plot(bsxfun(@plus, lfp.data(10000:40000,:), fliplr([1:32])*2000),'k')

rfp = Trial.load('lfp',[57,64]);



figure,
hold('on');
plot(lfp.data(1:5e4,3),'k')
plot(diff(rfp.data(1:5e4,:),1,2)+10000,'g');
plot(lfp.data(1:5e4,18)-10000,'r')


figure,
hold('on');
plot(nunity(lfp.data(1:5e4,3)),'k')
plot(nunity(lfp.data(1:5e4,18))-2,'r')
plot(nunity(diff(rfp.data(1:5e4,:),1,2))+2,'g');

rfp.resample(samplerate);
slfp = lfp.copy();
slfp.data = slfp(:,3);

rlfp = lfp.copy();
rlfp.resample(samplerate);
slfp = rlfp.copy();
slfp.data = [slfp(:,[3,9,18]),diff(rfp.data,1,2)];


defspec = struct('nFFT',2^9,'Fs',slfp.sampleRate,...
                 'WinLength',2^8,'nOverlap',2^8*.875,...
           'FreqRange',[1,40]);
[ys,fs,ts] = fet_spec(Trial,slfp,'mtcsdglong',false,slfp.sampleRate,defspec,[],true);

fet = fet_mis(Trial);

ny = 7;
figure,
for chan = 1:4
    subplot(ny,1,chan);
    imagesc(ts,fs,log10(ys(:,:,chan,chan))');
    colormap('jet');
    axis('xy');
end
chan = chan+1;
subplot(ny,1,chan); chan = chan+1;
plot([1:size(xyz,1)]./xyz.sampleRate,xyz(:,:,3))
subplot(ny,1,chan); chan = chan+1;
imagesc([1:size(fet,1)]./fet.sampleRate, 1:size(fet,2), nunity(fet.data)');
subplot(ny,1,chan);
plotSTC(stc,1);
linkx();


figure,
imagesc(ts,fs,angle(ys(:,:,1,4))');
axis('xy');
colormap('hsv');

thpow = log10(mean(ys(:,WithinRanges(fs,[6,10]),chan,chan),2)) ...
             - log10(mean(ys(:,fs<3|WithinRanges(fs,[12,15]),chan,chan),2));

CheckEegStates(fullfile(Trial.spath,Trial.name),[],[],[],67,[],'display',0);


tp = log10(mean(ys(:,WithinRanges(fs,[6,10]),chan,chan),2));
dp = log10(mean(ys(:,fs<3|WithinRanges(fs,[12,15]),chan,chan),2));

ind = nniz([tp,dp]);
figure,
hist2( [tp(ind), dp(ind)], 50,50);



res = spkp(20);

dzh = circshift(xyz(:,'head_back',3),-1)-circshift(xyz(:,'head_back',3),1);

adzh = abs(dzh);

dzh_sign = sign(dzh);

ladzh = clip(log10(adzh),0,inf);

sladzh = dzh_sign.*ladzh;


rper = ThreshCross(xyz(:,'head_back',3),150,20);

zh = xyz(:,'head_back',3);
figure,
plot(zh);
Lines(rper(:,1),[],'r');
Lines(rper(:,2),[],'k');

zh_segs = GetSegs(zh, rper(:,1)-120, 240,0);

figure,
hold('on')
plot(mean(zh_segs'))
plot(mean(zh_segs')-std(zh_segs'))
plot(mean(zh_segs')+std(zh_segs'))


xy_dhs = xyz(:,'head_back',[1,2])- xyz(:,'spine_upper',[1,2]);

[TH, R] = cart2pol(xy_dhs(:,1), xy_dhs(:,2));


figure
plot(R)
% mean height of rear



figure
hist2([xyz(:,'head_back',3), ...
       sladzh],50,50)

samplerate = 250;
bhvstates = 'theta-groom-sit-rear';
spkmode = 'all';

Trial = Trials{20};
interneurons = unitsInts{20};


% X, Y, HDir, Theta

spki =                          ...
    Trial.load('spk',           ...
        samplerate,             ...
        bhvstates,              ...
        interneurons,           ...
        spkmode);

