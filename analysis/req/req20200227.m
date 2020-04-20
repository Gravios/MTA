%req20200227


MjgER2016_load_data();

sampleRate = 250;
t = 20;
Trial = Trials{20};
%tclus = units{20};

% LOAD subject position object
xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);

% COMPUTE polar coordinates of horizontal position
% $$$ dmcd = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
% $$$ dmca = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
% $$$                  atan2(diff(xyz(:,{'hcom','nose'},2),1,2),diff(xyz(:,{'hcom','nose'},1),1,2)));

% LOAD lfp
set(Trial.lfp,'filename',[Trial.name,'.lfp']);
lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
% COMPUTE theta LFP phase
phz = lfp.phase([6,12]);
set(phz,'data',unwrap(phz.data));
phz.resample(xyz);
phz = mod(phz.data+pi,2*pi)-pi;
clear('lfp');

if ismember(t,[3,4,5]),
    phz = phz + 2 * 0.78;
elseif ismember(t,[17:23]),
    phz = phz + 0.78;
end
phz(phz<0) = phz(phz<0) + 2*pi;

states = {'theta','rear','groom','sit','turn','walk','pause'};


stc  = Trial.load('stc','msnn_ppsvd_raux');

stcm =  stc2mat(stc,xyz,states);

spk = Trial.spk.create(Trial,sampleRate,'',[],'');


figure();
for u = 1:181,
res = spk(u);
thetaPhzAwk = phz(res(logical(stcm(res,1))&any(logical(stcm(res,[5,6,7])),2)));
thetaPhzRem = phz(res(logical(stcm(res,1))&logical(stcm(res,4))));
clf();
subplot(121);
rose(thetaPhzAwk)
subplot(122);
rose(thetaPhzRem)
title(num2str(u));
waitforbuttonpress();
end


rd = load(Trial,'lfp',76);
lm = load(Trial,'lfp',82);

defspec = struct('nFFT',2^7,'Fs',rd.sampleRate,...
                            'WinLength',2^6,'nOverlap',2^6*.875,...
                            'FreqRange',[40,100]);
   
[yrd,frd,trd,prd] = fet_spec(Trial,rd,'mtchglong',false,[],defspec);
[ylm,flm,tlm,plm] = fet_spec(Trial,lm,'mtchglong',false,[],defspec);


yrdm = copy(yrd);
yrdm.data = max(yrdm.data,[],2);

ylmm = copy(ylm);
ylmm.data = max(ylmm.data,[],2);

resample(yrdm,xyz);
resample(ylmm,xyz);



mylmm = median(log10(abs(ylmm.data)));
myrdm = median(log10(abs(yrdm.data)));

figure();
for u = 1:181,
res = spk(u);
thetaAwk = res(logical(stcm(res,1))&any(logical(stcm(res,[5,6,7])),2));
thetaRem = res(logical(stcm(res,1))&logical(stcm(res,4)));
clf();

if numel(thetaAwk)<10|numel(thetaRem)<10,
    continue;
end

subplot(321);
rose(phz(thetaAwk))
subplot(323); hold('on');
bar(linspace(3.5,6.5,30),histc(log10(abs(ylmm(thetaAwk))),linspace(3.5,6.5,30)),'histc');
Lines(mylmm,[],'r');
Lines(median(log10(abs(ylmm(thetaAwk)))),[],'g');

subplot(325); hold('on');
bar(linspace(3.5,6.5,30),histc(log10(abs(ylmm(thetaRem))),linspace(3.5,6.5,30)),'histc');
Lines(mylmm,[],'r');
Lines(median(log10(abs(ylmm(thetaRem)))),[],'g');

subplot(322);
rose(phz(thetaRem))
title(num2str(u));
subplot(324); hold('on');
bar(linspace(3.5,6.5,30),histc(log10(abs(yrdm(thetaAwk))),linspace(3.5,6.5,30)),'histc');
Lines(myrdm,[],'r');
Lines(median(log10(abs(yrdm(thetaAwk)))),[],'g');
subplot(326); hold('on');
bar(linspace(3.5,6.5,30),histc(log10(abs(yrdm(thetaRem))),linspace(3.5,6.5,30)),'histc');
Lines(myrdm,[],'r');
Lines(median(log10(abs(yrdm(thetaRem)))),[],'g');

waitforbuttonpress();
end

fb = cwtfilterbank('SignalLength',size(rd,1),'SamplingFrequency',rd.sampleRate,...
                   'FrequencyLimits',[40 120]);
[rdw,f] = cwt(WhitenSignal(rd.data),'FilterBank',fb);
[lmw,f] = cwt(WhitenSignal(lm.data),'FilterBank',fb);
rdw = rdw';
lmw = lmw';

figure,
imagesc(1:rd.size(1),f,log10(abs(rdw)))



rdg = copy(rd);
rdg.data = filter(copy(rd),'ButFilter',4,[40,50],'bandpass').data;
for b = [45,50,55,60,65,70,75,80,85, 90, 95,100,105,110,115;...
         55,60,65,70,75,80,85,90,95,100,105,110,115,120,125],
    rdg.data = cat(2,rdg.data,filter(copy(rd),'ButFilter',4,[b(1),b(2)],'bandpass').data);
end

lmg = copy(lm);
lmg.data = filter(copy(lm),'ButFilter',4,[40,50],'bandpass').data;
for b = [45,50,55,60,65,70,75,80,85, 90, 95,100,105,110,115;...
         55,60,65,70,75,80,85,90,95,100,105,110,115,120,125],
    lmg.data = cat(2,lmg.data,filter(copy(lm),'ButFilter',4,[b(1),b(2)],'bandpass').data);
end

lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
% COMPUTE theta LFP phase
phz = lfp.phase([6,12]);


rdg.data = Shilbert(rdg.data);
lmg.data = Shilbert(lmg.data);



spk = Trial.spk.create(Trial,rd.sampleRate,'',[],'');

stc  = Trial.load('stc','msnn_ppsvd_raux');
stcm =  stc2mat(stc,rd,states);
cmap = cool(16);

figure();

out = [];
for u = 1:181,
    res = spk(u);
    thetaAwk = res(logical(stcm(res,1))&any(logical(stcm(res,[5,6,7])),2));
    thetaRem = res(logical(stcm(res,1))&logical(stcm(res,4)));
    if numel(thetaAwk)<10|numel(thetaRem)<10,
        continue;
    end
    out(u,:) = [circ_r(angle(lmw(thetaAwk,5)))./circ_r(angle(rdw(thetaAwk,14))),...
                circ_r(angle(lmw(thetaRem,5)))./circ_r(angle(rdw(thetaRem,14)))];
end

figure,plot(out(units{t},1),out(units{t},2),'.');

    subplot(321);
    hold('on');
    out = [];
% $$$     for b = 1:16
% $$$         out(:,b) = log10(accumarray(discretize(phz(thetaAwk),linspace(-pi,pi,9)),angle(lmg(thetaAwk,b)),[8,1],@circ_r)./accumarray(discretize(phz(thetaAwk),linspace(-pi,pi,9)),angle(rdg(thetaAwk,b)),[8,1],@circ_r));
% $$$     end

    (accumarray(discretize(phz(thetaAwk),linspace(-pi,pi,9)),angle(lmg(thetaAwk,7)),[8,1],@circ_r)./accumarray(discretize(phz(thetaAwk),linspace(-pi,pi,9)),angle(rdg(thetaAwk,2)),[8,1],@circ_r))
    (accumarray(discretize(phz(thetaRem),linspace(-pi,pi,9)),angle(lmg(thetaRem,7)),[8,1],@circ_r)./accumarray(discretize(phz(thetaRem),linspace(-pi,pi,9)),angle(rdg(thetaRem,2)),[8,1],@circ_r))
    
    imagesc(out([5,6,7,8,1,2,3,4],:)');
    caxis([-1.5,1.5]);
    axis('tight');
    colormap('jet');
    subplot(322);
    hold('on');
    out = [];
    for b = 1:16
        out(:,b) = log10(accumarray(discretize(phz(thetaRem),linspace(-pi,pi,9)),angle(lmg(thetaRem,b)),[8,1],@circ_r)./accumarray(discretize(phz(thetaRem),linspace(-pi,pi,9)),angle(rdg(thetaRem,b)),[8,1],@circ_r));
    end
    imagesc(out([5,6,7,8,1,2,3,4],:)');
    caxis([-1.5,1.5]);
    colormap('jet');    
    axis('tight');    
    title(num2str(u));
% $$$     plot(circ_mean(angle(rdg(thetaAwk,:))),fbins);
% $$$     subplot(323);
% $$$     plot(circ_mean(angle(rdg(thetaRem,:))),fbins);
% $$$     subplot(322);
% $$$     plot(circ_mean(angle(lmg(thetaAwk,:))),fbins);
% $$$     subplot(324);
% $$$     plot(circ_mean(angle(lmg(thetaRem,:))),fbins);
% $$$ 
% $$$     ForAllSubplots('xlim([-pi,pi]);');
    waitforbuttonpress();
end
