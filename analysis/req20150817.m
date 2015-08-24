
stc = Trial.stc.copy;


wf = mean(log10(vmv(:,1:2)),2);
Mwf = MTADxyz('data',wf,'sampleRate',trajSampleRate);

eds = linspace(-6,6,200);
hfig = figure;hold on;
ind = Trial.stc{'a'};
ha = bar(eds,histc(Mwf(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'w'};
hs = bar(eds,histc(Mwf(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;
Lines(wft,[],'k');
stc.load(Trial,'auto_wbhr');
ind = stc{'w'};
hs = bar(eds,histc(Mwf(ind),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .8;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;
Lines(wft,[],'k');

title('Step one: Threshold Coherent Movement');
xlabel('mean(mean(distance)*mean(variance))')
ylabel('Count');


wdist_score = [];
for i= stc{'w'}.data',
    wdist_score(end+1) = sum(sqrt(sum(sq(diff(fxyz(i(1):i(2),:))).^2,2)));
end
wdist_score = log10(wdist_score);



Maf = MTADxyz('data',circ_mean(atrajMean,[],2).*mean(atrajVarD,2),'sampleRate',trajSampleRate);
Maf.filter('ButFilter',3,3,'low');

figure,plot(Maf.data)
figure,plot(diff(Maf.data))
Lines(Trial.stc{'n',Maf.sampleRate}(:)-4,[],'g');

Maf.data = log10(abs(Maf.data));



eds = linspace(-12,1,200);
hfig = figure;hold on;
ind = Trial.stc{'a-n'};
ha = bar(eds,histc(Maf(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'n'};
hs = bar(eds,histc(Maf(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;
Lines(wft,[],'k');




figure,plot(af)
Lines(Trial.stc{'n',trajSampleRate}(:),[],'g');
Lines(sturn(:),[],'m');
Lines([],-1.6,'k');


tangc = [];
sturn = MTADepoch(afp-(stc{'r',trajSampleRate}+[-.25,.25])).data;
for a = MTADepoch(afp-(stc{'r',trajSampleRate}+[-.25,.25])).data'
    tangc(end+1) = circ_dist(ang(a(1),1,4,1),ang(a(2),1,4,1));
end
figure,hist(log10(abs(tangc)),100),

tangc = [];
for a = stc{'n',trajSampleRate}.data',
    tangc(end+1) = circ_dist(ang(a(1),1,4,1),ang(a(2),1,4,1));
end
figure,hist(log10(abs(tangc)),100),


Mtn = MTADxyz('data',log10(abs(circ_dist(circshift(ang(:,1,4,1),-4),circshift(ang(:,1,4,1),4)))),...
              'sampleRate',trajSampleRate);
figure,plot(mtn.data)
Lines(Trial.stc{'n',trajSampleRate}(:),[],'g');

eds = linspace(-6,.5,200);
hfig = figure;hold on;
ind = Trial.stc{'a-n-w-r'};
ha = bar(eds,histc(Mtn(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'n'};
hs = bar(eds,histc(Mtn(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;


sfet = [circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
        circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
        circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
        circ_dist(ang(:,5,7,1),ang(:,4,5,1))];


Msf = MTADxyz('data',abs(sum(sfet,2)),'sampleRate',Trial.xyz.sampleRate);

eds = linspace(0,pi*2,200);
hfig = figure;hold on;
ind = Trial.stc{'a-m'};
ha = bar(eds,histc(Msf(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'m'};
hs = bar(eds,histc(Msf(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;

Mang = MTADxyz('data',ang(:,1,4,3).*cos(ang(:,1,4,2)),'sampleRate',Trial.xyz.sampleRate);

ads = linspace(40,160,100);
eds = linspace(0,pi*2,100);

figure,
ind = stc{'a-m-r'}
hist2([Mang(ind),Msf(ind)],ads,eds);caxis([0,100])
figure,
ind = stc{'m'}
hist2([Mang(ind),Msf(ind)],ads,eds);caxis([0,100])
figure,
ind = stc{'r'}+[-.5,0];
hist2([Mang(ind),Msf(ind)],ads,eds);caxis([0,100])



figure,
ind = stc{'a'}
hist2([Mang(ind),ang(ind,3,4,2)],ads,linspace(-1,1.5,100));caxis([0,100])
figure,
ind = stc{'r'}
hist2([Mang(ind),ang(ind,3,4,2)],ads,linspace(-1,1.5,100));caxis([0,100])



Mang = MTADxyz('data',ang(:,3,4,2),'sampleRate',Trial.xyz.sampleRate);
figure,plot(nanmean(diff(Mang.segs(stc{'r'}(:,1)-60,120,nan)),2))
figure,plot(nanmean(diff(Mang.segs(stc{'r'}(:,2)-60,120,nan)),2))

figure,plot(nanmean(Mang.segs(stc{'r'}(:,1)-60,120,nan),2))

figure,plot(mean(bsxfun(@times,Mang.segs(1:Mang.size(1),120,nan), ...
                               nanmean(Mang.segs(stc{'r'}(:,1)-60,120,nan),2))))

Lines(stc{'r'}(:,1),[],'g');

hold on,plot(mean(bsxfun(@times,Mang.segs(1:Mang.size(1),120,nan), ...
                               nanmean(Mang.segs(stc{'r'}(:,2)-60,120,nan),2))))


rd = (mean(bsxfun(@times,Mang.segs(1:Mang.size(1),120,nan), ...
                   nanmean(Mang.segs(stc{'r'}(:,1)-60,120,nan),2)))-...
     mean(bsxfun(@times,Mang.segs(1:Mang.size(1),90,nan), ...
                   nanmean(Mang.segs(stc{'r'}(:,2)-60,120,nan),2))))';
figure,plot(rd)
Lines(stc{'r'}(:,1)-60,[],'r');
Lines(stc{'r'}(:,2)-60,[],'k');


rd = (mean(bsxfun(@times,Mang.segs(1:Mang.size(1),90,nan), ...
                   nanmean(Mang.segs(stc{'r'}(:,1)-45,90,nan),2)))-...
     mean(bsxfun(@times,Mang.segs(1:Mang.size(1),90,nan), ...
                   nanmean(Mang.segs(stc{'r'}(:,2)-45,90,nan),2))))';
figure,plot(rd)
Lines(stc{'r'}(:,1)-45,[],'r');
Lines(stc{'r'}(:,2)-45,[],'k');


[mid,miv] = LocalMinimaN(-abs(rd),-.0002,40);


[mir,mirv] = LocalMinimaN(-abs(rd),-.2,30);

figure,hist(log10(-miv),100)
Lines(-.699,[],'k');


rper = stc{'r'}(:)-45;
mrd = repmat(mir(:,1),[1,size(rper)])-repmat(rper',[length(mir),1]);


mrdl = any(abs(mrd)<60,2);
for r = 1:numel(rper),
   mrd
end
   
eds = linspace(-4,0,200);
figure,hold on
ha = bar(eds,histc(log10(-miv(~mrdl)),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
hs = bar(eds,histc(log10(-miv(mrdl)),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;


fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
% $$$ fvel.data(fvel.data<0)=.1;


Mvel = MTADxyz('data',fvel(:,1),'sampleRate',Trial.xyz.sampleRate);

s = 'w+n';

win = 90;
hwin = round(win/2);
rd = (mean(bsxfun(@times,Mvel.segs(1:Mvel.size(1),win,0), ...
                   nanmean(Mvel.segs(stc{s}(:,1)-hwin,win,0),2)))-...
     mean(bsxfun(@times,Mvel.segs(1:Mvel.size(1),win,0), ...
                   nanmean(Mvel.segs(stc{s}(:,2)-hwin,win,0),2))))';
figure,plot(rd)
Lines([stc{s}(:,1)]-hwin,[],'r');
Lines([stc{s}(:,2)]-hwin,[],'k');

figure,plot(abs(rd).*circshift(Mvel.data,-hwin)),
Lines(stc{s}(:,1)-60,[],'r');
Lines(stc{s}(:,2)-60,[],'k');



Mas = MTADxyz('data',circ_dist(ang(:,1,3,1),ang(:,2,4,1)),'sampleRate',Trial.xyz.sampleRate);



%% ANT Vel ang dist between high and low walk
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
rpow = MTADxyz('data',rhm(:,fs<15&fs>6),'sampleRate',rhm.sampleRate);
rpow.resample(xyz);


stc = Trial.stc.copy;
stc.load(Trial,'auto_wbhr');

edv = linspace(-7,-2,100);
eda = linspace(-1,1.5,100);


figure,
subplot(131)
ind = stc{'w'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);
subplot(132)
ind = stc{'h'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);
subplot(133)
ind = stc{'l'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);

figure,
subplot(131)
ind = stc{'w'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);
subplot(132)
ind = stc{'c'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);
subplot(133)
ind = stc{'p'};
hist2([rpow(ind,1),ang(ind,5,7,2)],edv,eda)
ylabel('ang(rad)'),xlabel('speed log10(cm/s)')
title(ind.label);
%% END ANT Vel ang dist between high and low walk
