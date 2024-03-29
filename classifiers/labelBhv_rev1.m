

Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120317');
%hostPath = '/gpfs01/sirota/homes/gravio/figures/labelBhv_rev1/';
%hostPath = '/storage/gravio/figures/labelBhv_rev1/';


% XYZ Positions of Markers
xyz = Trial.load('xyz');

% COM Body Center of Mass
rbb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
xyz.addMarker('hcom',[.7,0,.7],{{'head_back','head_front',[0,0,255]}},xyz.com(rbb));

rba = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'});
xyz.addMarker('acom',[.7,0,.7],{{'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front',[0,0,255]}},xyz.com(rba));


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');


xyz.addMarker('fhcom',[.7,0,.7],{{'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front',[0,0,255]}},fxyz(:,10,:));

fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');


% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

%figure,plot(circ_dist(circshift(ang(:,10,12,1),-2),circshift(ang(:,10,12,1),2)))

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% FVEL Filtered marker speeds in XY plane
fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0)=.1;


fac = fvel.copy;
fac.data = [diff(fac.data);zeros([1,fvel.size(2)])];


fvel.data = log10(fvel.data);

% UVEL 
uvel = xyz.vel([],[3]);
uvel.filter('ButFilter',3,2.5,'low');
uvel.data(uvel.data<0)=.1;

uac = uvel.copy;
uac.data = [diff(uac.data);zeros([1,uvel.size(2)])];

uvel.data = log10(uvel.data);

%% End Var Setup




%% Definition of Rearing
% Rat stands on its hind legs and elevates the head about 10 cm.
% 
% Features
%
%     1.  Pitch of the line created by the end junctions of the thoracie
%         vertebrae (upper spine)
%
%     2.  Change of upper spine pitch with respect to time
%
%     3.  Height of the head relative to marker above the junction
%         between the Sacral and caudal vertebrae (lower spine)
%
%     4.  Horizontal distance between the upper and lower spine
%


% Rear: Feature 1

% Time Series
eind = 1; % Event index
rind = [Trial.stc{'r'}(eind,1)-120:Trial.stc{'r'}(eind,2)+120];
ts = rind./fang.sampleRate;
figure,plot(ts,fang(rind,'spine_middle','spine_upper',2))
xlabel('Time (s)')
ylabel('Pitch Body Middle to Upper (radians)')
title ('Rear: Feature 1');
ylim([-.8,pi/2])
Lines(ts([121,end-120]),[],'r');

% pdfs for rear vs not rear
figure,hold on
eds = linspace(-.8,pi/2,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(fang(ind,3,4,2),eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(fang(ind,3,4,2),eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('Pitch Body Middle to Upper (radians)');
ylabel('p_{pitch}');
title ('PDF of Rear vs not Rear');
legend({'Not Rear','Rear'});



% Rear: Feature 2

% Time Series
eind = 1; % Event index
rind = [Trial.stc{'r'}(eind,1)-120:Trial.stc{'r'}(eind,2)+120];
ts = rind./fang.sampleRate;
vfang = MTADang('data',diff(fang(:,'spine_middle','spine_upper',2)).*fang.sampleRate,'sampleRate',fang.sampleRate);
figure,plot(ts,vfang(rind))
xlabel('Time (s)')
ylabel('d(Pitch)/dt Body Middle to Upper (rad/sec)')
title ('Rear: Feature 2');
ylim([-6,6])
Lines(ts([121,end-120]),[],'r');

% pdfs for rear vs not rear
figure,hold on
eds = linspace(-6,6,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(vfang(ind),eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(vfang(ind),eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('d(Pitch)/dt Body Middle to Upper (rad/sec)')
ylabel('p_{d(pitch)/dt}');
title ('PDF of Rear vs not Rear');
legend({'Not Rear','Rear'});



% Rear: Feature 2 alt

% Time Series
eind = 1; % Event index
rind = [Trial.stc{'r'}(eind,1)-120:Trial.stc{'r'}(eind,2)+120];
ts = rind(1:end)./fang.sampleRate;
vfang = MTADang('data',log10(abs(diff(fang(:,'spine_middle','spine_upper',2)).*fang.sampleRate)),'sampleRate',fang.sampleRate);
figure,plot(ts,vfang(rind))
xlabel('Time (s)')
ylabel('log10(abs(d(Pitch)/dt)) Body Middle to Upper (rad/sec)')
title ('Rear: Feature 2 alt');
ylim([-3,1])
Lines(ts([121,end-120]),[],'r');

% pdfs for rear vs not rear
figure,hold on
eds = linspace(-5,1,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(vfang(ind),eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(vfang(ind),eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('log10(abs(d(Pitch)/dt)) Body Middle to Upper (rad/sec)')
ylabel('p_{d(pitch)/dt}');
title ('PDF of Rear vs not Rear');
legend({'Not Rear','Rear'},'location','northwest');





% Rear: Feature 3

% Time Series
eind = 1; % Event index
rind = [Trial.stc{'r'}(eind,1)-120:Trial.stc{'r'}(eind,2)+120];
ts = rind./fang.sampleRate;
figure,plot(ts,(fxyz(rind,5,3)-fxyz(rind,1,3))./10)
xlabel('Time (s)')
ylabel('Head Height - Lower Body Height (cm)')
title ('Rear: Feature 3');
ylim([0,25])
Lines(ts([121,end-120]),[],'r');

% pdfs for rear vs not rear
figure,hold on
eds = linspace(0,25,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc((fxyz(ind,5,3)-fxyz(ind,1,3))/10,eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc((fxyz(ind,5,3)-fxyz(ind,1,3))/10,eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('Head Height - Lower Body Height (cm)')
ylabel('p_{d(Height)}');
title ('PDF of Rear vs not Rear');
legend({'Not Rear','Rear'});






% Rear: Feature 4

% Time Series
eind = 1; % Event index
rind = [Trial.stc{'r'}(eind,1)-120:Trial.stc{'r'}(eind,2)+120];
ts = rind./fang.sampleRate;
figure,plot(ts,fang(rind,1,4,3).*cos(fang(rind,1,4,2))./10)
xlabel('Time (s)')
ylabel('Horizontal Distance between Head and Lower Body (cm)')
title ('Rear: Feature 3');
ylim([0,18])
Lines(ts([121,end-120]),[],'r');

% pdfs for rear vs not rear
figure,hold on
eds = linspace(0,18,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(fang(ind,1,4,3).*cos(fang(ind,1,4,2))./10,eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(fang(ind,1,4,3).*cos(fang(ind,1,4,2))./10,eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('Horizontal Distance between Head and Lower Body (cm)')
ylabel('p_{dist}');
title ('PDF of Rear vs not Rear');
legend({'Not Rear','Rear'});







%% Definition of Grooming
% Rat stands on its hind legs and elevates the head about 10 cm.
% 
% Features
%
%     1.  Horizontal distance between the upper and lower spine
%
%     2.  Height of the lower spine
%
%     3.  
%         
%
%


% Groom: Figure 1

% Time Series
eind = 5; % Event index
rind = [Trial.stc{'m'}(eind,1)-240:Trial.stc{'m'}(eind,2)+240];
ts = rind./fang.sampleRate;
figure,plot(ts,fang(rind,1,4,3).*cos(fang(rind,1,4,2))./10)
xlabel('Time (s)')
ylabel('Horizontal Distance between Head and Lower Body (cm)')
title ('Groom: Feature 1');
ylim([5,17])
Lines(ts([241,end-240]),[],'r');


% PDF
figure,hold on
eds = linspace(5,17,100);
ind = Trial.stc{'a-r-m'};
hn = bar(eds,histc(fang(ind,1,4,3)./10,eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'m'};
hs = bar(eds,histc(fang(ind,1,4,3)./10,eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('Horizontal Distance between Head and Lower Body (cm)')
ylabel('p_{dist}');
title ('PDF of Groom vs not Groom and Rear');
legend({'Not Rear','Rear'},'location','northwest');



% Groom: Figure 2

% Time Series
eind = 5; % Event index
rind = [Trial.stc{'m'}(eind,1)-240:Trial.stc{'m'}(eind,2)+240];
ts = rind./fang.sampleRate;
figure,plot(ts,fxyz(rind,1,3)./10)
xlabel('Time (s)')
ylabel('Horizontal Distance between Head and Lower Body (cm)')
title ('Groom: Feature 1');
ylim([0,7])
Lines(ts([241,end-240]),[],'r');

% pdf
figure,hold on
eds = linspace(0,7,100);%ind = Trial.stc{'a-r-m'}; linspace(5,17,100);
ind = Trial.stc{'a-r-m'};
hn = bar(eds,histc(fxyz(ind,1,3)./10,eds)./sum(diff(ind.data,1,2)),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'m'};
hs = bar(eds,histc(fxyz(ind,1,3)./10,eds)./sum(diff(ind.data,1,2)),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
xlabel('Horizontal Distance between Head and Lower Body (cm)')
ylabel('p_{dist}');
title ('PDF of Groom vs not Groom and Rear');
legend({'Not Rear','Rear'},'location','northeast');


ind = Trial.stc{'a-m-r'};
pa = histc(fxyz(ind,1,3)./10,eds)./sum(diff(ind.data,1,2));
ind = Trial.stc{'m'};
pm = histc(fxyz(ind,1,3)./10,eds)./sum(diff(ind.data,1,2));
 
-nansum(pa.*log2(pa))
-nansum(pm.*log2(pm))

sqrt(sum((pa-pm).^2))*100

pa'*pm

sum(sum(pa*pm'))
figure,imagesc(pa*pm')

ind = Trial.stc{'a-m-r'};
pa = normpdf(eds,nanmean(fxyz(ind,2,3)./10),nanvar(fxyz(ind,2,3)./10));
ind = Trial.stc{'m'};
pm = normpdf(eds,nanmean(fxyz(ind,2,3)./10),nanvar(fxyz(ind,2,3)./10));

sum(pa.*log2(pa./pm))

ffxyz = fxyz.copy;
ffxyz.resample(16);
xls = GetSegs(circshift(ffxyz(:,1,3),-8),1:size(ffxyz,1),16,0)./10;

eds = linspace(0,7,50);
pa = zeros([size(xls,2),numel(eds)]);
pm = zeros([size(xls,2),numel(eds)]);
for i = 1:size(xls,2),
    pa(i,:) = normpdf(eds,nanmean(xls(1:8,i)),nanstd(xls(1:8,i)))';
    pm(i,:) = normpdf(eds,nanmean(xls(9:16,i)),nanstd(xls(9:16,i)))';
end
pa = bsxfun(@rdivide,pa,sum(pa,2))+eps;
pm = bsxfun(@rdivide,pm,sum(pm,2))+eps;

dkl = log10(nansum(pa.*log2(pa./pm),2));
dkl = log10(nansum(pm.*log2(pm./pa),2));
figure,plot(abs(log10(clip(dkl,0.001,10000))))
Lines(Trial.stc{'w',ffxyz.sampleRate}(:),[],'m');


figure,
ind = Trial.stc{'s'};
hist2(fxyz(ind,[1,3],3)./10,linspace(0,7,100),linspace(5,13,100));
caxis([0,200])


figure,
ind = Trial.stc{'s'};
hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2))./10,fxyz(ind,1,3)./10],linspace(5,17,100),linspace(0,7,100));
caxis([0,200])

figure,hold on
eds = linspace(.2,1.5,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(fang(ind,1,4,2),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(fang(ind,1,4,2),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;


figure,hist(diff(Trial.stc{'r'}.data,1,2)./xyz.sampleRate,40);
min_dur = min(diff(Trial.stc{'r'}.data,1,2)./xyz.sampleRate); % .7341
min_dur = .5;
% Head is elevated by about 10cm





%% Definition of Turn
% Type 1: Fore arm assisted 
% 
% Type 2: Rotation with hindlimb 
%

figure,plot(circ_dist(circshift(fang(:,1,4,1),-1),fang(:,1,4,1))),Lines(Trial.stc{'n'}(:),[],'g');


%% Definition of Shake
% Type 1: Head shake
%
% Type 2: Body shake
%


%% Inter animal JPDF labeled contours

fet = xyz.copy;
fet.data = [fang(:,1,4,2),[log10(abs(diff(fang(:,3,4,2))));0]];
edgs    = {linspace(0,1.4,75)};
edgs(2) = {linspace(-6,-1,75)};
edc = edgs;
[edc{:}] = get_histBinCenters(edc);
[X,Y] = meshgrid(edc{:});

sts = 'rwnms';
stc = 'rcymg';
hfig = figure(2);

% JPDF - Head/Body speed
%ind = Trial.stc{'a'};
ind = nniz(fet);
b = fet(ind,:);
hist2(b,edgs{1},edgs{2});
xlabel('log10 body pitch (radians)');
ylabel('log10(abs(d(BMBU_p_i_t_c_h)/dt)) log10(rad/sec)');
title({'JPDF of log10(abs(d(BMBU_p_i_t_c_h)/dt)) VS BMBU_p_i_t_c_h',...
       [Trial.filebase ': overlayed with jg05-20120317 labeled states']});

% Feature and states of jg05-20120317 
if strcmp(Trial.name,'jg05-20120317'),
    ofet = fet.copy;
    ostc = MTATrial('jg05-20120317').load('stc');
end

hold on,
for i = 1:numel(sts),
    b = ofet(ostc{sts(i)},:);
    o = hist2(b,edgs{1},edgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[10,10],'linewidth',1.5,'Color',stc(i))
end
caxis([0,200])
legend({'rear','walk','turn','groom','sit'},'location','SouthEast')
hfig.Position  = [100   100   782   629];

saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150724',...
    [Trial.filebase '-BMBUpVSdBMBUdt_R-stc-jg05-20120317.eps']),'eps2')
saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150724',...
    [Trial.filebase '-BMBUpVSdBMBUdt_R-stc-jg05-20120317.png']),'png')


%% END Inter animal JPDF labeled contours 


%% FF
fet = xyz.copy;
fet.data = cos(fang(ind,1,m2,2)).*fang(ind,1,m2,3);
fet.data = fang(:,1,4,2);

figure,hold on,m2 = 4;
eds = linspace(20,250,100);
eds = linspace(20,180,100);
eds = linspace(0,1.4,100);
ind = Trial.stc{'a-r-m'};
hn = bar(eds,histc(fet(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha =  0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(fet(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha =  0;
ind = Trial.stc{'m'};
hs = bar(eds,histc(fet(ind),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .4;
hs.EdgeAlpha =  0;
xlabel('BLBU_r_x_y (mm)');
legend({'EE','Rear','Groom'});




figure,
eds = linspace(20,170,100);
ads = linspace(-0.8,1.5,100);
ind = Trial.stc{'r'};
hist2([fang(ind,3,4,2),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,500])

figure,% Body length xy VS lower spine marker speed
eds = linspace(20,170,100);
ads = linspace(-0.5,2,100);
ind = Trial.stc{'a'};
hist2([fvel(ind,1),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,200])

figure,
eds = linspace(20,170,80);
ads = linspace(-5,.5,80);
ind = Trial.stc{'m'};
hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-20)))),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);
caxis([0,200])


%% FF
vangp = MTADxyz('data',log10(abs(diff(fang(:,3,4,2)))),'sampleRate',xyz.sampleRate);
figure,
ads = linspace(-5,-1.5,100);
eds = linspace(20,170,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
title(ind.label);xlabel('log10(abs(d(BMBU_p_i_t_c_h)/dt))');ylabel('d_x_y(BLBU)/dt');
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
title(ind.label);xlabel('log10(abs(d(BMBU_p_i_t_c_h)/dt))');ylabel('d_x_y(BLBU)/dt');
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
title(ind.label);xlabel('log10(abs(d(BMBU_p_i_t_c_h)/dt))');ylabel('d_x_y(BLBU)/dt');
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),cos(fang(ind,1,4,2)).*fang(ind,1,4,3)],ads,eds);caxis([0,200])
title(ind.label);xlabel('log10(abs(d(BMBU_p_i_t_c_h)/dt))');ylabel('d_x_y(BLBU)/dt');
suptitle([Trial.filebase])

 

figure,
ads = linspace(-4,1,100);
eds = linspace(-1,1.8,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),fang(ind,3,4,2)],ads,eds);caxis([0,200])

figure,
ads = linspace(-4,1,100);
eds = linspace(60,260,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([vangp(ind),fxyz(ind,5,3)],ads,eds);caxis([0,200])


figure,
ads = linspace(20,170,100);
eds = linspace(60,260,100);
subplot(1,4,1); 
ind = Trial.stc{'a'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,2); 
ind = Trial.stc{'m'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,3); 
ind = Trial.stc{'r'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])
subplot(1,4,4); 
ind = Trial.stc{'a-r-m'};hist2([cos(fang(ind,1,4,2)).*fang(ind,1,4,3),fxyz(ind,5,3)],ads,eds);caxis([0,200])






figure,hold on
eds = linspace(-.5,2,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(ang(ind,1,4,2),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
ind = Trial.stc{'r'};
hs = bar(eds,histc(ang(ind,1,4,2),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;




%% Detect Rearing Periods

figure,plot(nunity(bcom(:,3)))


v = xyz.vel(1:8,[1,2]);


%% Find best low pass filter for detecting walk periods
% $$$ dp = []
% $$$ frange = .5:.2:4;
% $$$ for i = frange,
% $$$ vt = v.copy;
% $$$ vt.filter('ButFilter',3,i,'low');
% $$$ vt.data(vt.data<0) = .001;
% $$$ vta = log10(vt(Trial.stc{'a-w'},:));
% $$$ vtw = log10(vt(Trial.stc{'a&w'},:));
% $$$ nw = nniz(vtw);
% $$$ na = nniz(vta);
% $$$ dp(end+1,:) = (mean(vtw(nw,:))-mean(vta(na,:)))./(.5*sqrt(var(vtw(nw,:))+var(vta(na,:))));
% $$$ end


% $$$ [~,fid] = max(dp(:,1));

% $$$ v.filter('ButFilter',3,frange(fid),'low');
v.filter('ButFilter',3,2.5,'low');
% $$$ figure,plot(v.data)
% $$$ Lines(Trial.stc{'w'}(:),[],'b');

v.data(v.data<.1) = .1;
v.data = log10(v.data);

pv =MTADxyz('data',prod(v(:,[1,4,7]),2),'sampleRate',xyz.sampleRate);
pv.data(pv<0)= 1e-15;
pv.data = log10(pv.data);
figure,plot(pv.data);
Lines(Trial.stc{'w'}(:),[],'b');



figure,hold on
eds = linspace(-15,2.5,100);
ind = Trial.stc{'a-w-r-n'};
bar(eds,histc(pv(ind,1),eds),'histc');
ind = Trial.stc{'w+n'};
h = bar(eds,histc(pv(ind,1),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;


%lb pitch vs lb vel
peds = linspace(-10,2.5,100);
veds = linspace(-.5,2,100);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([ang(ind,1,2,2),v(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w+n'};
figure,hist2([ang(ind,1,2,2),v(ind,1)],peds,veds);caxis([0,500])

% 
peds = linspace(110,170,100);
veds = linspace(-10,2.5,100);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([ang(ind,1,4,3),pv(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w+n'};
figure,hist2([ang(ind,1,4,3),pv(ind,1)],peds,veds);caxis([0,500])


%% Turing overlap 
tfet = Trial.xyz.copy;
tfet.data = circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1));
tfet.filter('ButFilter',3,2.5,'low');
tfet.data = log10(abs(tfet.data));

afet = Trial.xyz.copy;
afet.data = circ_dist(ang(:,5,7,1),circshift(ang(:,5,7,1),-1));
afet.filter('ButFilter',3,2.5,'low');
afet.data = log10(abs(afet.data));



figure,hold on
eds = linspace(-7,-1,70);
ind = Trial.stc{'a-r-n-w'};
hn = bar(eds,histc(tfet(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'n'};
hs = bar(eds,histc(tfet(ind),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;


afet = Trial.xyz.copy;
afet.data = fang(:,1,4,3);
ads = linspace(110,170,70);

figure,hold on
ind = Trial.stc{'w'};
hn = bar(eds,histc(afet(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'n'};
hs = bar(eds,histc(afet(ind),eds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;


zfet = Trial.xyz.copy;
zfet.data = fxyz(:,1,3);
zds = linspace(20,60,170);

zfet = Trial.xyz.copy;
zfet.data = circ_dist(fang(:,1,3,1),fang(:,1,4,1));
zds = linspace(-1.3,1.3,170);

figure,hold on
ind = Trial.stc{'w'};
hn = bar(zds,histc(zfet(ind),zds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'n'};
hs = bar(zds,histc(zfet(ind),zds),'histc');
hs.FaceColor = 'g';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;



ind = Trial.stc{'n'};
figure,hist2([tfet(ind),afet(ind)],eds,ads);

hold on
for g = ind.data',
    plot(tfet(round(mean(g))),afet(round(mean(g(2)))),'*y');
    %plot(tfet(g(1)),afet(g(1)),'*c');
    %plot(tfet(g(2)),afet(g(2)),'*y');
end


saveas(gcf,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150724',...
    [Trial.filebase '-turing_jpdf_langS_BMBUdist_turn.png']),'png')

eds = linspace(-.5,2,100);
ind = Trial.stc{'n'};
figure,hist2([fvel(ind,1),fvel(ind,7)],eds,eds);

hold on
plot(fvel(ind(10,:),1),fvel(ind(10,:),7),'m');
plot(fvel(ind(10,1),1),fvel(ind(10,1),7),'*c');
plot(fvel(ind(10,2),1),fvel(ind(10,2),7),'*y');

plot(fvel(ind(25,:),1),fvel(ind(25,:),7),'m');
plot(fvel(ind(25,1),1),fvel(ind(25,1),7),'*c');
plot(fvel(ind(25,2),1),fvel(ind(25,2),7),'*y');

plot(fvel(ind(15,:),1),fvel(ind(15,:),7),'m');
plot(fvel(ind(15,1),1),fvel(ind(15,1),7),'*c');
plot(fvel(ind(15,2),1),fvel(ind(15,2),7),'*y');

hold on
for g = ind.data',
    %plot(fvel(g',1),fvel(g',7),'m');
    plot(fvel(g(1),1),fvel(g(1),7),'*c');
    plot(fvel(g(2),1),fvel(g(2),7),'*y');
end

saveas(gcf,fullfile(Trial.spath,[Trial.filebase '-turing_hbVel-turn_start_stop_more.png']),'png')


%figure,plot(tfet.data)
%Lines(Trial.stc{'n'}(:),[],'g');

fet = xyz.copy;
%fet.data = [fvel(:,1),tfet(:,1)];
fet.data = [fang(:,1,4,2),[log10(abs(diff(fang(:,3,4,2))));0]];
edgs    = {linspace(0.2,1.7,75)};
edgs(2) = {linspace(-7,-1,75)};
% $$$ edgs    = {linspace(-.5,2,75)};
% $$$ edgs(2) = {linspace(-7,-1,75)};
edc = edgs;
[edc{:}] = get_histBinCenters(edc);
[X,Y] = meshgrid(edc{:});

sts = 'rwnms';
stc = 'rcymg';
hfig = figure(2);
clf
% JPDF - Head/Body speed
ind = Trial.stc{'a-r'};
ind = nniz(fet);
b = fet(ind,:);
hist2(b,edgs{1},edgs{2});
xlabel('log10 body pitch (radians)');
ylabel('log10(abs(d(BMBU_p_i_t_c_h)/dt)) log10(rad/sec)');
title({'JPDF of log10(abs(d(BLBU_p_i_t_c_h)/dt)) VS BMBU_p_i_t_c_h',...
       [Trial.filebase ': overlayed with jg05-20120317 labeled states']});

% Feature and states of jg05-20120317 
if strcmp(Trial.name,'jg05-20120317'),
    ofet = fet.copy;
    ostc = MTATrial('jg05-20120317').load('stc');
end

hold on,
for i = 1:numel(sts),
    b = ofet(ostc{sts(i)},:);
    o = hist2(b,edgs{1},edgs{2});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[10,10],'linewidth',1.5,'Color',stc(i))
end
caxis([0,200])
legend({'rear','walk','turn','groom','sit'},'location','SouthEast')
hfig.Position  = [100   100   782   629];

saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150724',...
    [Trial.filebase '-BMBUpVSdBMBUdt_R-stc-jg05-20120317.eps']),'eps2')
saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/p20150724',...
    [Trial.filebase '-BMBUpVSdBMBUdt_R-stc-jg05-20120317.png']),'png')


%%tahuasthensthaoesnuhoae

peds = linspace(100,160,70);
veds = linspace(-2,1,70);
ind = Trial.stc{'a-w-r-n'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,500])
ind = Trial.stc{'w'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'n'};
figure,hist2([fang(ind,1,4,3).*cos(fang(ind,1,4,2)),pv(ind,1)],peds,veds);caxis([0,200])



peds = linspace(-4,0,70);
ind = Trial.stc{'n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'w'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])
ind = Trial.stc{'a-w-r-n'};
figure,hist2([log10(abs(circ_dist(fang(ind,1,4,1),circshift(fang(ind,1,4,1),-10)))),pv(ind,1)],peds,veds);caxis([0,200])

veds = linspace(50,160,70);
peds = linspace(-6,-1,70);
ind = Trial.stc{'n'};
tfet = xyz.copy;
tfet.data = log10(abs(circ_dist(fang(:,1,4,1),circshift(fang(:,1,4,1),-1))));
figure,hist2([tfet(ind),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])
ind = Trial.stc{'w'};
figure,hist2([tfet(ind),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])
ind = Trial.stc{'a-w-r-n'};
figure,hist2([tfet(ind),fang(ind,1,4,3).*cos(fang(ind,1,4,2))],peds,veds);caxis([0,200])


fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0)=.1;
fvel.data = log10(fvel.data);

ind = Trial.stc{'a-w-n-r-m'};
figure,
hist2([log10(abs(fvel(ind,1))),log10(abs([diff(fvel(ind,1));0]))],-.5:.05:2,-4:.05:.5)

ind = Trial.stc{'a-w'};
figure,
hist2([log10(abs(fvel(ind,1))),log10(abs(fvel(ind,7)))],-.5:.05:2,-.5:.05:2)
caxis([0,500])


figure,hold on
eds = linspace(-.5,2,100);
ind = Trial.stc{'a-w-r'};
bar(eds,histc(v(ind,1),eds),'histc');
ind = Trial.stc{'w'};
h = bar(eds,histc(v(ind,1),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;

figure,plot(v(:,1).^2./abs(v(:,1)-v(:,7)));
Lines(Trial.stc{'w'}(:),[],'m');

vco = v.data;
vco = bsxfun(@minus,v.data,mean(v.data)


v.data(v.data<0) = .001;
v.data = log10(v.data);



[bys,bfs,bts,bphi,bfst] = fet_spec(Trial,bcom);

figure,imagesc(bts,bfs,log10(bys(:,:,1,1))'),axis xy,colormap jet

figure,
for i = 1:bys.size(3),
sp(i) = subplot(bys.size(3),1,i);
imagesc(bts,bfs,log10(bys(:,:,i,i))')
axis xy
colormap jet
end
linkaxes(sp,'xy');
Lines(Trial.stc{'w',1}(:),[],'b');
Lines(Trial.stc{'r',1}(:),[],'r');


figure,imagesc(bts,bfs,bys(:,:,1,3)'),axis xy, colormap jet
figure,imagesc(bts,bfs,bys(:,:,2,3)'),axis xy, colormap jet

figure,imagesc(bts,bfs,(bys(:,:,1,3).*bys(:,:,1,2).*bys(:,:,2,3))'),axis xy, colormap jet


figure,imagesc(bts,bfs,(bys(:,:,1,3)+bys(:,:,1,2)+bys(:,:,2,3))'),axis xy, colormap jet
caxis([2,3])
Lines(Trial.stc{'m',1}(:),[],'m');
Lines(Trial.stc{'k',1}(:),[],'r');

figure,imagesc(bts,bfs,log10(bfst(:,:,3))'),axis xy, colormap jet

figure,plot(bcom(:,1,3))
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'r'}(:),[],'r');


xyz.filter('ButFilter',3,2,'low');
a = create(MTADang,Trial,xyz);
h = Trial.transformOrigin(xyz,'head_back','head_front',{'head_left','head_right'});

r = xyz.copy;
r.data = h.roll;


%Groom
ind = Trial.stc{'m'};
figure,
hist2([circ_dist(a(ind,1,3,1),a(ind,4,7,1)),r(ind)-.2],linspace(-pi,pi,70),linspace(-pi/2,pi/2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);

ind = Trial.stc{'w'};
figure,
hist2([a(ind,1,4,3),v(ind,1)],linspace(110,170,70),linspace(-.5,2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);





ind = Trial.stc{'m'};
figure,
hist2([a(ind,1,4,3),a(ind,2,7,3)],linspace(110,170,70),linspace(70,240,70));
colormap jet
caxis([0,300]);

caxis([0,800]);

ind = Trial.stc{'a-m'};
figure,
hist2([a(ind,1,4,3),a(ind,4,7,3)],linspace(100,170,70),linspace(70,170,70));
colormap jet
caxis([0,300]);

ind = Trial.stc{'m'};
figure,
hist2([a(ind,1,4,3),circ_dist(a(ind,1,2,1),a(ind,1,4,1))],linspace(100,170,70),linspace(-pi/2,pi/2,70));
colormap jet
caxis([0,300]);

caxis([0,800]);



figure,hold on
eds = linspace(100,160,100);
ind = Trial.stc{'a-m-s'};
bar(eds,histc(a(ind,1,4,3),eds),'histc');
ind = Trial.stc{'m'};
h = bar(eds,histc(a(ind,1,4,3),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;




figure,hold on
eds = linspace(-pi,pi,100);
ind = Trial.stc{'a-m'};
bar(eds,histc(circ_dist(a(ind,1,2,1),a(ind,1,4,1)),eds),'histc');
ind = Trial.stc{'m'};
h = bar(eds,histc(circ_dist(a(ind,1,2,1),a(ind,1,4,1)),eds),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .4;


%% WALK Features

xyz = Trial.load('xyz');
afet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-5)-circshift(xyz(:,:,[1,2]),5);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
afet.data = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));



mag = zeros([afet.size(1),1]);
for i= 1:afet.size(1),
mag(i) = PPC(afet(i,[1:5,7]));
end

msync = Trial.xyz.sync.copy;
msync.data = msync.sync.data;


man = MTADfet(Trial.spath,Trial.filebase,...
               mag,...
               Trial.xyz.sampleRate,...
               msync,...
               msync(1),...
               [],[],[],'lower_spine_trajectory_yaw_PPC','lsppc','y');
man.save;




save(fullfile(Trial.spath,...
    [Trial.filebase '-walk_fet_ppc.mat']),'mag')

load(fullfile(Trial.spath,...
    [Trial.filebase '-walk_fet_ppc.mat']))


man = Trial.xyz.copy;
man.data = mag;
man.filter('ButFilter',3,1.5,'low');

% $$$ figure,plot(ma)
% $$$ hold on,plot(man.data)

hfig = figure,hold on
eds = linspace(-.2,1,170);
ind = Trial.stc{'a-w-r'};
hn = bar(eds,histc(man(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'w'};
hs = bar(eds,histc(man(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;
title({'Filtered Pair Wise Phase Consistency (PPC)',['of trajectories ' ...
                    'markers along the head and spine']})
ylabel('count');
xlabel('PPC');
legend({'a-r-w','w'})

saveas(hfig,fullfile(hostPath,['PPC-',Trial.filebase,'.eps']),'epsc')
saveas(hfig,fullfile(hostPath,['PPC-',Trial.filebase,'.png']),'png')


figure
sts = {'a','n','a-r-w','w'};
for i = 1:numel(sts)
subplot(2,2,i)
ind = Trial.stc{sts{i}};
eds = linspace(-.2,1,100);
vds = linspace(-.5,2,100);
hist2([man(ind),fvel(ind,1)],eds,vds);
caxis([0,100])    
title(sts{i})
xlabel('ppc')
ylabel('log10 LS speed')
end



hold on,
for i = 1:80
ind = Trial.stc{'w'}(i,:);
plot(man(ind(1)),fvel(ind(1),1),'*g')
plot(man(ind(end)),fvel(ind(end),1),'*r')
plot(man(ind),fvel(ind,1),'.w')
pause(.1)
end


dman = man.copy;
dman.data = [diff(man.data);0];

figure
ind = Trial.stc{'w'};
eds = linspace(-.2,1,100);
vds = linspace(-.015,.015,100);
hist2([man(ind),dman(ind,1)],eds,vds);
caxis([0,200])    


hold on,
for i = 1:80
ind = Trial.stc{'w'}(i,:);
plot(man(ind(1)),dman(ind(1),1),'*g')
plot(man(ind(end)),dman(ind(end),1),'*r')
%plot(man(ind),dman(ind,1),'.w')
pause(.1)
end



figure
ind = Trial.stc{'a-r-w'};
eds = linspace(-.015,.015,100);
vds = linspace(-1.5,1.5,100);
hist2([dman(ind),fac(ind,1)],eds,vds);
caxis([0,100])    

hold on,
for i = 1:80
ind = Trial.stc{'w'}(i,:);
plot(dman(ind(1)),fac(ind(1),1),'*g')
plot(dman(ind(end)),fac(ind(end),1),'*r')
%plot(fvel(ind),fac(ind,1),'.w')
pause(.1)
end


tff = MTADxyz('data',...
              log10(abs(circ_dist(circshift(fang(:,1,4,1),-20),circshift(fang(:,1,4,1),20)))),'sampleRate',xyz.sampleRate);


figure
sts = {'a','n','a-r-w-n-m-k','w'};
for i = 1:numel(sts)
subplot(2,2,i)
ind = Trial.stc{sts{i}};
%eds = linspace(-4,1,100);
eds = linspace(-.2,1,100);
%vds = linspace(-.2,1,100);
vds = linspace(-.8,2,100);
%hist2([tff(ind),man(ind,1)],eds,vds);
%hist2([tff(ind),fvel(ind,1)],eds,vds);
hist2([man(ind),fvel(ind,1)],eds,vds);
caxis([0,100])    
title(sts{i})
%xlabel('ang speed log10')
xlabel('ppc')
%ylabel('ppc')
ylabel('BL log10 speed')
end

hist2([tff(ind),man(ind,1)],eds,vds);
caxis([0,100])    

figure,
ind = Trial.stc{'w'};
dind = diff(ind.data,1,2);
ind = ind(dind>60,:);
eds = linspace(-.8,2,100);
vds = linspace(-.2,1,100);
hist2([fvel(ind,1),man(ind,1)],eds,vds);
caxis([0,50])    
title(sts{i})




figure
ind = Trial.stc{'a'};
eds = linspace(110,170,100);
vds = linspace(-.2,1,100);
hist2([fang(ind,1,4,3),man(ind,1)],eds,vds);
caxis([0,50])    
figure,ind = Trial.stc{'a-r-w'}; hist2([fang(ind,1,4,3),man(ind,1)],eds,vds); caxis([0,50])    
figure,ind = Trial.stc{'w'}; hist2([fang(ind,1,4,3),man(ind,1)],eds,vds); caxis([0,50])    

figure
ind = Trial.stc{'a'};
eds = linspace(110,280,100);
vds = linspace(-.2,1,100);
hist2([fang(ind,1,7,3),man(ind,1)],eds,vds);
caxis([0,50])    
figure
ind = Trial.stc{'w'};
hist2([fang(ind,1,7,3),man(ind,1)],eds,vds);
caxis([0,50])    




figure
ind = Trial.stc{'a'};
eds = linspace(40,120,100);
vds = linspace(60,150,100);
hist2([fxyz(ind,2,3),fang(ind,2,4,3)],eds,vds);
caxis([0,150])    

figure
ind = Trial.stc{'a'};
eds = linspace(-1,2,100);
vds = linspace(-.6,1.6,100);
hist2([fvel(ind,11),fang(ind,3,4,2)],eds,vds);
caxis([0,250]) 



%% Distribution of behaviors for speed of ACOM

aph = .4;
hfig = figure,hold on

ind = Trial.stc{'s'};
hr = bar(eds,histc(fvel(ind,11),eds),'histc');
hr.FaceColor = 'm';
hr.FaceAlpha = .2;
hr.EdgeAlpha = 0;

ind = Trial.stc{'m'};
ha = bar(eds,histc(fvel(ind,11),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
hs = bar(eds,histc(fvel(ind,11),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = aph;
hs.EdgeAlpha = 0;

ind = Trial.stc{'r'};
hr = bar(eds,histc(fvel(ind,11),eds),'histc');
hr.FaceColor = 'y';
hr.FaceAlpha = .7;
hr.EdgeAlpha = 0;

ind = Trial.stc{'n'};
hr = bar(eds,histc(fvel(ind,11),eds),'histc');
hr.FaceColor = 'g';
hr.FaceAlpha = .8;
hr.EdgeAlpha = 0;

legend({'s','m','w','r','n',});
title('Center of mass Spine & Head Speed');
ylabel('Sample Count');
xlabel('Speed 1og10(cm/s)');

saveas(hfig,fullfile(hostPath,['ACOM_speed-',Trial.filebase,'.eps']),'epsc')
saveas(hfig,fullfile(hostPath,['ACOM_speed-',Trial.filebase,'.png']),'png')


%% Distribution of behaviors for speed of Lower Spine
aph = .4;
hfig = figure,hold on

ind = Trial.stc{'s'};
hr = bar(eds,histc(fvel(ind,1),eds),'histc');
hr.FaceColor = 'm';
hr.FaceAlpha = .2;
hr.EdgeAlpha = 0;

ind = Trial.stc{'m'};
ha = bar(eds,histc(fvel(ind,1),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
hs = bar(eds,histc(fvel(ind,1),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = aph;
hs.EdgeAlpha = 0;

ind = Trial.stc{'r'};
hr = bar(eds,histc(fvel(ind,1),eds),'histc');
hr.FaceColor = 'y';
hr.FaceAlpha = .7;
hr.EdgeAlpha = 0;

ind = Trial.stc{'n'};
hr = bar(eds,histc(fvel(ind,1),eds),'histc');
hr.FaceColor = 'g';
hr.FaceAlpha = .8;
hr.EdgeAlpha = 0;

legend({'s','m','w','r','n',});
title('Center of mass Spine & Head Speed');
ylabel('Sample Count');
xlabel('Speed 1og10(cm/s)');

saveas(hfig,fullfile(hostPath,['LS_speed-',Trial.filebase,'.eps']),'epsc')
saveas(hfig,fullfile(hostPath,['LS_speed-',Trial.filebase,'.png']),'png')

%% Distribution of behaviors for speed of Lower Spine
aph = .4;
hfig = figure,hold on
eds = linspace(-.2,1,100);

ind = Trial.stc{'s'};
hr = bar(eds,histc(man(ind,1),eds),'histc');
hr.FaceColor = 'm';
hr.FaceAlpha = .2;
hr.EdgeAlpha = 0;

ind = Trial.stc{'m'};
ha = bar(eds,histc(man(ind,1),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
hs = bar(eds,histc(man(ind,1),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = aph;
hs.EdgeAlpha = 0;

ind = Trial.stc{'r'};
hr = bar(eds,histc(man(ind,1),eds),'histc');
hr.FaceColor = 'y';
hr.FaceAlpha = .7;
hr.EdgeAlpha = 0;

ind = Trial.stc{'n'};
hr = bar(eds,histc(man(ind,1),eds),'histc');
hr.FaceColor = 'g';
hr.FaceAlpha = .8;
hr.EdgeAlpha = 0;

legend({'s','m','w','r','n',});
title('Center of mass Traj PPC & Head Speed');
ylabel('Sample Count');
xlabel('PPC');

saveas(hfig,fullfile(hostPath,['TRAJ_PPC-',Trial.filebase,'.eps']),'epsc')
saveas(hfig,fullfile(hostPath,['TRAJ_PPC-',Trial.filebase,'.png']),'png')

ads = linspace(-.2,1,100);
eds = linspace(-.8,2,100);
figure,
subplot(121)
ind = Trial.stc{'a-r-w-n'};
hist2([man(ind),fvel(ind,1)],ads,eds)
caxis([0,100])
subplot(122)
ind = Trial.stc{'n'};
hist2([man(ind),fvel(ind,1)],ads,eds)
caxis([0,100])



ffet = fxyz.copy;
ffet.data = nunity(ffet.data);

figure,hold on,
eds= linspace(-6,6,200);
ind = Trial.stc{'a-w-r'};
ha = bar(eds,histc(ffet(ind,1,3),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
hs = bar(eds,histc(ffet(ind,1,3),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


figure
subplot(121)
ind = Trial.stc{'a-r-w'};
eds = linspace(-.8,2,100);
vds = linspace(-5,5,100);
hist2([fvel(ind,1),ffet(ind,1,3)],eds,vds);
caxis([0,250])    
subplot(122)
ind = Trial.stc{'w'};
hist2([fvel(ind,1),ffet(ind,1,3)],eds,vds);
caxis([0,250])    



mxyz = xyz.copy;
mxyz.filter('ButFilter',3,14,'low');
mang = create(MTADang,Trial,mxyz);

daf = Trial.xyz.copy;
daf.data = sq(mang(:,1,1:4,3).*sin(mang(:,1,1:4,2)));

daf.data = bsxfun(@minus,xyz(:,2:4,3),xyz(:,1,3))
daf.filter('ButFilter',3,8,'low');
daf.data = diff(daf.data);
p
figure,plot(daf.data),
Lines(Trial.stc{'k'}(:),[],'m');
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'n'}(:),[],'g');

sig = diff(nunity(xyz(:,4,3)-xyz(:,1,3)));

figure,plot(diff(circ_dist(mang(:,1,2,1),mang(:,1,4,1))))


figure,plot(log10(sq(sum(GetSegs(circshift(diff(circ_dist(mang(:,1,3,1),mang(:,2,4,1))),20),1:mang.size(1),40,nan).^2))))
Lines(Trial.stc{'k'}(:),[],'m');
Lines(Trial.stc{'w'}(:),[],'c');
Lines(Trial.stc{'n'}(:),[],'g');


sfet = Trial.xyz.copy;
sfet.data = log10(sq(sum(GetSegs(circshift(diff(circ_dist(mang(:,1,3,1),mang(:,2,4,1))),20),1:mang.size(1),40,nan).^2)))';
wfet = Trial.xyz.copy;
wfet.data = log10(sq(sum(GetSegs(circshift(diff(mxyz(:,1,3)),20),1:mang.size(1),40,nan).^2)))';

eds= linspace(-6,0,200);

figure,hold on
ind = Trial.stc{'a-k-w-r'};
noise = sfet(ind);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = sfet(ind);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


dp = abs((nanmean(signal)-nanmean(noise))/(.5*sqrt(nanvar(signal)+nanvar(noise))));


ads = linspace(-.2,1,100);
eds = linspace(-5,3,100);

figure,
subplot(121)
ind = Trial.stc{'a-r-w-k-n'};
hist2([man(ind),wfet(ind)],ads,eds)
caxis([0,100])
subplot(122)
ind = Trial.stc{'n'};
hist2([man(ind),wfet(ind)],ads,eds)
caxis([0,100])


ads = linspace(-.8,2,100);
eds = linspace(-5,3,100);

figure,
subplot(121)
ind = Trial.stc{'a-r-w'};
hist2([fvel(ind,1),wfet(ind)],ads,eds)
caxis([0,100])
subplot(122)
ind = Trial.stc{'w'};
hist2([fvel(ind,1),wfet(ind)],ads,eds)
caxis([0,100])

ads = linspace(-.8,2,100);
eds = linspace(-.2,1,100);

figure,
subplot(121)
ind = Trial.stc{'a-r-w-k-n'};
hist2([fvel(ind,1),man(ind)],ads,eds)
caxis([0,100])
subplot(122)
ind = Trial.stc{'w'};
hist2([fvel(ind,1),man(ind)],ads,eds)
caxis([0,100])


eds= linspace(-.2,1.2,200);
figure,hold on
ind = Trial.stc{'a-r'};
noise = circ_dist(ang(ind,1,3,2),ang(ind,2,4,2));
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'r'};
signal = circ_dist(ang(ind,1,3,2),ang(ind,2,4,2));
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


dp = abs((nanmean(signal)-nanmean(noise))/(.5*sqrt(nanvar(signal)+nanvar(noise))));



bb =bsxfun(@minus,[circshift(xyz(:,1,:),-30),circshift(xyz(:,1,:),30)],xyz(:,1,:));

% $$$ %figure,plot((sqrt(sum(diff(sq(cross(bb(:,1,:),bb(:,2,:),3))).^2,2))))
% $$$ %Lines(Trial.stc{'w'}(:),[],'m');
% $$$ %hold on,plot(sqrt(sum(sq(diff(bb(:,1,:)).^2),2)).*100)
% $$$ 
% $$$ [mode,wsig,sampleRate,,overwrite] = ...
% $$$      DefaultArgs(varargin,{'raw',true,120,...
% $$$ 
% $$$ defspec = struct('nFFT',2^7,'Fs',fet.sampleRate,...
% $$$                  'WinLength',2^6,'nOverlap',2^6*.875,...
% $$$                  'FreqRange',[1,20]);                    
% $$$ sg = MTADxyz('data',sqrt(sum(diff(sq(cross(bb(:,1,:),bb(:,2,:),3))).^2,2)),'sampleRate',xyz.sampleRate);
% $$$ [ys,fs,ts] = fet_spec(Trial,sg,'mtchglong',false);
% $$$ 
% $$$ %figure,imagesc(ts,fs,log10(ys)'),axis xy,colormap jet
                    
ya = ys.copy;
ya.data = log10(nanmean(ys.data,2));

eds= linspace(-9,4,200);
figure,hold on
ind = Trial.stc{'a-w-r'};
noise = ya(ind);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .3;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = ya(ind);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .3;
hs.EdgeAlpha = 0;
 
                    


% $$$ dps = [];
% $$$ for sh = 5:5:70;
% $$$     dp = [];
% $$$     for fv  = .5:.25:4,

sh = 40;
fv = 2.5;


        bb =bsxfun(@minus,[circshift(xyz(:,1,:),-sh),circshift(xyz(:,1,:),sh)],xyz(:,1,:));
        mfet = Trial.xyz.copy;
        mfet.data = sqrt(sum(diff(sq(cross(bb(:,1,:),bb(:,2,:),3))).^2,2));
        mfet.filter('ButFilter',3,fv,'low');
        mfet.data(mfet.data<0) = 1e-9;
        mfet.data = log10(mfet.data);
% $$$ mfet.data = diff(mfet.data);
% $$$ mfet.data = log10(mean((mfet.segs(1:mfet.size(1),35,nan).^2))');

eds= linspace(0,3,200);
eds= linspace(-4,3,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = mfet(ind);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .3;
ha.EdgeAlpha = 0;

ind = Trial.stc{'a&w'};
signal = mfet(ind);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .3;
hs.EdgeAlpha = 0;

dp = abs((nanmean(signal)-nanmean(noise))/(.5*sqrt(nanvar(signal)+nanvar(noise))));
% $$$         dp(end+1) = abs((nanmean(signal)-nanmean(noise))/(.5*sqrt(nanvar(signal)+nanvar(noise))));
% $$$     end
% $$$     dps = cat(1,dps,dp);
% $$$ end

man = Trial.load('fet','lsppc');
man.filter('ButFilter',3,2,'low');

ads = linspace(-4,3,100);
eds = linspace(-.2,1,100);

figure,
subplot(121)
ind = Trial.stc{'a-r-w-n'};
hist2([mfet(ind,1),man(ind)],ads,eds)
caxis([0,100])
subplot(122)
ind = Trial.stc{'a&n'};
hist2([mfet(ind,1),man(ind)],ads,eds)
caxis([0,100])


eds = linspace(-.8,2,100);
figure,
ind = Trial.stc{'r'};
hist2([mfet(ind,1),fvel(ind,1)],ads,eds)
caxis([0,200])

ads = linspace(-.2,1,100);
eds = linspace(-.8,2,100);
figure,
subplot(121)
ind = Trial.stc{'a-w-r'};
hist2([man(ind,1),fvel(ind,1)],ads,eds)
caxis([0,200])
subplot(122)
ind = Trial.stc{'w'};
hist2([man(ind,1),fvel(ind,1)],ads,eds)
caxis([0,200])

Lines(Trial.stc{'w'}(:),[],'m');





figure,plot(cross(circshift(xyz(:,1,:),-sh)-circshift(xyz(:,1,:),sh),fxyz(:,4,:)-fxyz(:,1,:)));
sh = 1;
figure,plot(sqrt(sum(cross(circshift(xyz(:,1,:),-sh)-circshift(xyz(:,1,:),sh),fxyz(:,4,:)-fxyz(:,1,:)).^2,3)))
hold on,plot(sqrt(sum(cross(circshift(xyz(:,1,:),-sh)-circshift(xyz(:,1,:),sh),fxyz(:,7,:)-fxyz(:,1,:)).^2,3)))
Lines(Trial.stc{'w'}(:),[],'m');

figure,plot(dot(circshift(xyz(:,1,[1,2]),-sh)-circshift(xyz(:,1,[1,2]),sh),fxyz(:,4,[1,2])-fxyz(:,1,[1,2]),3))

%% GOOD FET walk and turn
nfet = MTADxyz('data',dot(repmat(circshift(xyz(:,1,[1,2]),-sh)-circshift(xyz(:,1,[1,2]),sh),[1,5,1]),fxyz(:,[2:5,7],[1,2])-repmat(fxyz(:,1,[1,2]),[1,5,1]),3),'sampleRate',Trial.xyz.sampleRate);
nfet.filter('ButFilter',3,2.4,'low');
% $$$ hold on,plot(nfet.data);
% $$$ Lines(Trial.stc{'w'}(:),[],'m');

nfet.data = log10(abs(nanmean(nfet.data,2))+1).*sign(nanmean(nfet.data,2));

%nfet.data = log10(var(nfet.data,[],2)).*sign(nanmean(nfet.data,2));

figure,
ind = Trial.stc{'w'};
hist2([nfet(ind,1),log10(fet(ind,15))],linspace(-4,4,100),linspace(-6,0,100));
caxis([0,100])

figure,
ind = Trial.stc{'a-r-n-w'};
hist2([nfet(ind,1),log10(mfet(ind))],linspace(-4,5,100),linspace(-4,5,100));
caxis([0,300])


eds= linspace(-5,6,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = nfet(ind,1);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = nfet(ind,1);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;



dp = abs((nanmean(signal)-nanmean(noise))/(.5*sqrt(nanvar(signal)+nanvar(noise))));

%% BUTT WAG

Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');

nfet = MTADfet(Trial.spath,Trial.filebase,...
               sqrt(sum((xyz(:,1:4,[1,2])-fxyz(:,1:4,[1,2])).^2,3)),...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'ButtWag','bw','b');


nfet.filter('ButFilter',3,8,'low');
nfet.data = [diff(nfet.data);0,0,0,0];

dsa =  struct('nFFT',2^7,'Fs',nfet.sampleRate,...
              'WinLength',2^6,'nOverlap',2^6*.875,...
                            'FreqRange',[1,40]);
[rhm,fs,ts,phi,fst] = fet_spec(Trial,nfet,'mtchglong',false,'defspec',dsa);

figure,
imagesc(ts,fs,log10(rhm(:,:,1,1)')),axis xy,colormap jet,
Lines(Trial.stc{'w',1}(:),[],'m');

figure,
imagesc(ts,fs,phi(:,:,1,3)'),axis xy,colormap hsv,
Lines(Trial.stc{'w',1}(:),[],'m');

figure,
eds= linspace(-8,.5,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = log10(rhm(ind,5));
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = log10(rhm(ind,5));
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


nfet = MTADfet(Trial.spath,Trial.filebase,...
               [diff(sqrt(sum((xyz(:,1,:)-fxyz(:,1,:)).^2,3)));0],...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'ButtWag','bw','b');
nfet.filter('ButFilter',3,8,'low');
nfet.data = diff(nfet.data);

ns = MTADxyz('data',log10(sq(mean(nfet.segs(1:nfet.size(1),50,nan).^2))),'sampleRate',xyz.sampleRate);


eds= linspace(-9,.7,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = ns(ind,1);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = ns(ind,1);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


eds = linspace(-8,.7,100);
ads = linspace(-.2,1,100);

figure,
subplot(221)
ind = Trial.stc{'a'};
hist2([ns(ind),man(ind)],eds,ads);
caxis([0,100])
subplot(222)
ind = Trial.stc{'w'};
hist2([ns(ind),man(ind)],eds,ads);
caxis([0,100])
subplot(223)
ind = Trial.stc{'a-w-r-n'};
hist2([ns(ind),man(ind)],eds,ads);
caxis([0,100])


wag = MTADxyz('data',circ_dist(ang(:,1,4,1),fang(:,1,4,1))-circ_dist(ang(:,1,3,1),fang(:,1,3,1)),'sampleRate',xyz.sampleRate);
wag.filter('ButFilter',3,8,'low');


%% Shake discrimination
Nang = MTADxyz(Trial.spath,Trial.filebase,...
               circ_dist(ang(:,1,3,1),ang(:,2,4,1)),...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'bshake','bs','s');

Nang.filter('ButFilter',3,[1,20],'bandpass');


figure,plot(diff(Nang.data)),Lines(Trial.stc{'w'}(:),[],'m');

dsa =  struct('nFFT',2^7,'Fs',nfet.sampleRate,...
              'WinLength',2^6,'nOverlap',2^6*.875,...
                            'FreqRange',[1,20]);
[rhm,fs,ts,phi,fst] = fet_spec(Trial,Nang,'mtchglong',true,'defspec',dsa);

figure,
imagesc(ts,fs,log10(rhm(:,:,1,1)'))
axis xy,colormap jet,
nLines(Trial.stc{'m',1}(:),[],'m');

msk = ones([rhm.size]);
msk(:,[1:5,16:20]) = -1;

[scpow,scm,scs] = nunity(mean(rhm.data.*msk,2));

figure,hold on
Lines(Trial.stc{'k',1}(:),[0,100],'m');
plot(ts,scpow)

%% END Shake discrimination

figure,plot(diff(wag.data))
hold on,plot(nfet.data/20)
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');



eds= linspace(-.7,2,200);
figure,hold on
ind = Trial.stc{'a-w-r-n'};
noise = fvel(ind,1);
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;

ind = Trial.stc{'w'};
signal = fvel(ind,1);
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


tsh = 30;
tang = MTADang('data',circ_dist(circshift(fang(:,1,4,1),tsh),circshift(fang(:,1,4,1),-tsh)),...
               'sampleRate',ang.sampleRate);
tsh = 20;
tang = MTADang('data',circ_dist(circshift(ang(:,1,4,1),tsh),circshift(ang(:,1,4,1),-tsh)),...
               'sampleRate',ang.sampleRate);




eds= linspace(-3,.8,200);
figure,hold on
ind = Trial.stc{'w'};
noise = log10(abs(tang(ind)));
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;
ind = Trial.stc{'n'};         
signal = log10(abs(tang(ind)));
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;

ads = linspace(-3,.8,70);
eds = linspace(-.8,2,70);
figure,
ind = Trial.stc{'n'};
hist2([log10(abs(tang(ind,1))),fvel(ind,1)],ads,eds)
caxis([0,400])


tsh = 1;
afet = Trial.xyz.copy;
bfet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(xyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));

gfet = Trial.xyz.copy;
gfet.data = unwrap(circ_dist(afet(:,10),ang(:,1,4,1)));
gfet.filter('ButFilter',3,.1,'high');


eds= linspace(.5,2.6,200);
m = 1;
figure,hold on
ind = Trial.stc{'a-r'};
noise = log10(fxyz(ind,7,3));
ha = bar(eds,histc(noise,eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .5;
ha.EdgeAlpha = 0;
ind = Trial.stc{'r'};         
signal = log10(fxyz(ind,7,3));
hs = bar(eds,histc(signal,eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;

zfet = Trial.xyz.copy;
zfet.data = log10(sqrt(diff(fxyz(:,5,3)).^2));
eds =  linspace(1.5,2.6,100);
ads =  linspace(-4,1.5,100);
ind = Trial.stc{'a'};
figure
hist2([zfet(ind),log10(fxyz(ind,5,3))],ads,eds);
caxis([0,100])

rmax = [];
for r = Trial.stc{'r'}.data',
    rmax(end+1) = max(xyz(r',7,3));
end


%eds = linspace(1.4,2.5,200);
eds = linspace(-.8,1.5,200);
figure,bar(eds,histc(log10(rmax),eds),'histc');

sts = 'smwrn';
stc = 'bbbrb';

hfig = figure;hold on
for s = 1:numel(sts);
ind = Trial.stc{sts(s)};
smax = [];
for r = ind.data',
    smax(end+1) = max(ang(r',3,4,2));
    %smax(end+1) = log10(max(xyz(r',7,3)));
end
hr = bar(eds,histc(smax,eds),'histc');
hr.FaceColor = stc(s);
hr.FaceAlpha = .3;
hr.EdgeAlpha = 0;
end

legend({'s','m','w','r','n',});
title('Center of mass Spine & Head Speed');
ylabel('Sample Count');
xlabel('Speed 1og10(cm/s)');

saveas(hfig,fullfile(hostPath,['ACOM_speed-',Trial.filebase,'.eps']),'epsc')
saveas(hfig,fullfile(hostPath,['ACOM_speed-',Trial.filebase,'.png']),'png')

end


fvel = xyz.vel([],[1,2]);
fvel.filter('ButFilter',3,2.5,'low');
fvel.data(fvel.data<0)=.1;
fvel.data = log10(fvel.data);
%figure,plot(log10(fvel.data+1));

figure,hold on,
eds = linspace(-1,2,200);
ind = Trial.stc{'a-w-r'};
hr = bar(eds,histc(fvel(ind,1),eds),'histc');
hr.FaceColor = 'c';
hr.FaceAlpha = .3;
hr.EdgeAlpha = 0;
ind = Trial.stc{'w'};
hr = bar(eds,histc(fvel(ind,1),eds),'histc');
hr.FaceColor = 'r';
hr.FaceAlpha = .3;
hr.EdgeAlpha = 0;


figure,plot(diff(circ_dist(circshift(fang(:,1,3,1),-1),circshift(fang(:,1,3,1),1))))
figure,plot(circ_dist(ang(:,1,3,1),ang(:,2,4,1)))
afetcirc_dist(ang(:,1,3,1),ang(:,2,4,1))

figure,plot(circ_dist(circshift(ang(:,2,4,1),-1),circshift(ang(:,2,4,1),1)))
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

Nang = MTADxyz(Trial.spath,Trial.filebase,...
               ang(:,1,10,3),...
               xyz.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'bscratch','br','r');

Nang.filter('ButFilter',3,[1,20],'bandpass');
figure,plot(diff(Nang.data))

figure,plot(ang(:,10,12,2))


dsa =  struct('nFFT',2^8,'Fs',nfet.sampleRate,...
              'WinLength',2^7,'nOverlap',2^7*.875,...
                            'FreqRange',[1,20]);
[rhm,fs,ts,phi,fst] = fet_spec(Trial,Nang,'mtchglong',true,'defspec',dsa);

figure,
imagesc(ts,fs,log10(rhm(:,:,1,1)'))
axis xy,colormap jet,
Lines(Trial.stc{'m',1}(:),[],'m');

msk = ones([rhm.size]);
msk(:,[1:8,16:20]) = -1;

[scpow,scm,scs] = nunity(mean(rhm.data.*msk,2));
scpow = median(log10(rhm.data(:,8:15)),2);

figure,hold on
Lines(Trial.stc{'m',1}(:),[-6,-3],'m');
Lines([],-3,'m');
plot(ts,log10(rhm(:,10)))
plot(ts,scpow)


%% Start here

tsh = 1;
afet = Trial.xyz.copy;
afet.data = xyz(:,:,[1,2])-fxyz(:,:,[1,2]);
bfet = Trial.xyz.copy;
%afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(fxyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));

%bfet.filter('ButFilter',3,2.3,'low');


Mang = MTADfet(Trial.spath,Trial.filebase,...
               circ_dist(circshift(afet(:,[1:5,7]),-5),circshift(afet(:,[1:5,7]),5)),...
               afet.sampleRate,...
               Trial.sync.copy,...
               Trial.sync(1),...
               [],[],[],'mshake','mr','m');
% $$$ 
% $$$ dsa =  struct('nFFT',2^7,'Fs',Mang.sampleRate,...
% $$$               'WinLength',2^6,'nOverlap',2^6*.875,...
% $$$                             'FreqRange',[1,20]);
% $$$ 
% $$$ [rhm,fs,ts,phi,fst] = fet_spec(Trial,Mang,'mtchglong',true,'defspec',dsa);
% $$$ 
% $$$ figure,
% $$$ imagesc(ts,fs,log10(rhm(:,:,1,1)'))
% $$$ axis xy,colormap jet,
% $$$ 
% $$$ figure,plot(ts,skewness(rhm(:,:,1,1),[],2))
% $$$ Lines(Trial.stc{'w',1}(:),[],'m');
% $$$ 
% $$$ figure,plot(ts,log10(median(rhm(:,1:8,1,1),2)./median(rhm(:,9:end,1,1),2)))
% $$$ hold on,plot(ts,log10(median(rhm(:,1:8,2,2),2)./median(rhm(:,9:end,2,2),2)))
% $$$ hold on,plot(ts,log10(median(rhm(:,1:8,3,3),2)./median(rhm(:,9:end,3,3),2)))
% $$$ Lines(Trial.stc{'w',1}(:),[],'m');
% $$$ 
% $$$ figure,
% $$$ imagesc(ts,fs,rhm(:,:,1,3)')
% $$$ axis xy,colormap jet,
% $$$ Lines(Trial.stc{'w',1}(:),[],'m');


% $$$ Mangs = Mang.copy;
% $$$ Mangs.data = Mang.segs(1:Mang.size(1),40,nan);
% $$$ mag = zeros([Mangs.size(2),1]);
% $$$ for i = 1:Mangs.size(2),
% $$$     mag(i) = PPC(Mangs(:,i,[1]));
% $$$ end
% $$$ save(fullfile(Trial.spath,...
% $$$     [Trial.filebase '-walk_fet_another_ppc.mat']),'mag')

load(fullfile(Trial.spath,...
    [Trial.filebase '-walk_fet_another_ppc.mat']))
maw = Trial.xyz.copy;
maw.data = mag;
maw.filter('ButFilter',3,1.5,'low');

figure,plot(maw.data)
Lines(Trial.stc{'w'}(:),[],'m');



m = MTADxyz('data',circ_dist(afet(:,1),ang(:,1,4,1)),'sampleRate',Trial.xyz.sampleRate);
m.data = circ_dist(circshift(m.data,-5),circshift(m.data,5));
m.data = [diff(m.data);0];
wn = 40;
mv = m.copy;
mv.data = circshift(log10(1./permute(mean(m.segs(1:m.size(1),wn,nan).^2),[2,3,4,1])),-round(wn/2));


figure,hold on
eds = linspace(-1.5,3,100);
ind = Trial.stc{'a-w-n-r'};
hn = bar(eds,histc(mv(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'w'};
hs = bar(eds,histc(mv(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;


eds = linspace(-1.5,3,100);
edy = linspace(-.2,1,100);
ind = Trial.stc{'a'};
%ind = Trial.stc{'w'};
%ind = Trial.stc{'a-w'};
%ind = Trial.stc{'m'};
figure,
hist2([mv(ind),man(ind,1)],eds,edy)
caxis([0,50])


eds = linspace(-1.5,3,100);
edy = linspace(-.8,2,100);
ind = Trial.stc{'a'};
ind = Trial.stc{'w'};
% $$$ ind = Trial.stc{'a-w'};
% $$$ ind = Trial.stc{'m'};
figure,
hist2([mv(ind),fvel(ind,1)],eds,edy)
caxis([0,50])


eds = linspace(-.8,2,100);
edy = linspace(-.05,1,100);
ind = Trial.stc{'a-w'};
ind = Trial.stc{'a'};
ind = Trial.stc{'n'};
ind = Trial.stc{'m'};
ind = Trial.stc{'r'};
ind = Trial.stc{'w'};
figure,
hist2([fvel(ind,1),maw(ind,1)],eds,edy)
caxis([0,100])

eds = linspace(-.8,2,100);
edy = linspace(-.05,1,100);
ind = Trial.stc{'a-w'};
ind = Trial.stc{'a'};
ind = Trial.stc{'n'};
ind = Trial.stc{'m'};
ind = Trial.stc{'r'};
ind = Trial.stc{'w'};
figure,
hist2([fvel(ind,1),man(ind,1)],eds,edy)
caxis([0,100])


eds = linspace(-.2,1,100);
edy = linspace(-.2,1,100);
ind = Trial.stc{'a-w'};
ind = Trial.stc{'n'};
ind = Trial.stc{'w'};
ind = Trial.stc{'a-w-n'};
figure,
hist2([man(ind,1),maw(ind,1)],eds,edy)
caxis([0,100])



eds = linspace(-8,-1,100);
edy = linspace(-.05,1,100);
ind = Trial.stc{'a'};
ind = Trial.stc{'a-w'};
ind = Trial.stc{'n'};
ind = Trial.stc{'w'};
ind = Trial.stc{'a-w-n'};
figure,
hist2([ns(ind,1),maw(ind,1)],eds,edy)
caxis([0,100])
p

eds = linspace(-5,.5,100);
edy = linspace(-.2,1,100);
ind = Trial.stc{'a'};
nnnind = Trial.stc{'a-n'};
ind = Trial.stc{'n'};
ind = Trial.stc{'w'};
ind = Trial.stc{'a-w-n'};
figure,
hist2([tff(ind,1),man(ind,1)],eds,edy)
caxis([0,100])




fxyz = xyz.copy;
fxyz.filter('ButFilter',3,5,'low');


dmat = markerDiffMatrix(fxyz);

%mx = mean(fxyz(:,[1,3],:),2);


tmat = [];
for i = 1:3,
    cmat = sq(cross(dmat(:,12,i,:),dmat(:,12,4,:),4));
    nmat = MTADxyz('data',sqrt(sum(cmat.^2,2)),'sampleRate',xyz.sampleRate);
    tmat = cat(2,tmat,nmat.data);
end


%figure,
figure,plot(sum(diff(tmat),2))
%hold on,plot(diff(nmat.data))
Lines(Trial.stc{'w'}(:),[],'m');
%Lines(Trial.stc{'r'}(:),[],'r');

figure,hold on
eds = linspace(10,10000,100);
ind = Trial.stc{'a-r'};
hn = bar(eds,histc(nmat(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'r'};
hs = bar(eds,histc(nmat(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;






figure,hold on
eds = linspace(0,2,100);
ind = Trial.stc{'a-w-n-r'};
hn = bar(eds,histc(log10(bfet(ind,1)+1),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .6;
hn.EdgeAlpha = 0;
ind = Trial.stc{'w'};
hs = bar(eds,histc(log10(bfet(ind,1)+1),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeAlpha = 0;




edx = linspace(-.8,2,100);
edy = linspace(0,1.5,100);

ind = Trial.stc{'a'};
figure,hist2([fvel(ind,1),log10(bfet(ind,1)+1)],edx,edy);
caxis([0,200])



edx = linspace(-2,1,100);
edy = linspace(-2,1,100);
ind = Trial.stc{'a'};
figure,hist2([bfet(ind,1),bfet(ind,7)],edx,edy);
caxis([0,200])

xs = circshift(xyz.data,-10)-circshift(xyz.data,10);

nn = MTADxyz('data',dot(sq(xs(:,1,:)),sq(bsxfun(@rdivide,dmat(:,2,11,:),sqrt(sum(dmat(:,2,11,:).^2,4)))),2),'sampleRate',xyz.sampleRate);;
nns = sign(nn.data);
nn.data = log10(abs(nn.data)+1).*nns;


nh = MTADxyz('data',dot(sq(xs(:,7,:)),sq(bsxfun(@rdivide,dmat(:,2,11,:),sqrt(sum(dmat(:,2,11,:).^2,4)))),2),'sampleRate',xyz.sampleRate);;
nhs = sign(nh.data);
nh.data = log10(abs(nh.data)+1).*nhs;

figure,hold on
eds = linspace(-2,2.8,200);
ind = Trial.stc{'a-r-w'};
hn = bar(eds,histc(nn(ind),eds),'histc');
hn.FaceColor = 'c';
hn.FaceAlpha = .4;
hn.EdgeAlpha = 0;
ind = Trial.stc{'w'};
hs = bar(eds,histc(nn(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeAlpha = 0;


edx = linspace(-2,2.5,70);
edy = linspace(-8,-1,70);
ind = Trial.stc{'w'};
figure,
hist2([nn(ind),ns(ind,1)],edx,edy);
caxis([0,200])


