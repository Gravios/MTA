

%% figure 1
Trial = MTATrial('Ed10-20140812');
Trial = MTASession('jg05-20120317');
Trial = MTASession('er06-20130612');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;

xyz = Trial.load('xyz');
xyz.filter(gtwin(.3,xyz.sampleRate));
ang = Trial.ang.copy;
ang.create(Trial,xyz);
figPath = '/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Figures/Figure_1';
figPath = '/home/gravio/Documents/man2015-jgEd-MoCap/Figures/Figure_1';

%% Fig:1:A Recording setup ( Maze and Cameras )
hfig = figure(83832);
pic_path = '';
img = imread(pic_path);
imagesc(img);

%% Fig:1:B Rat with Markers
hfig = figure(83832);
pic_path = '';
img = imread(pic_path);
imagesc(img);



%% Fig:1:C Marker Skeleton Recontruction



%% Fig:1:D Head motion and feature types
    Trial = MTATrial('jg05-20120317');
    fs = []; ts = [];

    xyz = Trial.load('xyz');
    % create a ridgid body model
    rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
    % find the center of mass of the model
    hcom = xyz.com(rb);

    % add coordinates of the model's center of mass to the xyz object
    xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
                  ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));

    xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);

    nm = cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_right',:)-hcom);
    nm = bsxfun(@rdivide,nm,sqrt(sum((nm).^2,3))).*20+hcom;

    xyz.addMarker('htx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);

    nm = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
    nm = bsxfun(@rdivide,nm,sqrt(sum((nm).^2,3))).*20+hcom;
    
    xyz.addMarker('hrx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);

    nm = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
    nm = bsxfun(@rdivide,nm,sqrt(sum((nm).^2,3))).*20+hcom;
    
    xyz.addMarker('hbx',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},nm);
    
    ind = 1000;
    figure,plot3(nm(ind,1,1),nm(ind,1,2),nm(ind,1,3),'.m')
    hold on,plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.b')
    hold on,plot3(xyz(ind,5,1),xyz(ind,5,2),xyz(ind,5,3),'.b')
    hold on,plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.g')
    hold on,plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.r')
    
    
    xyz.filter(gausswin(5)./sum(gausswin(5)));

    ang = create(Trial.ang.copy,Trial,xyz);
    ang.data(~nniz(ang(:,1,2,1)),:,:,:)=0;
    bang = [];
% $$$     bang = [bang,ButFilter(ang(:,'hbx','fhcom',3),3,[2,55]./(ang.sampleRate/2),'bandpass')];
% $$$     bang = [bang,ButFilter(ang(:,'hrx','fhcom',3),3,[2,55]./(ang.sampleRate/2),'bandpass')];
% $$$     bang = [bang,ButFilter(ang(:,'htx','fhcom',3),3,[2,55]./(ang.sampleRate/2),'bandpass')];

    bang = [bang,ang(:,'hbx','fhcom',3)];
    bang = [bang,ang(:,'hrx','fhcom',3)];
    bang = [bang,ang(:,'htx','fhcom',3)];
    bang = [zeros([1,size(bang,2)]);diff(bang)];
    bang = [bang,circ_dist(ang(:,'hbx','hcom',2),ang(:,'hbx','fhcom',2))];
    bang = [bang,circ_dist(ang(:,'hbx','hcom',1),ang(:,'hbx','fhcom',1))];
    bang = [bang,fet_roll(Trial,[],'raw')];
     
    bfet = Trial.xyz.copy;
    bfet.data = bang;

    [ys,fs,ts,phi,fst] = fet_spec(Trial,bfet,'mtchglong',true,'overwrite',true);

    xyz.resample(ys);
    wang = create(Trial.ang.copy,Trial,xyz);    
    
    c = 4;
    f = 35;
    nind = nniz(wang(:,1,2,1));
    figure,hist2([wang(nind,5,7,2),log10(ys(nind,f,c,c))],linspace(-1.5,1.5,100),linspace(-7,-2,100));
    s = 'w';
    nind = Trial.stc{s};
    figure,hist2([wang(nind,5,7,2),log10(ys(nind,f,c,c))],linspace(-1.5,1.5,100),linspace(-7,-2,100));
    
    
    nc = ys.size(3);
    figure,sp = [];
    for i = 1:nc,
        sp(i) = subplot(nc,1,i);imagesc(ts,fs,log10(ys(:,:,i,i))'),axis xy
        if i<4,caxis([-7,-2]),else,caxis([-9,-4]),end
    end
    linkaxes(sp,'xy');

    ind = 1:3:ys.size(1);
    ind = resample(Trial.stc{'w'}.cast('TimeSeries'),ys)&nniz(ys(:,1,1,1));

    ind = Trial.stc{'a'};    
    figure,plot3(log10(ys(ind,40,1,1)),log10(ys(ind,40,2,2)),log10(ys(ind,40,3,3)),'.');
    line([-9,0],[-9,0],[-9,0])    
    c = 2;
    f = 35;
    figure,hist2(log10([ys(:,35,1,1),ys(:,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'m'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'w'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    nind = Trial.stc{'r'};
    figure,hist2(log10([ys(nind,35,1,1),ys(nind,f,c,c)]),linspace(-7,-2,100),linspace(-7,-2,100));

    
    vfet = vel(Trial.load('xyz')
    [ys,fs,ts,phi,fst] = fet_spec(Trial,bfet,'mtchglong',true,'overwrite',true);
    
     pbins = linspace(-7,-3,100);
     pfd = histc(log10(ys(:,:,1,1)),pbins,1);
     figure,imagesc(pbins,fs,pfd'),axis xy


     hfig = figure(848283); 
     for s = 'arwms';
         clf(hfig);
         nind = Trial.stc{s};
         pfd = histc(log10(ys(nind,:,1,1)),pbins,1);
         imagesc(pbins,fs,pfd'),axis xy
         title(['PSD PDF by freq for ' nind.label ' peridos']);
         xlabel('Power (Roentgens)');
         ylabel('Frequency');
         saveas(hfig,fullfile(figPath,['Fig1D-sts-' nind.label '.png']),'png');
         saveas(hfig,fullfile(figPath,['Fig1D-sts-' nind.label '.eps']),'eps2');
     end
     
     
     fet = fet_lgr(Trial);
     nind = nniz(fet);

     
     
 
     [isig] = fastica(cov(fet(nind,:)));
     nind = nniz(fet);figure,N = hist2([fet(nind,:)*isig(11,:)',fet(nind,:)*isig(2,:)'],linspace(-10,10,100),linspace(-15,15,100));



%% Fig:1:E
% marker error ridgid body marker distance PSD
Trial = MTATrial('jg05-20120317');
ang = create(Trial.ang.copy,Trial,Trial.load('xyz'));
ang.data = ang(Trial.stc{'a'},'head_back','head_front',3)-...
           median(ang(Trial.stc{'a'},'head_back','head_front',3));
pstr = struct('nFFT',2^11,'Fs',ang.sampleRate,...
              'WinLength',2^10,'nOverlap',2^10*.875,...
              'FreqRange',[0.1,40]);
[ys,fs,ts] = fet_spec(Trial,ang,'mtchglong',false,[],pstr);
pedges = linspace(-7,-1,125);
N = histc(log10(ys(nniz(ys),:)),pedges,1);
figure,imagesc(fs,pedges,N),axis xy;
xlabel('Frequency (Hz)')
ylabel('log10 mm^2/Hz')
title({'PSD of distance between two markers','of a ridgid body structure'})



Trial = MTATrial('jg05-20120317');
ang = create(Trial.ang.copy,Trial,Trial.load('xyz'));
bt = circ_dist(ang(:,'spine_lower','spine_middle',1),ang(:,'spine_middle','head_front',1));
vng = Trial.ang.copy;
vng.data = [0;diff(bt)];
[ys,fs,ts] = fet_spec(Trial,vng,'mtchglong',false);
figure,imagesc(ts,fs,log10(ys(:,:))'),axis xy;caxis([-9,-4])

[A,W] = fastica(log10(ys(Trial.stc{'a'},:))');

[A,W] = fastica(log10(ys({Trial.stc{'w'}},:))');

figure,hist2([log10(ys(Trial.stc{'a'},:))*A(:,166),log10(ys(Trial.stc{'a'},:))*A(:,95)],100,100)
v1 = linspace(-8,-4,100);
v2 = linspace(-90,-50,100);
figure,hist2([log10(ys(Trial.stc{'a'},:))*A(:,166),log10(ys(Trial.stc{'a'},:))*A(:,95)],v1,v2)
figure,hist2([log10(ys(Trial.stc{'m'},:))*A(:,166),log10(ys(Trial.stc{'m'},:))*A(:,95)],v1,v2)

sfet = (median(ys.data(:,30:60),2)+median(ys.data(:,78:100),2))./(mean(ys(:,1:25),2)+mean(ys(:,65:75),2)+mean(ys(:,105:166),2));
ss = ys.copy;
ss.data = log10(sfet);

figure,plot(ss.data)
Lines(Trial.stc{'m',ys.sampleRate}(:),[],'r')

vfet = vel(resample(Trial.load('xyz').filter(gtwin(.3,Trial.xyz.sampleRate)),ys.sampleRate),7);
vfet.data = log10(vfet.data);

nind = Trial.stc{'a'};
figure, hist2([vfet(nind),ss(nind)],linspace(-2,2,100),linspace(-2,1.5,100))



hfig = figure(838884);
edges = 43:.05:45;
%rTrial.stc{'w'}
rper = Trial.stc{'r'}.cast('TimeSeries');
wper = Trial.stc{'w'}.cast('TimeSeries');
vel = xyz.vel(7,[1,2]);
%ind = ~(rper.data|wper.data);
ind = vel<2;
N = histc(ang(ind,5,7,3),edges);
bar(edges,N,'histc')
title('Distance between the back and front head markers')
xlabel('Inter Marker Distance (mm)');
ylabel('Count');

vel = xyz.vel([1:9]);
vel.data(vel<0.01) = 0.01;
vel.data = log10(vel.data);

m = {'head_back','head_right'};mb = [];
vbins = -.5:.05:2;
[~,mb(:,1)] = histc(vel(:,m{1}),vbins);
[~,mb(:,2)] = histc(vel(:,m{2}),vbins);

mb = MTADxyz('data',mb,'sampleRate',xyz.sampleRate);

ind = resample(Trial.stc{'r'}.cast('TimeSeries'),xyz)&nniz(mb);


A = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@mean);
B = accumarray(mb(ind,:),ang(ind,m{1},m{2},3),repmat(numel(vbins),[1,2]),@std);
S = accumarray(mb(ind,:),ind(ind),repmat(numel(vbins),[1,size(mb,2)]),@sum);
A(S<100)=nan;
B(S<100)=nan;

% mean marker distance vs speed
figure,
subplot(121)
imagescnan({vbins,vbins,A'},prctile(A(nniz(A(:))),[5,95]),[],true,[0,0,0]),axis xy,
subplot(122)
imagescnan({vbins,vbins,B'},[0.4,5],[],true,[0,0,0]),axis xy

figure,
hist(ang(Trial.stc{'a'},'head_back','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_front',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_right','head_left',3),20:.02:37)
figure,
hist(ang(Trial.stc{'a'},'head_back','head_right',3),10:.02:37)


%% Fig:1:F Breathing (maybe)

Trial = MTASession('Ed10-20140815');

xyz = Trial.load('xyz').filter(gtwin(.25,Trial.xyz.sampleRate));;
ang = Trial.ang.copy;
ang.create(Trial,xyz);


ind = 866330:867110;
mar = {'pelvis_root','spine_upper'};

hfig = figure(8482838);
sp=[];
sp(1) = subplot(211);
plot(round((ind-ind(1))/xyz.sampleRate,2),ang(ind,mar{1},mar{2},3));
%ylim([131.9,132.6])
title('Distance Between the Pelvis and the Upper Spine')
ylabel('Distance (mm)')
sp(2) = subplot(212);
plot(round((ind-ind(1))/xyz.sampleRate,2),lfp.data(ind));
title('Nasal Cavity Pressure Sensor')
ylabel('NCP (mV)')
ylim([-4000,2800])
xlabel('Time (s)')
linkaxes(sp,'x')
xlim([.25,6]);

saveas(hfig,[fullfile(figPath,'Fig1F-alt1.png'),'png');
saveas(hfig,[fullfile(figPath,'Fig1F-alt1.eps'),'eps2');


%% new Figure2 
Trial = MTATrial('jg05-20120317');
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2';

%Features 
xyz = Trial.load('xyz').filter(gtwin(.1,Trial.xyz.sampleRate));

vl = vel(xyz,1:8,[1,2]);
vl.data = ButFilter(vl.data,3,4/(vl.sampleRate/2),'low');

[rhm,fs,ts] = fet_rhm(Trial,[],'wcsd');


ang = create(Trial.ang.copy,Trial,xyz);

sfet = Trial.ang.copy;
sfet.data = sum(abs([circ_dist(ang(:,1,3,1),ang(:,1,2,1)),...
                     circ_dist(ang(:,2,4,1),ang(:,2,3,1)),...
                     circ_dist(ang(:,3,5,1),ang(:,3,4,1)),...
                     circ_dist(ang(:,4,7,1),ang(:,4,5,1))]),2);

swg = Trial.ang.copy;
swg.data = circ_dist(ang(:,1,4,1),ang(:,2,4,1));


%% mutinfo crap
v = log10(vl.data);

edges = linspace(-2,2,64);
sbound = -30:30;
ixy = zeros([numel(sbound),size(v,2),size(v,2)]);

padding = [0,0];
vind = logical(subsref(cast(resample(Trial.stc{'a'}+padding,xyz),'TimeSeries'),substruct('.',{'data'})));
nind = numel(vind);

s = 1;
for m = 1:size(v,2)
for o = 1:size(v,2)
for shift = sbound
[out,xb,yb,p]=hist2([v(vind,m),circshift(v(vind,o),shift)],edges,edges);
pxy = out./nind;
px = histc(v(vind,m),xb);
px = px(1:end-1)/nind;
py = histc(circshift(v(vind,o),shift),yb);
py = py(1:end-1)/nind;
ixy(s,m,o) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
s = s+1;
end
s = 1;
end
end

[mixy,sixy] = max(ixy);
mixy = sq(mixy);
sixy = sq(sixy)-ceil(numel(sbound)/2);




%% Fig Parameters
sts = 'rwnms';
stc = 'rcymg';
edgs    = {linspace(-2,2,75)};
edgs(2) = {linspace(-2,2,75)};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
ns = 6;

smat = stc2mat(Trial.stc);;
Trial.stc{'w'}(1:3:end,1)
xpers = bsxfun(@plus,Trial.stc{'w',1}(1:3:end,1),[-15,15]);
xpers(xpers(:,1)<1,:) = [];
xpers((xpers(:,2)-Trial.sync(end))>0,:) = [];
s = 20;



hfig = figure(2);


% JPDF - Head/Body speed
subplot2(ns,4,[3,4],4);
b = log10([median(vl(nniz(vl),[1:2]),2),median(vl(nniz(vl),[5:8]),2)]);
hist2(b,edgs{1},edgs{2});
xlabel('log10 body speed (cm/s)');
ylabel('log10 head speed (cm/s)');
title('JPDF of log10 head and body speeds');

hold on,
for i = 1:numel(sts),
    b = log10([median(vl(Trial.stc{sts(i)},1:2),2),median(vl(Trial.stc{sts(i)},5:8),2)]);
    o = hist2(b,linspace(-1.5,2,75),linspace(-1.5,2,75));
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[20,20],'linewidth',1.5,'Color',stc(i))
end


ndag = double(~eye(8));
subplot2(ns,4,[1,2],4);
imagesc(mixy(:,:).*ndag);
colorbar
title('mutual information between marker speeds');
set(gca,'YtickMode','manual');
set(gca,'Ytick',1:8);
set(gca,'YtickLabelMode','manual');
set(gca,'YtickLabel',vl.model.ml('short'));

subplot2(ns,4,[5,6],4);
imagesc(sixy(:,:)/vl.sampleRate*1000);
colorbar
title('time lag of maximum mutual information (ms)')


ns = 6;
for s = 1:size(xpers,1),
    ind = round(xpers(s,:).*xyz.sampleRate);
    ind = ind(1):ind(2);
    i = 1;
    
    %figure(201) % SPEED head and body
    subplot2(ns,4,i,[1:3]);i=i+1;
    plot(ind/vl.sampleRate,[median(vl(ind,1:2),2),median(vl(ind,5:8),2)]),axis tight
    title('xy speed of head and body')
    %xlabel('Time (s)')
    ylabel('Speed (cm/s)');
    ylim([0,80])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(202) % DIRECTION head and body
    subplot2(ns,4,i,[1:3]);i=i+1;
    plot(ind/vl.sampleRate,[ang(ind,1,3,1),ang(ind,4,7,1)]),axis tight
    title('Direction of head and body')
    %xlabel('Time (s)')
    ylabel('Direction (radians)');
    ylim([-pi,pi])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(203)%  PITCH SLPR and SMSU
    subplot2(ns,4,i,[1:3]);i=i+1;
    plot(ind/vl.sampleRate,[ang(ind,1,2,2),ang(ind,3,4,2)]),axis tight
    title('Pitch of SLPR and SMSU')
    %xlabel('Time (s)')
    ylabel('Pitch (radians)');
    ylim([-pi/2,pi/2])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
    
    %figure(204)
    subplot2(ns,4,i,[1:3]);i=i+1;
    plot(ind/vl.sampleRate,[ang(ind,4,5,3)]),axis tight%,ang(ind,4,7,3)])
    title('Intermarker Distance of head and body')
    %xlabel('Time (s)')
    ylabel('Distance (mm)');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});

    
    %figure(205)
    subplot2(ns,4,i,[1:3]);i=i+1;
    plot(ind/vl.sampleRate,abs(sfet(ind,:))),axis tight%,ang(ind,4,7,3)])
    title('Cummulative Spine Angle')
    %xlabel('Time (s)')
    ylabel('Cummulative Angle (radians)');
    ylim([0,8])
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTickLabel',{});
   
    
    %figure(206) RHM Rhythmic Head Motion
    subplot2(ns,4,i,[1:3]);i=i+1;
    sind = round(xpers(s,:).*rhm.sampleRate);
    sind = sind(1):sind(2);
    title('RHM Rhythmic Head Motion')
    imagesc(ts(sind),fs,log10(rhm(sind,:))'),axis xy ,caxis([-5.2,-2.5])
    xlabel('Time (s)')
    
    
    saveas(hfig,fullfile(figPath,['Fig2-Features-sample_' num2str(s) '_' Trial.filebase '.png']),'png');
    saveas(hfig,fullfile(figPath,['Fig2-Features-sample_' num2str(s) '_' Trial.filebase '.eps']),'eps');


end

















%% Figure 2 Trajectories and behavioral labeling


Trial = MTATrial('jg05-20120317');
Stc = Trial.load('stc','hand_labeled_rev1'); 
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_3';
%figPath = '/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Figures/Figure_2';
%exPer = [26664,27100];
%exPer = [26000,26480];
%exPer = [25200,25580];
exPer = [55200,60000];
xyz = Trial.load('xyz').filter(gtwin(.25,Trial.xyz.sampleRate));
ang = create(Trial.ang.copy,Trial,xyz);
stateColors = 'brcgym';


hfig = figure(38239384);clf
set(hfig,'position',[836   110   775   792]);

%% Fig:2:A Skeleton examples


axes('Position',[0.13,0.5,0.775,0.45]);hold on;
set(gca,'CameraPositionMode', 'manual'                    ,...
	'YLim', [-200 200],...
        'XLim', [-300 400],...
	'ZLim', [0 300],...
	'CameraPosition',     [-1909.49 3535.01 1621.08],...
	'CameraTargetMode',   'manual'                    ,...
	'CameraTarget',       [50 0 150]                ,...
 	'CameraUpVectorMode', 'manual'                    ,...
 	'CameraUpVector',     [0 0 1]                   ,...
 	'CameraViewAngleMode','manual'                  ,...
 	'CameraViewAngle',    [6.88708]);
daspect([1,1,1])


Stc = Trial.load('stc','hand_labeled_rev1'); 
nper = Stc{'n'}&exPer;
wper = Stc{'w'}&exPer;
rper = Stc{'r'}&exPer;
sper = Stc{'s'}&exPer;

for i= [1:4,5,7],
p=plot3(xyz(exPer(1):exPer(2),i,1),xyz(exPer(1):exPer(2),i,2),xyz(exPer(1):exPer(2),i,3),'.k');
set(p,'MarkerSize',4)
end

pMode = [];%'surface';
plotSkeleton(xyz,exPer(1),pMode,ang);              % Skeleton @ Begining of trajectory
plotSkeleton(xyz,round(mean(nper.data)),pMode,ang);% Skeleton @ During Turn
plotSkeleton(xyz,round(mean(wper.data)),pMode,ang);% Skeleton @ During walk
plotSkeleton(xyz,exPer(2),pMode,ang);              % Skeleton @ end of trajectory



if ~sper.isempty, 
    for s = 1:sper.size(1),
        for i= [1:4,5,7],
            p=plot3(xyz(sper(s,:),i,1),xyz(sper(s,:),i,2),xyz(sper(s,:),i,3),'.c');
            set(p,'MarkerSize',4)
        end
    end
end
if ~nper.isempty, 
    for s = 1:nper.size(1),    
        for i= [1:4,5,7],
            p=plot3(xyz(nper(s,:),i,1),xyz(nper(s,:),i,2),xyz(nper(s,:),i,3),'.g');
            set(p,'MarkerSize',4)
        end
    end
end
if ~wper.isempty, 
    for s = 1:wper.size(1),
        for i= [1:4,5,7],
            p=plot3(xyz(wper(s,:),i,1),xyz(wper(s,:),i,2),xyz(wper(s,:),i,3),'.b');
            set(p,'MarkerSize',4)
        end,
    end
end
if ~rper.isempty,
    for s = 1:rper.size(1),
        for i= [1:4,5,7],
            p = plot3(xyz(rper(s,:),i,1),xyz(rper(s,:),i,2),xyz(rper(s,:),i,3),'.r');
            set(p,'MarkerSize',4)
        end
    end
end


zlim([0,300]);
print('-depsc','-r600',fullfile(figPath,['Fig2A_skeleton_test_20150218Z2032']))
% $$$ set(gca,'YTickLabelMode','manual');
% $$$ set(gca,'YTickLabel',{});
% $$$ set(gca,'XTickLabelMode','manual');
% $$$ set(gca,'XTickLabel',{});
% $$$ set(gca,'ZTickLabelMode','manual');
% $$$ set(gca,'ZTickLabel',{});

clf

%figure
% Fig:2:B - feature matrix
figName = 'Fig2B_feature_matrix';
fet = fet_lgr(Trial);

%subplot2(10,1,[1:4],1)
axes('Position',[0.13,0.25,0.775,0.2])
ts = (1:fet.size(1))./fet.sampleRate;
per = round(exPer./xyz.sampleRate)+[-10,10];
ind = ts>per(1)&ts<per(2);
ufet = nunity(fet);
imc = imagesc(ts(ind),1:fet.size(2),ufet(ind,:)');
caxis([-2,2]);
xlim(round(exPer./xyz.sampleRate)+[-10,10])
Lines(round(exPer./xyz.sampleRate),[],'k');
%xlabel('Time (s)')
flabels = {'speed SL'   ,...
           'speed SM'   ,...
           'speed HF'   ,...
           'height SL'  ,...
           'Z-diff SL_HF',...
           'pitch SL_PR',...
           'pitch SM_SU',...
           'dist SL_PR' ,...
           'dist PR_SM' ,...
           'dist SU_HB' ,...
           'dist SL_HB' ,...
           'av_SLSM_SMHF'};
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',1:numel(flabels),...
        'YTickLabel',flabels);
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
%saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
%saveas(hfig,fullfile(figPath,[figName,'.eps']),'eps2');




% Fig:2:C - Expert Labels
stateLabels = {'walk','rear'};
stateColors = 'br';

figName = 'Fig2C_expert_labels';
Stc = Trial.load('stc','hand_labeled_rev1');
%hfig = figure(10161);
%subplot2(10,1,5,1)
axes('Position',[0.13,0.20,0.775,0.04])
plotSTC(Stc,1,'patch',stateLabels,'br');
xlim(round(exPer./xyz.sampleRate)+[-10,10])
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'EXP'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
%title('Expert Labels');
%saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
%saveas(hfig,fullfile(figPath,[figName,'.eps']),'eps2');



% Fig:2:D - Empirical Labels
figName = 'Fig2D_empirical_labels';
stateColors = 'br';
Stc = Trial.load('stc','auto_wbhr');
%hfig = figure(10171);
%subplot2(10,1,6,1)
axes('Position',[0.13,0.15,0.775,0.04])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+[-10,10])
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'EMP'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
%title('Empirical Labels');
%saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
%saveas(hfig,fullfile(figPath,[figName,'.eps']),'eps2');



% Fig:2:E - LGR Model Labels
figName = 'Fig2E_LGR_labels';
Stc = Trial.load('stc','LGR_wrsnkm');
%hfig = figure(10181);
%subplot2(10,1,7,1)
axes('Position',[0.13,0.10,0.775,0.04])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+[-10,10])
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'LGR'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
%title('Logistic Regression Labels');
%saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
%saveas(hfig,fullfile(figPath,[figName,'.eps']),'eps2');


% Fig:2:F - LDA Model Labels
figName = 'Fig2F_LDA_labels';
Stc = Trial.load('stc','LDA_wrsnkm');
%hfig = figure(10191);
axes('Position',[0.13,0.05,0.775,0.04])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+[-10,10])
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'LDA'});
xlabel('Time (s)');
%saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
%saveas(hfig,fullfile(figPath,[figName,'.eps']),'eps2');


export_fig(fullfile(figPath,['Fig2A-F_complete_test_20150218Z1347']),'eps','-r300')
print('-depsc','-r300',fullfile(figPath,['Fig2A-F_complete_test_20150218Z1218']))


% Fig:2:G - LDA Model Classifier scores
figName = 'Fig2G_LDA_clsfr';
[Stc,d_state] = bhv_lda(Trial,false);
hfig = figure(10194);
plot((1:size(d_state,1))/30,d_state)
xlim(round(exPer./xyz.sampleRate)+[-10,10])
ylim([-100,10])
ylabel('A.U.');
xlabel('Time (s)');
saveas(hfig,fullfile(figPath,[figName,'.png']),'png');
saveas(hfig,fullfile(figPath,[figName,'.eps']),'epsc');


%% NEED CONFUSION MATRIX STATS
%Newest lgr model as of 2014.02.17
[Stc,d_state] = bhv_lgr(Trial,false,'model_name','MTAC_testLGRnewfet');

% Training
bhv_lda(Trial,true,'fet',fet,'model_name','MTAC_testLDAnewfet');

% Classification
[Stc,d_state] = bhv_lda(Trial,'fet','fet_lgr','model_name','MTAC_testLDAnewfet');









%% Figure 3 JPDFs
figPath = '/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Figures/Figure_3';

Trial = MTATrial('jg05-20120310');

% Vars of interest 
%Trial.load('stc','hand_labeled_rev1');
xyz = Trial.load('xyz');
ang = Trial.load('ang');
vxy = vel(Trial.load('xyz').filter(gtwin(.25,xyz.sampleRate)),{'spine_lower','head_front'},[1,2]);
vxy.data(vxy.data<.001) = 0.001;
rol = fet_roll(Trial,[],'default');


% Fig:3:A - Rearing from everything else
tag = 'Rearing-Everything';
v1 = Trial.ang.copy;
v1.data = ang(:,2,3,3);
v2 = Trial.ang.copy;
v2.data = ang(:,3,4,2);
bhv_JPDF(Trial,v1,v2,70,70,[],[],...
         'Distance(Pelvis,Spine Middle) (mm)',...
         'Upper Spine Pitch (rad)',{'a-r','r'},tag)


tag = 'pelvis2spineM_vs_headroll';
v1 = Trial.ang.copy;
v1.data = ang(:,2,3,3);pp
v2 = rol.copy;
bhv_JPDF(Trial,v1,v2,70,70,[],[],...
         'Distance(Pelvis,Spine Middle) (mm)',...
         'Head roll (rad)','rwhl',tag)


% Fig:3:B - walk from everything else
tag = 'speed_SpineL_vs_HeadF';
v1 = vxy.copy;
v1.data = log10(vxy(:,1));
v2 = vxy.copy;
v2.data = log10(vxy(:,2));
bhv_JPDF(Trial,v1,v2,70,70,[-2,2],[-2,2],...
         'Speed (Lower Spine) log10(mm)',...
         'Speed (Head Front) log10(mm)',{'a-r','w','a-w-r'},...
         tag)


% Fig:3:C - walk to high and low walk
tag = 'Walk_High-low';
v1 = Trial.ang.copy;
v1.data = ang(:,5,7,3);
v2 = Trial.ang.copy;
v2.data = ang(:,4,5,3);
bhv_JPDF(Trial,v1,v2,70,70,[],[],...
         'Distance(SpineU,HeadB) (mm)',...
         'Head roll (rad)',{'w','h','l'},tag)


tag = 'HeightSL_vs_headPitch';
v1 = Trial.xyz.copy;
v1.data = xyz(:,1,3);
v2 = Trial.ang.copy;
v2.data = ang(:,5,7,2);
bhv_JPDF(Trial,v1,v2,70,70,[],[],...
         'Hight(SpineL) (mm)',...
         'Head Pitch (rad)','whl',tag)




% Fig:3:D - d-prime metrics


dprime = (mean(SelectPeriods(xyz(:,7,3),Trial.stc{'r'}.data,'c',0))+mean(xyz(Trial.stc{'r'},7,3)))...
          ./sqrt((var(SelectPeriods(xyz(:,7,3),Trial.stc{'r'}.data,'c',0))+var(xyz(Trial.stc{'r'},7,3))).*.5);



[~,dstates] = bhv_lgr(Trial,false);
figure,hist(log10(dstates(:,2)),1000)
figure,hist2(log10(dstates(:,[1,6])),linspace(-10,0,70),linspace(-8,0,70)),caxis([0,600])






%% Figure 5 RHM (rythmic head motion) feature versus NCP (nasal cavity pressure)
figPath = '/gpfs01/sirota/homes/gravio/Documents/Manuscripts/Vicon_Methods_2015/Figures/Figure_5';

Trial = MTATrial('Ed10-20140820','all','rov');

%generate features

[rhm,fs,ts] = fet_rhm(Trial,[],'csd');
ncp = fet_ncp(Trial,[],'csd',66);
%plot features with linked axes

drhm = fet_rhm(Trial,[],'raw');
dncp = fet_ncp(Trial,[],'raw',66);


% Fig:5:A - Spectrums of the rhythmic head motion and the nasal
%             cavity pressure sensor 
figFileName = 'RHM_NCP_spec_ex2';
hfig = figure(3929439);
sp(1) = subplot(211);
imagesc(ts,fs,log10(rhm.data)'),axis xy,caxis([-4.7,-2.5])
title('Rhythmic Head Motion (RHM)')
ylabel('frequency (Hz)')
xlabel('Time (s)')
sp(2) = subplot(212);
imagesc(ts,fs,log10(ncp.data)'),axis xy,caxis([4.8,7])
title('Nasal Cavity Pressure (NCP)')
ylabel('frequency (Hz)')
xlabel('Time (s)')
linkaxes(sp,'xy');
xlim([600,860])
%xlim([700,860])
%xlim([2991,3151])
saveas(hfig,fullfile(figPath,[figFileName '.png']),'png');
saveas(hfig,fullfile(figPath,[figFileName '.eps']),'eps2');


% Fig:5:B - Mean Coherence as a function of head pitch and
%             frequency 
Trial = MTATrial('Ed10-20140812');
figFileName = ['bhv_rhm_ncp_distrb_' Trial.filebase '_ex1'];
hfig = bhv_rhm_ncp_distrb(Trial,[],[],66);
saveas(hfig,fullfile(figPath,[figFileName '.png']),'png');
saveas(hfig,fullfile(figPath,[figFileName '.eps']),'eps2');
saveas(hfig,fullfile(figPath,[figFileName '.pdf']),'pdf');





lfp = Trial.lfp.copy;
lfp.load(Trial,[8,35]);

specParms = struct('nFFT',2^9,...
                    'Fs',lfp.sampleRate,...
                    'WinLength',2^8,...
                    'nOverlap',2^8*.875,...
                    'FreqRange',[20,150]);

[ys,fs,ts] = fet_spec(Trial,lfp,'mtchglong',true,lfp.sampleRate,specParms);



%%Testing Grounds: Enter at Own Risk

Stc = Trial.load('stc','hand_labeled_rev1');

xyz = Trial.load('xyz').filter(gtwin(.25,Trial.xyz.sampleRate));


xyz.addMarker('fhcom',[.7,1,.7],...
              {{'head_back','head_front',[0,0,1]}},...
              xyz.com(xyz.model.rb({'head_left','head_front','head_right'})));

xyz.addMarker('fbcom',[.7,0,.7],...
              {{'spine_lower','pelvis_root',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root'})));


xyz.addMarker('fscom',[.7,0,.7],...
              {{'spine_middle','spine_upper',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_middle','spine_upper'})));

fps = xyz.vel({'spine_lower','spine_middle','head_front'},[1,2]);
fps = xyz.vel({'fbcom','fscom','fhcom'},[1,2]);
fps.data(fps.data<.0001) = .0001;
fps.data = log10(fps.data);

[ys,fs,ts,ps] = fet_spec(Trial,fps,'mtchglong',true);



figure,imagesc(ts,fs,ys(:,:,1,7)');axis xy,
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k',[],3);


figure,imagesc(ts,fs,log10(ys(:,:,1,1))');axis xy,
Lines(Trial.stc{'w',ys.sampleRate}(:),[],'k',[],3);


ys(~nniz(ys(:)))=1e-10;
ind = nniz(ys(:,1,1,1));
ind = Stc{'r',ys}.cast('TimeSeries');
ind.data(1:10) = 0;

figure,hist2(sq([nanmean(log10(ys(ind(1:end-5)==1,fs<4,1,1)),2),...
                 nanmean(log10(ys(ind(1:end-5)==1,fs<4,3,3)),2)]),linspace(-8,.5,100),linspace(-7,.5,100))


figure,hist2([fps(:,1),fps(:,3)],linspace(-3,1.5,100),linspace(-3,1.5,100)),caxis([0,1000])

vind = Trial.stc{'w'};
figure,hist2([fps(vind,1),fps(vind,2)],linspace(-3,1.5,100),linspace(-3,1.5,100)),caxis([0,1000])



%% PPCA figure




TrialList = {'jg05-20120309','jg05-20120310','jg05-20120311','jg05-20120317'};
%TrialList = {'jg05-20120310','jg05-20120317'};
nt = numel(TrialList);

states = {'a','r','w'};
nsts = numel(states);

Trial = MTATrial(TrialList{1});
fet = fet_lgr(Trial);
Stc = Trial.stc.copy;
syncCount = diff(Trial.stc.sync([1,end]));
for t = 2:nt
    Trial = MTATrial(TrialList{t});
    Trial.stc.updateMode('auto_wbhr');
    Trial.stc.load;
    fet.data = cat(1,fet.data,...
                   subsref(fet_lgr(Trial),substruct('.',{'data'})));    
    for s = 1:nsts,
        sti = Stc.gsi(states{s});
        nstate = Trial.stc.states{Trial.stc.gsi(states{s})}.copy;
        nstate.resample(Stc.states{sti}.sampleRate);
        nstate.data = nstate.data+round(syncCount*Stc.states{sti}.sampleRate);
        Stc.states{sti}.data = cat(1,Stc{states{s}}.data,nstate.data);
    end
    syncCount = syncCount + diff(Trial.stc.sync([1,end]));
end



figure,
for t = 0:nt
    if t>0,
        Trial = MTATrial(TrialList{t});
        fet = fet_lgr(Trial);
        Stc = Trial.load('stc','auto_wbhr');
    end

    fet.data = nunity(fet.data);

    nind = nniz(fet);
    [U,S,V] = svd(cov(fet(nind,:)));

    subplot2(5,nt+1,1,t+1);
    hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],...
          linspace(-8,5,100),linspace(-6,4.5,100)),
    if t<1,caxis([0,3000]);else caxis([0,300]);end

    nind = Stc{'a-r'};
    subplot2(5,nt+1,2,t+1);
    hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],...
          linspace(-8,5,100),linspace(-6,4.5,100)),
    if t<1,caxis([0,3000]);else caxis([0,300]);end

    [U,S,V] = svd(cov(fet(nind,:)));
    subplot2(5,nt+1,3,t+1);
    hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],...
          linspace(-8,5,100),linspace(-6,4.5,100))
    if t<1,caxis([0,3000]);else caxis([0,300]);end

    nind = Stc{'a-r-w'};
    subplot2(5,nt+1,4,t+1);
    hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],...
          linspace(-8,5,100),linspace(-6,4.5,100))
    if t<1,caxis([0,3000]);else caxis([0,300]);end

    [U,S,V] = svd(cov(fet(nind,:)));
    subplot2(5,nt+1,5,t+1);
    hist2([fet(nind,:)*V(:,1),fet(nind,:)*V(:,2)],...
          linspace(-8,5,100),linspace(-6,4.5,100))
    if t<1,caxis([0,3000]);else caxis([0,300]);end

end

% $$$ figure
% $$$ per = Trial.stc{'w'};
% $$$ hist2([fet(per,:)*V(:,1),fet(per,:)*V(:,2)],linspace(-8,5,100),linspace(-6,4.5,100)),caxis([0,300])





 