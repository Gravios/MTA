%MjgEd2016_supfig_SlowRespiration
%
% Exploration for slow timescale features coherent with respiration
% as detected by nasal cavity pressure sensor (ncp).
%

% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/ownCloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Suplementary';
% --------------------------------------------------------------------------------------


Trial = MTATrial.validate('Ed01-20140707');
Trial = MTATrial.validate('Ed01-20140709');
Stc = Trial.load('stc','NN0317');
embeddingWindow = 64;


% LOAD xyz data
xyz = Trial.load('xyz');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle',});
bcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(bcom,4,[.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},bcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsm',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_middle',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsu',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_upper',:),4,[1.5]./(xyz.sampleRate/2),'low'));

rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(xyz.sampleRate/2),'low'));
clear('hcom');
clear('bcom');

ang = create(MTADang,Trial,xyz);
ang.data(~nniz(xyz(:,1,1)),:,:,:) = 0;
%fxyz = xyz.copy();
%fxyz.filter('ButFilter',3,2.4);
%vxy = fxyz.vel([],[1,2]);
%vxy.data(vxy.data<=1e-3)=1e-3;

% traj NORM body vector projection
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
tvec = [];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
end
nind = nniz(tvec);

% body vector
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'fsl',[1,2]);
umvec = bsxfun(@rdivide,bsxfun(@times,permute([1,-1],[1,3,2]),mvec(:,1,[2,1])),sqrt(sum(mvec.^2,3)));

unvec = [];
rotationAngles = deg2rad([0,45,90,135]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

walkFetRot = [];
for t = rotationAngles;
    for m = 1:numel(tmar),
        walkFetRot(nind,t==rotationAngles,m) = nunity(dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3));
    end
end


% SVD on XY translational vectors
nz = nniz(xyz);
wfet = xyz.copy;
wfet.data= zeros([size(xyz,1),size(walkFetRot,2)*size(walkFetRot,3)]);
wfet.data(nz,:) = [reshape(walkFetRot(nz,:),[],size(walkFetRot,2)*size(walkFetRot,3))];
wfs = wfet.segs([],embeddingWindow);
wfs = circshift(wfs,embeddingWindow/2,2);
wfs = MTADxyz('data',reshape(permute(wfs,[2,1,3]),size(wfs,2),[]),'sampleRate',xyz.sampleRate);
wfs.data(isnan(wfs.data(:)))=0;
[Uw,Sw,Vw] = svd(wfs(Stc{'w+n'},:),0);
wts = (1:embeddingWindow)./wfet.sampleRate;


% COMPUTE timeseries score for first 10 eigenvectors
fetW = MTADxyz('data',wfs.data*Vw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,fetW.data(:,i) = wfs.data*Vw(:,i);end

% REDUCED pc 
rVw = Vw;
for i = 1:16,
    rVw(i:embeddingWindow:end,:) = 0;
end
for i = 48:64,
    rVw(i:embeddingWindow:end,:) = 0;
end

rfetW = MTADxyz('data',wfs.data*rVw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:10,rfetW.data(:,i) = wfs.data*rVw(:,i);end

fetW.sync = xyz.sync.copy;
fetW.origin = xyz.origin;
rfetW.sync = xyz.sync.copy;
rfetW.origin = xyz.origin;




hang = xyz.copy;
hang.data = ButFilter(ang(:,'fbcom','spine_upper',3),4,[0.5,20]./(xyz.sampleRate/2),'bandpass');
hang.data = ButFilter(ang(:,'fsm','spine_upper',3),4,[0.5,20]./(xyz.sampleRate/2),'bandpass');


rhmOri = fet_rhm(Trial,[],'mta');
rhmExp = fet_rhm_exp(Trial,[],'mta');
ncp = fet_ncp(Trial,[],'mta');
%lfp = fet_ncp(Trial,[],'mta',[65:2:80]);

swayPhase = fetW.phase([1,5]);
rhmOriPhase = rhmOri.phase([0.5,15]);
rhmExpPhase = rhmExp.phase([0.5,15]);
%rhmPhase = rhm.phase([5,12]);
%lfpPhase = lfp.phase([6,12]);
ncpPhase = ncp.phase([0.5,15]);

hswPhase = hang.phase([.5,15]);
%bswPhase = bang.phase([.5,15]);


% PRINT timeseries of RBM and NCP signals with phase
states = {'walk','turn','pause','rear','sit'};
sclr = 'bgcrm';
ts = [1:size(xyz,1)]'./xyz.sampleRate;
sp = [];
hfig = figure(20170322);clf
set(hfig,'Units','centimeters',...
         'PaperPositionMode', 'auto')
sp(end+1)=subplot2(10,4,[1:8],[1:4]);
hold on
plot(ts,nunity(ncp.data)+6,'b')
plot(ts,ncpPhase.data+3,'k')
plot(ts,nunity(hang.data)*10-3,'m')
plot(ts,hswPhase.data,'c')
plot(ts,rhmOriPhase.data+11)
plot(ts,nunity(rhmOri.data)*20+15)
legend({'ncp','ncp phase','rbm','rbm phase','rhm','rhm phase'})
sp(end+1)=subplot2(10,4,[9,10],[1:4]);
plotSTC(Stc,1,'text',states,sclr,[],false);
linkaxes(sp,'x');
FigName = ['SlowRespiration_NCPxRBM_',Trial.filebase];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



% PRINT expanded RBMxNCP
hfig = bhv_rhm_ncp_distrb(Trial,[],[],[],'NN0317');
FigName = ['SlowRespiration_NCPxRBM_coherence',Trial.filebase];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



ind = Trial.stc{'s'};
figure
hist2([ncpPhase(ind,1),circshift(hswPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))



ind = Stc{'w'};
figure
%hist2([swayPhase(ind,4),circshift(hswPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
%hist2([swayPhase(ind,4),circshift(rhmExpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
hist2([swayPhase(ind,4),circshift(rhmOriPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))


% NCP 
figure,
subplot(131);
ind = Stc{'w'};
hist2([swayPhase(ind,4),circshift(rhmPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
subplot(132);
hist2([swayPhase(ind,4),circshift(ncpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
subplot(133);
hist2([rhmPhase(ind,1),circshift(ncpPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))


% LFP 
chan = 3;
shift = 0;
figure,
for chan = 1:8,
    %subplot2(8,3,chan,1);
ind = Stc{'w'};
hist2([swayPhase(ind,3),circshift(rhmPhase(ind),0)],linspace(-pi,pi,30),linspace(-pi,pi,30))
subplot2(8,3,chan,2);
%hist2([swayPhase(ind,3),circshift(lfpPhase(ind,chan),shift)],linspace(-pi,pi,30),linspace(-pi,pi,30))
%subplot2(8,3,chan,3);
hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift)],linspace(-pi,pi,30),linspace(-pi,pi,30))
end



% LFP 
chan = 3;
shift = -180:10:-20;
shift = -150:30:150;

ind = Stc{'w'};
figure,
for chan = 1:8,
    for shiftInd = 1:numel(shift)
        subplot2(8,numel(shift),chan,shiftInd);
%hist2([swayPhase(ind,3),circshift(rhmPhase(ind),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
%hist2([swayPhase(ind,3),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
    end
end

figure,chan = 5;shiftInd = 6;
%hist2([swayPhase(ind,4),circshift(rhmPhase(ind),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30))
hist2([rhmPhase(ind,1),circshift(lfpPhase(ind,chan),shift(shiftInd))],linspace(-pi,pi,30),linspace(-pi,pi,30));





rhmlpf = rhm.copy;
rhmlpf.filter('ButFilter',3,5,'low');

ind = Stc{'w'};
figure,plot(rhmlpf(ind),fetW(ind,3),'.');

figure,plot(nunity(rhmlpf(:)))
hold on,plot(fetW(:,3))

