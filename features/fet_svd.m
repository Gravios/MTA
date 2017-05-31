function [fet] = fet_svd(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_svd(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});


defargs = struct('newSampleRate', 12,                       ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');


% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_mis','m');                  

xyz = Trial.load('xyz','trb');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1.5]./(xyz.sampleRate/2),'low'));
clear('hcom');
ang = create(MTADang,Trial,xyz);



% Tranlational movements relative to body
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
tvec = [];cvec = [];,zvec=[];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    cvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    zvec(:,m,:) = circshift(xyz(:,tmar{m},[3]),-shft)-circshift(xyz(:,tmar{m},[3]),shft);
end
zvec = nunity(zvec);
zvec(~nniz(zvec(:)))=0;

unvec = [];
rotationAngles = deg2rad([0,45,90,135]);
mvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

nind = nniz(tvec);
for t = rotationAngles;
    for m = 1:numel(tmar),
        walkFetRot(nind,t==rotationAngles,m) = dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3);
    end
end



embeddingWindow = round(0.5*xyz.sampleRate);
embeddingWindow = embeddingWindow+mod(embeddingWindow,2);

%% Walk
% compute walk features (no shifts needed)
nz = nniz(xyz);
% svd
bhvPeriods = Trial.stc{'w+n'};
wfet = xyz.copy;
wfet.data= zeros([size(xyz,1),size(walkFetRot,2)*size(walkFetRot,3)]);
wfet.data(nz,:) = [reshape(walkFetRot(nz,:),[],size(walkFetRot,2)*size(walkFetRot,3)),zvec(nz,:)];
wfs = wfet.segs([],embeddingWindow);
wfs = circshift(wfs,embeddingWindow/2,2);
wfs = MTADxyz('data',reshape(permute(wfs,[2,1,3]),size(wfs,2),[]),'sampleRate',xyz.sampleRate);
wfs.data(isnan(wfs.data(:))) = 0;
[Uw,Sw,Vw] = svd(wfs(bhvPeriods,:),0);

%figure,for i = 1:25,subplot(5,5,i);imagesc(reshape(Vw(:,i),[],size(wfet,2))'),end

fetW = MTADxyz('data',wfs.data*Vw(:,1),'sampleRate',xyz.sampleRate);
for i = 1:25,fetW.data(:,i) = wfs.data*Vw(:,i);end

%% Rear
% compute rear features (shift required)
bhvPeriods = Trial.stc{'r+m'};
rfet = xyz.copy;
rfet.data= zeros([size(xyz,1),numel(tmar)]);
rfet.data(nz,:) = xyz(nz,tmar,3);
rfs = rfet.segs([],embeddingWindow);
rfs = circshift(rfs,embeddingWindow/2,2);
rfs = MTADxyz('data',reshape(permute(rfs,[2,1,3]),size(rfs,2),[]),'sampleRate',xyz.sampleRate);
[Ur,Sr,Vr] = svd(rfs(bhvPeriods,:),0);

%figure,for i = 1:25,subplot(5,5,i);imagesc(reshape(Vr(:,i),[],size(rfet,2))'),end

fetR = MTADxyz('data',rfs.data*Vr(:,1),'sampleRate',xyz.sampleRate);
for i = 1:15,fetR.data(:,i) = rfs.data*Vr(:,i);end


fet.data = [fetW.data,fetR.data];




% $$$ 
% $$$ %% Groom
% $$$ bhvPeriods = Trial.stc{'m'};
% $$$ mfet = xyz.copy;
% $$$ mfet.data= [];
% $$$ for i = 1:numel(tmar),
% $$$     for j = 1:i+1,
% $$$         for k = 1:j+1
% $$$             if k>j && k<=numel(tmar),
% $$$                 mfet.data(nz,end+1) = nunity(circ_dist(ang(nz,tmar{i},tmar{j},1),ang(nz,tmar{i},tmar{k},1)));
% $$$             end
% $$$         end
% $$$     end
% $$$ end
% $$$ for i = 1:numel(tmar),
% $$$     for j = 1:i+1,
% $$$         if j>i && j<=numel(tmar),
% $$$             mfet.data(nz,end+1) = nunity(ang(nz,tmar{i},tmar{j},2));
% $$$         end
% $$$     end
% $$$ end
% $$$ mfs = mfet.segs([],embeddingWindow);
% $$$ mfs = circshift(mfs,embeddingWindow/2,2);
% $$$ mfs = MTADxyz('data',reshape(permute(mfs,[2,1,3]),size(mfs,2),[]),'sampleRate',xyz.sampleRate);
% $$$ bhvPeriods.cast('TimeSeries',xyz);
% $$$ bhvPeriods = bhvPeriods.data&nniz(mfs);
% $$$ [Um,Sm,Vm] = svd(mfs(bhvPeriods,:),0);
% $$$ 
% $$$ %figure,for i = 1:15,subplot(3,5,i);imagesc(reshape(Vm(:,i),[],size(mfet,2))'),end
% $$$ 
% $$$ fetM = MTADxyz('data',mfs.data*Vm(:,1),'sampleRate',xyz.sampleRate);
% $$$ for i = 1:10,fetM.data(:,i) = mfs.data*Vm(:,i);end
% $$$ 
% $$$ 
% $$$ figure,
% $$$ hold on,
% $$$ plot(fetM(Trial.stc{'a-m-s-w-n'},1),fetM(Trial.stc{'a-m-s-w-n'},2),'.m')
% $$$ plot(fetM(Trial.stc{'m'},1),fetM(Trial.stc{'m'},2),'.b')


% $$$ [xyz,ss] = preproc_xyz(Trial,procOpts);
% $$$ xyz.resample(newSampleRate);
% $$$ ss.resample(xyz);
