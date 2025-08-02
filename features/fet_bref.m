function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,                 ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQI'});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              Trial.xyz.sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'body referenced position and motion','fet_bref','b');                  

% PREPROC xyz
xyz = preproc_xyz(Trial,procOpts);
%xyz = Trial.load('xyz','trb');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[5]./(xyz.sampleRate/2),'low'));
% $$$ rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
% $$$ hcom = xyz.com(rb);
% $$$ xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
% $$$ clear('hcom');




% Tranlational movements relative to body
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
tvec = [];cvec = [];,zvec=[];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    cvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    dzvec(:,m,:) = circshift(xyz(:,tmar{m},[3]),-shft)-circshift(xyz(:,tmar{m},[3]),shft);
end
%zvec = nunity(zvec);
dzvec(~nniz(dzvec(:)))=0;

unvec = [];
rotationAngles = deg2rad([0,90]);
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2]);
%mvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

nind = nniz(tvec);
dwalkFetRot = zeros([size(nind,1),numel(rotationAngles),numel(tmar)]);
for t = rotationAngles;
    for m = 1:numel(tmar),
        dwalkFetRot(nind,t==rotationAngles,m) = dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3);
    end
end


shft = 0;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
tvec = [];cvec = [];,zvec=[];
for m = 1:numel(tmar),
    tvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,'bcom',[1,2]),shft);
    cvec(:,m,:) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,'bcom',[1,2]),shft);
    zvec(:,m,:) = xyz(:,tmar{m},[3]);
end
%zvec = nunity(zvec);
zvec(~nniz(zvec(:)))=0;

unvec = [];
rotationAngles = deg2rad([0,90]);
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2]);
%mvec = xyz(:,'fbcom',[1,2])-xyz(:,'fsl',[1,2]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

nind = nniz(tvec);
walkFetRot = zeros([size(nind,1),numel(rotationAngles),numel(tmar)]);
for t = rotationAngles;
    for m = 1:numel(tmar),
        walkFetRot(nind,t==rotationAngles,m) = dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3);
    end
end





% CAT feature
fet.data = [ reshape(walkFetRot,size(xyz,1),[]),zvec,reshape(dwalkFetRot,size(xyz,1),[]),dzvec ];
% $$$ defSpec = struct('nFFT',2^9,'Fs',fet.sampleRate,...
% $$$                  'WinLength',2^8,'nOverlap',2^8-4,...
% $$$                  'FreqRange',[1,15]);
% $$$ for s = 1:size(fet,2)
% $$$     tfet = fet.copy;
% $$$     tfet.data = tfet.data(:,s); 
% $$$     [sfet{s},fs,ts] = fet_spec(Trial,tfet,'mtchglong',false,'defspec',defSpec);
% $$$ end
% $$$ 
% $$$ sfet = cf(@(f) f.data,sfet);
% $$$ fet.data = cat(2,sfet{:});
% $$$ fet.data(:,1:2:end) = [];
% $$$ fet.sampleRate = 1./diff(ts(1:2));
% $$$ xyz.resample(fet);

fet.data(~nniz(xyz),:)=0;

fet.resample(newSampleRate);

featureTitles = {};
featureDesc = {};
if nargout>1,

end


%---------------------------------------------------------------------------------------------------