function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref_rev6(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});


Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,    ...
                 'normalize'    , false,                   ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});

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
              [],'TimeSeries',[],'body referenced position and motion',mfilename,'b');

% PREPROC xyz
xyz = preproc_xyz(Trial,procOpts);
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[0.1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);
clear('hcom');
xyz.filter('ButFilter',3,20,'low');


% Translational movements relative to body
shft = 3;
tmar = {'spine_lower','pelvis_root','spine_middle','spine_upper'};
tvec = zeros([size(xyz,1),numel(tmar),2]);
dzvec = zeros([size(xyz,1),numel(tmar),1]);
for m = 1:numel(tmar),
    tvec(:,m,1:2) = circshift(xyz(:,tmar{m},[1,2]),-shft)-circshift(xyz(:,tmar{m},[1,2]),shft);
    dzvec(:,m,1)  = circshift(xyz(:,tmar{m},[3]),  -shft)-circshift(xyz(:,tmar{m},[3]),  shft);
end
dzvec(~nniz(dzvec(:)))=0;

unvec = [];
rotationAngles = deg2rad([0,90]);
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2]);
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



fldwalkFetRot = MTADxyz('data',cat(2,dwalkFetRot,permute(dzvec,[1,3,2])),'sampleRate',xyz.sampleRate);
fldwalkFetRot.filter('ButFilter',3,2,'low');

fmdwalkFetRot = MTADxyz('data',cat(2,dwalkFetRot,permute(dzvec,[1,3,2])),'sampleRate',xyz.sampleRate);
fmdwalkFetRot.filter('ButFilter',3,[2,8],'bandpass');


tvec = zeros([size(xyz,1),numel(tmar),2]);
zvec = zeros([size(xyz,1),numel(tmar),1]);
for m = 1:numel(tmar),
    tvec(:,m,1:2) = xyz(:,tmar{m},[1,2])-xyz(:,'bcom',[1,2]);
    zvec(:,m,1)   = xyz(:,tmar{m},[3]);
end
zvec(~nniz(zvec(:)))=0;

unvec = [];
rotationAngles = deg2rad([0,90]);
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2]);
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

fmdwalkSegs = GetSegs([reshape(fmdwalkFetRot.data,size(xyz,1),[])], ...  x
                    1:size(dwalkFetRot,1),                         ...  start points
                    round(xyz.sampleRate/2),                     ...  segments' lengths
                    0                                            ...  If not complete
);

% CAT feature
fet.data = [ reshape(walkFetRot,size(xyz,1),[]),zvec,...
             reshape(fldwalkFetRot.data,size(xyz,1),[]),...
             permute(rms(fmdwalkSegs),[2,3,1]) ];
fet.data(~nniz(xyz),:)=0;

fet.filter('ButFilter',3,2.4,'low');

fet.resample(newSampleRate);

if normalize,
    [popMean,popStd] = load_normalization_parameters(mfilename);
    fet.unity(@nan,popMean,popStd);
end

    

featureTitles = {};
featureDesc = {};
if nargout>1,

end


%---------------------------------------------------------------------------------------------------