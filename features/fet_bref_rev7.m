function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref_rev7(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {{}});


Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', Trial.xyz.sampleRate,    ...
                 'normalize'    , false,                   ...
                 'procOpts'     , {{'SPLINE_SPINE_HEAD_EQD'}});

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

xyz.filter('ButFilter',3,2.5,'low');



% Translational movements relative to body
shft = 3;
tmar = {'spine_lower','spine_middle','spine_upper','head_back','head_front'};
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


fldwalkFetRot = MTADxyz('data',cat(2,dwalkFetRot,permute(dzvec,[1,3,2])),'sampleRate',xyz.sampleRate);
fldwalkFetRot.filter('ButFilter',5,1.5,'low');


% CAT feature
fet.data = [ reshape(walkFetRot,size(xyz,1),[]),zvec,...
             reshape(fldwalkFetRot.data,size(xyz,1),[]),...
             reshape(fldwalkFetRot.data,size(xyz,1),[]).^2];
fet.data(~nniz(xyz),:)=0;


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