function [drzScore,drzCenter,pfds] = compute_hrz(Trial,varargin)
%function [drzScore,drzCenter] = compute_hrz(Trial,varargin)
%
%  
%
%  colapse 2d place field 
%
%  varargin:
%
%    units - NumericArray: contains integer identifiers of units
%
%    pft -   MTAApfs: MTA placefield object
%
%    pfstats - Struct: contains information on placefield statistics
%
%    filtCutOffFreq - Numeric: frequency cutoff for position filtering
%
%    marker - String: name of marker used if feature is empty
%
%    interpPar - Struct: contains interpolation parameters for higher resolution placefields
%
%    feature - NumericMatrix[Nx2]: trajectories to be compressed to 1 dimension
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                  [],                                                   ...
                 'pft',                    [],                                                   ...
                 'pfstats',                [],                                                   ...
                 'filtCutOffFreq',         2.5,                                                  ...
                 'marker',                 'nose',                                               ...
                 'interpPar',              struct('bins',{{linspace(-500,500,200)',              ...
                                                           linspace(-500,500,200)'}},            ...
                                                  'nanMaskThreshold', 0,                         ...
                                                  'methodNanMap',     'linear',                  ...
                                                  'methodRateMap',    'linear'),                 ...
                 'feature',                '',                                                   ...
                 'pfm',                    [],                                                   ...
                 'maxRateScaleFactor',     1,                                                    ...
                 'sampleRate',             []                                                    ...
);
[units,pft,pfstats,filtCutOffFreq,marker,interpPar,feature,pfm,maxRateScaleFactor,sampleRate] =  ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

if isempty(pfm),
    pfm = pft;
end

xyz = preproc_xyz(Trial,'trb');
if ~isempty(feature)
    xyz.resample(feature);
else    
    xyz.resample(sampleRate);
end
if isempty(feature),
% LOAD xyz as feature set
% SELECT xy plane of specified marker as feature    
    feature = preproc_xyz(Trial,'trb');
    feature.resample(sampleRate);    
    nind = nniz(feature);
    feature.filter('ButFilter',3,filtCutOffFreq,'low');
    feature.data = sq(feature(:,marker,[1,2]));
    feature.data(~nind,:) = 0; % BUG: patch
elseif isa(feature,'MTADxyz'),
    feature = copy(feature);    
    feature.resample(sampleRate);    
    nind = nniz(feature);
    feature.filter('ButFilter',3,filtCutOffFreq,'low');
    feature.data = sq(feature(:,marker,[1,2]));
    feature.data(~nind,:) = 0; % BUG: patch
else    
    feature.resample(sampleRate);    
    nind = nniz(feature);
end


% SET bins for estimation of 2D placefield centers
if isempty(interpPar),
    bins = pft.adata.bins;
else
    bins = interpPar.bins;
end


% CONSTRUCT map of expected rate given position
maxRate      = nan([numel(units),1]);
drzCenter    = nan([numel(units),2]);
wpmr = zeros(feature.size(1),numel(units));
for u = 1:numel(units),
    featureCA = mat2cell(feature(nind,:),sum(nind),ones([1,size(feature,2)]));
    rateMap = pfm.plot(units(u),'mean',false,[],false,0.99,false,[]);
    maxRate(u) = max(rateMap(:));
    rateMap(isnan(rateMap)) = 0;
    interpolatedFeatureCA = cell([1:size(feature,2)]);
    wpmr(nind,u) = interpn(pfm.adata.bins{:},rateMap,featureCA{:},'linear');
end
wpmr(wpmr<0) = 0;

for u = 1:numel(units),
    [~,mxp] = max(reshape(pft.plot(units(u),'mean',false,[],false,0.99,false,interpPar),[],1));
    mxp = Ind2Sub(cellfun(@numel,bins),mxp);
    drzCenter(u,:) = [bins{1}(mxp(:,1)), bins{2}(mxp(:,2))];
end

% COMPUTE the trajectory heading
pfds = zeros([size(feature,1),1]);
pfdd = zeros([size(feature,1),1]);
peakPatchRate = maxRate'.*maxRateScaleFactor;
for unit = units
    pfhxy = cat(2,permute(feature.data(nind,[1,2]),[1,3,2]),xyz(nind,'hcom',[1,2]));
    pfhxy = cat(2,pfhxy,permute(repmat(drzCenter(unit==units,:),[sum(nind),1]),[1,3,2]));
    pfhxy = MTADxyz([],[],pfhxy,feature.sampleRate);
% SUBSTRACT reference trajectory from second trajectory    
    cor = cell(1,2);
    [cor{:}] = cart2pol(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2));
    cor = cell2mat(cor);    
% SUBSTRACT placefield center from positions        
    por = cell(1,2);
    [por{:}] = cart2pol(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2));
    por = cell2mat(por);
% TRANSFORM from Cartesian to polar coordinates    
    pfds(nind,unit==units) = circ_dist(cor(:,1),por(:,1));
    pfdd(nind,unit==units) = por(:,2);
    
end
pfd = zeros(size(pfds));
pfd(abs(pfds(:))>=pi/2)=-1;
pfd(abs(pfds(:))<pi/2)=1;

% CALCULATE DRZ 
%drz = pfd.*(1-bsxfun(@rdivide,wpmr,mrt'));
drzScore = (pfd).*(1-bsxfun(@rdivide,wpmr,peakPatchRate));

% $$$ if MTA_DIAGNOSTIC_STATE
% $$$     figDir = fullfile(Session.spath,'figures',mfilename);
% $$$     create_directory(figDir);
% $$$     ind = [Trial.stc{'theta'}];
% $$$     for u = units,
% $$$         plot(Pfs,92);
% $$$         scatter(feature(ind,'nose',1),feature(ind,'nose',2),10,wpmr(ind,30));
% $$$         scatter(feature(ind,'nose',1),feature(ind,'nose',2),10,drz(ind,u=units));
% $$$     end
% $$$ end



