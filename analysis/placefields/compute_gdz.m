function [gdz,fieldCenter,pfds] = compute_gdz(Trial,varargin)
% function [ddz,fieldCenter,pfds] = compute_gdz(Trial,varargin)
%
%  place field center distance 
%
%  colapse 2d place field using the trajectory heading to asign toward and away from field as 
%  -1 and +1 resepectively.
%
%
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
                 'marker',                 'hcom',                                               ...
                 'interpPar',              struct('bins',{{linspace(-500,500,200)',              ...
                                                           linspace(-500,500,200)'}},            ...
                                                  'nanMaskThreshold', 0,                         ...
                                                  'methodNanMap',     'linear',                  ...
                                                  'methodRateMap',    'linear'),                 ...
                 'feature',                '',                                                   ...
                 'sampleRate',             [],                                                   ...
                 'sigma',                  [150]                                                 ...
);
[units,pft,pfstats,filtCutOffFreq,marker,interpPar,feature,sampleRate,sigma] =                   ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

% SET bins for estimation of placefield centers
if isempty(interpPar),
    bins = pft.adata.bins;
else,
    bins = interpPar.bins;
end

if isempty(feature),
% LOAD xyz as feature set
% SELECT xy plane of specified marker as feature
    feature = preproc_xyz(Trial,'trb');
    feature.resample(sampleRate);            
    nind = nniz(feature);    
    feature.filter('ButFilter',3,filtCutOffFreq,'low');
    feature.data = sq(feature(:,marker,[1,2]));
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


% ESTIMATE placefield centers
fieldCenter = nan([numel(units),2]);
for u = 1:numel(units),
    rateMap = pft.plot(units(u),'mean',false,[],false,0.99,false,interpPar);
    [~,mxp]  = max(rateMap(:));
    mxp = Ind2Sub(cellfun(@numel,bins),mxp);
    fieldCenter(u,:) = [bins{1}(mxp(:,1)), bins{2}(mxp(:,2))];
end

% COMPUTE the trajectory heading
pfds = [];
pfdd = [];
for unit = units
% CONCATENATE features and placefield center
    pfhxy = cat(2,permute(feature.data,[1,3,2]),circshift(permute(feature.data,[1,3,2]),round(feature.sampleRate/5)));
    pfhxy = cat(2,pfhxy,permute(repmat(fieldCenter(unit==units,:),[size(feature,1),1]),[1,3,2]));
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
    pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
    pfdd(:,unit==units) = por(:,2);    
end

% MAP trajectory heading 
%    -1 -> towards the placefield center
%     1 -> outwards from placefield center
pfd = zeros(size(pfds));
pfd(abs(pfds(:))>=pi/2)=-1;
pfd(abs(pfds(:))<pi/2)=1;
% CALCULATE DDZ 
gdz = 1./(sigma.*sqrt(2.*pi)).*exp(-0.5.*(pfdd./sigma).^2);
gdz = pfd.*(1-bsxfun(@rdivide,gdz,max(gdz)));


