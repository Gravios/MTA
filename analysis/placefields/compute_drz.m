function [drzScore,drzCenter] = compute_drz(Trial,varargin)
%function [drzScore,drzCenter] = compute_drz(Trial,varargin)
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
                 'feature',                ''                                                    ...
);
[units,pft,pfstats,filtCutOffFreq,marker,interpPar,feature]=DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

if isempty(feature),
% LOAD xyz as feature set
% SELECT xy plane of specified marker as feature    
    feature = preproc_xyz(Trial,'trb');
    feature.filter('ButFilter',3,filtCutOffFreq,'low');
    feature.data = sq(feature(:,marker,[1,2]));
end


% SET bins for estimation of placefield centers
if isempty(interpPar),
    bins = pft.adata.bins;
else
    bins = interpPar.bins;
end


% FIND bin index given position
wpmr = zeros(feature.size(1),numel(units));
[~,indx] = min(abs( repmat(bins{1}(:)',feature.size(1),1)...
                    -repmat(feature(:,1),1,numel(bins{1}))),...
               [],2);
[~,indy] = min(abs( repmat(bins{2}(:)',feature.size(1),1)...
                    -repmat(feature(:,2),1,numel(bins{2}))),...
               [],2);

% CONSTRUCT map of expected rate given position
% ESTIMATE placefield centers
maxRate      = nan([numel(units),1]);
drzCenter    = nan([numel(units),2]);
rateMapIndex = sub2ind(cellfun(@numel,bins),indx,indy);
for u = 1:numel(units),
    rateMap = pft.plot(units(u),'mean',false,[],false,0.99,false,interpPar);
    wpmr(:,u) = rateMap(rateMapIndex);
    [maxRate(u),mxp]  = max(rateMap(:));
    mxp = Ind2Sub(cellfun(@numel,bins),mxp);
    drzCenter(u,:) = [bins{1}(mxp(:,1)), bins{2}(mxp(:,2))];
end

% COMPUTE the trajectory heading
pfds = [];
pfdd = [];
peakPatchRate = [];
peakPatchRate = maxRate';
for unit = units
    pfhxy = cat(2,permute(feature.data,[1,3,2]),circshift(permute(feature.data,[1,3,2]),round(feature.sampleRate/5)));
    pfhxy = cat(2,pfhxy,permute(repmat(drzCenter(unit==units,:),[size(feature,1),1]),[1,3,2]));
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
pfd = zeros(size(pfds));
pfd(abs(pfds(:))>=pi/2)=-1;
pfd(abs(pfds(:))<pi/2)=1;

% CALCULATE DRZ 
%drz = pfd.*(1-bsxfun(@rdivide,wpmr,mrt'));
drzScore = pfd.*(1-bsxfun(@rdivide,wpmr,peakPatchRate));

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



