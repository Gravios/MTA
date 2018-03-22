function [drzScore,drzCenter] = compute_drz(Trial,varargin)
%function [drzScore,drzCenter] = compute_drz(Trial,varargin)
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
                                                  'methodRateMap',    'linear')                  ...
);
[units,pft,pfstats,filtCutOffFreq,marker,interpPar] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%global MTA_DIAGNOSTICS_PATH

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

%[mrt,mrp] = pft.maxRate(units,'mode','first');
%[mrt,mrp] = pft.maxRate(units);
xyz = preproc_xyz(Trial,'trb');
xyz.filter('ButFilter',3,filtCutOffFreq,'low');


if isempty(interpPar),
    bins = pft.adata.bins;
else
    bins = interpPar.bins;
end


% Get the mean firing rate for each xy position along trajectory 
wpmr = zeros(xyz.size(1),numel(units));
[~,indx] = min(abs( repmat(bins{1}',xyz.size(1),1)...
                    -repmat(xyz(:,marker,1),1,numel(bins{1}))),...
               [],2);
[~,indy] = min(abs( repmat(bins{2}',xyz.size(1),1)...
                    -repmat(xyz(:,marker,2),1,numel(bins{2}))),...
               [],2);

rateMapIndex = sub2ind(cellfun(@numel,bins),indx,indy);
maxRate = [];
for u = 1:numel(units),
    rateMap = pft.plot(units(u),'mean',false,[],false,0.99,false,interpPar);
    wpmr(:,u) = rateMap(rateMapIndex);

    [maxRate(u),mxp]  = max(rateMap(:));
    mxp = Ind2Sub(cellfun(@numel,bins),mxp);
    drzCenter(u,:) = [bins{1}(mxp(:,1)), ...
                       bins{2}(mxp(:,2))];
end

% Get the rat's heading 
pfds = [];
pfdd = [];
peakPatchRate = [];
peakPatchRate = maxRate;
for unit = units
    %peakPatchRate(1,end+1) = mean(pfstats.peakPatchRate(8,:,pfstats.cluMap==unit),'omitnan');
    %pfhxy = xyz(:,{'spine_middle','head_back'},:);
    pfhxy = cat(2,xyz(:,marker,:),circshift(xyz(:,marker,:),round(xyz.sampleRate/5)));
    pfhxy = cat(2,pfhxy,permute(repmat([drzCenter(unit==units,:),0],[size(xyz,1),1]),[1,3,2]));
    %pfhxy = cat(2,pfhxy,permute(repmat([fliplr(sq(mean(pfstats.peakPatchCOM(8,:,pfstats.cluMap==unit,:),'omitnan'))'),0],[size(xyz,1),1]),[1,3,2]));
    pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
    
    cor = cell(1,3);
    [cor{:}] = cart2pol(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
    cor = cell2mat(cor);
    
    por = cell(1,3);
    [por{:}] = cart2pol(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    por = cell2mat(por);
    
    pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
    pfdd(:,unit==units) = por(:,2);
    
end
pfd = zeros(size(pfds));
pfd(abs(pfds(:))>=pi/2)=-1;
pfd(abs(pfds(:))<pi/2)=1;

% Calculate DRZ 
%drz = pfd.*(1-bsxfun(@rdivide,wpmr,mrt'));
drzScore = pfd.*(1-bsxfun(@rdivide,wpmr,peakPatchRate));

% $$$ if MTA_DIAGNOSTIC_STATE
% $$$     figDir = fullfile(Session.spath,'figures',mfilename);
% $$$     create_directory(figDir);
% $$$     ind = [Trial.stc{'theta'}];
% $$$     for u = units,
% $$$         plot(Pfs,92);
% $$$         scatter(xyz(ind,'nose',1),xyz(ind,'nose',2),10,wpmr(ind,30));
% $$$         scatter(xyz(ind,'nose',1),xyz(ind,'nose',2),10,drz(ind,u=units));
% $$$     end
% $$$ end



