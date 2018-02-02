function drz = compute_drz(Trial,varargin)
%function drz = compute_drz(Trial,varargin)
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
                 'marker',                 'nose'                                                ...
);
[units,pft,pfstats,filtCutOffFreq,marker] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

%[mrt,mrp] = pft.maxRate(units,'mode','first');
[mrt,mrp] = pft.maxRate(units);
xyz = preproc_xyz(Trial,'trb');
xyz.filter('ButFilter',3,filtCutOffFreq,'low');

% Get the mean firing rate for each xy position along trajectory 
wpmr = zeros(xyz.size(1),numel(units));
[~,indx] = min(abs( repmat(pft.adata.bins{1}',xyz.size(1),1)...
                    -repmat(xyz(:,marker,1),1,numel(pft.adata.bins{1}))),...
               [],2);
[~,indy] = min(abs( repmat(pft.adata.bins{2}',xyz.size(1),1)...
                    -repmat(xyz(:,marker,2),1,numel(pft.adata.bins{2}))),...
               [],2);

rateMapIndex = sub2ind(pft.adata.binSizes',indx,indy);
for unit = units,
    %rateMap = pft.plot(unit,1);
    rateMap = pft.plot(unit,'mean');
    wpmr(:,unit==units) = rateMap(rateMapIndex);
end

% Get the rat's heading 
pfds = [];
pfdd = [];
peakPatchRate = [];
peakPatchRate = mrt';
for unit = units
    %peakPatchRate(1,end+1) = mean(pfstats.peakPatchRate(8,:,pfstats.cluMap==unit),'omitnan');
    %pfhxy = xyz(:,{'spine_middle','head_back'},:);
    pfhxy = cat(2,xyz(:,marker,:),circshift(xyz(:,marker,:),round(xyz.sampleRate/5)));
    pfhxy = cat(2,pfhxy,permute(repmat([mrp(unit==units,:),0],[size(xyz,1),1]),[1,3,2]));
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
drz = pfd.*(1-bsxfun(@rdivide,wpmr,peakPatchRate));


