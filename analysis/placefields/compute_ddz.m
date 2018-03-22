function [ddz,pfds] = compute_ddz(Trial,varargin)

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

%if isempty(pfstats),    pfstats = compute_pfstats_bs(Trial);    end

dims = 2;

%[mrt,mrp] = pft.maxRate(units);

maxRatePos = nan([numel(units),2]);

if isempty(interpPar),
    bins = pft.adata.bins;
else
    bins = interpPar.bins;
end

for u = 1:numel(units),
    rateMap = pft.plot(units(u),'mean',false,[],false,0.99,false,interpPar);
    [~,mxp]  = max(rateMap(:));
    mxp = Ind2Sub(cellfun(@numel,bins),mxp);
    maxRatePos(u,:) = [bins{1}(mxp(:,1)), ...
                       bins{2}(mxp(:,2))];
end

xyz = preproc_xyz(Trial,'trb');
xyz.filter('ButFilter',3,filtCutOffFreq,'low');


% Get the rat's heading 
pfds = [];
pfdd = [];
for unit = units
    pfhxy = cat(2,xyz(:,{'nose'},:),circshift(xyz(:,{'nose'},:),round(xyz.sampleRate/5)));
    %pfhxy = xyz(:,{'head_back','head_front'},:);
    pfhxy = cat(2,pfhxy,permute(repmat([maxRatePos(unit==units,:),0],[size(xyz,1),1]),[1,3,2]));
    %pfhxy = cat(2,pfhxy,permute(repmat([fliplr(sq(mean(pfstats.peakPatchCOM(8,:,pfstats.cluMap==unit,:)))'),0],[size(xyz,1),1]),[1,3,2]));    

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
ddz = pfd.*pfdd;


