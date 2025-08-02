function [pfs] = compute_bhv_ratemaps_crossval(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps_crossval(Trial,varargin);
% 
% tag [cv1,cv2];
%
% Compute spatially restricted behavior ratemaps
%
% CODE STRUCTURE
% 
%  1. ATTEMPT to load existing MTAApfs object from filesytem
%  2. IF all units exist return object
%  2. ELSE Compute
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                         [],                                            ...
                 'get_featureSet',                @fet_HB_pitchB,                                ...
                 'sampleRate',                    16,                                            ...
                 'pft',                           [],                                            ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [0.1,0.1],          ...
                                                         'SmoothingWeights', [1.8,1.8],          ...
                                                         'numIter',          1,                  ...
                                                         'boundaryLimits',   [-2,0.8;-0.8,2],    ...
                                                         'halfsample',       false),             ...
                 'threshRate',                    0.8,                                           ...
                 'threshDist',                    250,                                           ...
                 'tag',                           'cv1',                                         ...                 
                 'stcMode',                       'msnn_ppsvd_raux',                             ...
                 'states',                        {{'theta','rear','hloc','hpause',              ...
                                                   'lloc','lpause','groom','sit'}},              ...
                 'overwrite',                     false                                          ...
);
[units,get_featureSet,sampleRate,pft,pfsArgs,threshRate,threshDist,tag,stcMode,states,overwrite]=...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({tag,mfilename(),func2str(get_featureSet),sampleRate,threshRate,threshDist,    ...
                   pfsArgs.states,pfsArgs.binDims,pfsArgs.numIter,pfsArgs.boundaryLimits});
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------


% ATTEMPT to load existing MTAApfs object from filesytem
pfs = [];
if ~overwrite && MTAApfs.exist(Trial,pfsTag)
% LOAD existing MTAApfs object from filesytem    
    pfs = MTAApfs(Trial,'tag',pfsTag);
    if all(ismember(units,pfs.data.clu)),
% RETURN MTAApfs object        
        return;
    end
end



if isempty(pft),
% LOAD spatial placefields
    pft = pfs_2d_theta(Trial,units);
end

% LOAD Feature space
fet  = get_featureSet(Trial,sampleRate);
% LOAD distance metrics
xyz  = resample(preproc_xyz(Trial,'trb'),sampleRate);
drz  = compute_drz(Trial,units,pft,'feature',xyz);
ddz  = compute_ddz(Trial,units,pft,'feature',xyz);

% LOAD Spike clusters
spk  = create(copy(Trial.spk),Trial,sampleRate,'',units,'deburst');

% STCM : State Collection Matrix
stcm = stc2mat(Trial.load('stc',stcMode),xyz,states);

% APER : theta periods without sit or groom
FFLAG = '';
switch tag
  case 'cv1'
    FFLAG = 'last';
  case 'cv2'
    FFLAG = 'first';
end



% ASSIGN name and units to MTAApfs computational arguments
pfsArgs.tag   = pfsTag;
pfsArgs.units = units;
pfsArgs.spk   = spk;



% MAYBE overwrite all values
pfsArgs.overwrite = overwrite;        
    
u = 1;        
% RESTRICT periods by drz and ddz
pfsArgs.units  = units(u);

stcmSubset = stcm;
if ~isempty(spk.per),
    sper = spk.per.copy();
    sper.data = sper.data(spk.perInd(units(u),:),:);
    cast(sper,'TimeSeries',xyz,'absolute');
    resample(sper,xyz);
    stcmSubset(~sper.data) = 0;
end
for s = 2:6
    stcmSubset(find(stcmSubset(:,s),round(sum(nniz(stcmSubset(:,s)))/2),FFLAG),s) = 0;
end
aper = stcmSubset(:,1)==1 & any(stcmSubset(:,2:6),2);

pfsArgs.states = MTADepoch([],                                           ...
                           [],                                           ...
                           ThreshCross(aper & abs(drz(:,u))<threshRate   ...
                                            & ddz(:,u)<threshDist,       ...
                                       0.5,0),                           ...
                           fet.sampleRate,fet.sync.copy(),               ...
                           fet.origin,'TimePeriods','sts',[],'tdrz','d');

pfsArgs.autoSaveFlag = false;
pfsArgs.xyzp = fet;
pfsArgsArray = struct2varargin(pfsArgs);
pfs = Trial;
pfs = MTAApfs(pfs,pfsArgsArray{:});
pfs.save();
for u = 1:numel(units);
% RESTRICT periods by drz and ddz
    drzState{u} = MTADepoch([],                                           ...
                            [],                                           ...
                            ThreshCross(aper & abs(drz(:,u))<threshRate   ...
                                             & ddz(:,u)<threshDist,       ...
                                        0.5,0),                           ...
                            fet.sampleRate,fet.sync.copy(),               ...
                            fet.origin,'TimePeriods','sts',[],'tdrz','d');
    pfsArgs.units  = units(u);
    pfsArgs.states = drzState{u};
    pfsArgsArray = struct2varargin(pfsArgs);
    pfs = MTAApfs(pfs,pfsArgsArray{:});
end
pfs.save();

% $$$ if nargout==1,
% $$$     return;
% $$$ else
% $$$     auxData.periods     = {drzState};
% $$$     auxData.sampleRate  = sampleRate;
% $$$     auxData.get_featureSet = fun2str(get_featureSet);
% $$$     auxData.threshRate  = threshRate;
% $$$     auxData.threshDist  = threshDist;    
% $$$ end


% END MAIN -----------------------------------------------------------------------------------------