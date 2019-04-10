function [pfs] = compute_ratemaps_crossval(Trial,varargin);
% function [pfs,auxData] = compute_bhv_ratemaps(Trial,varargin);
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
                 'get_featureSet',                @fet_xy,                                       ...
                 'sampleRate',                    16,                                            ...
                 'pfsArgs',                       struct('states',           'theta-groom-sit',  ...
                                                         'binDims',          [20,20],            ...
                                                         'SmoothingWeights', [3.5,3.5],          ...
                                                         'numIter',          1,                  ...
                                                         'boundaryLimits',   [-500,500;-500,500],...
                                                         'halfsample',       false),             ...
                 'tag',                           'cv1',                                         ...   
                 'stcMode',                       'msnn_ppsvd_raux',                             ...
                 'states',                        {{'theta','rear','hloc','hpause',              ...
                                                   'lloc','lpause','groom','sit'}},              ...
                 'overwrite',                     false                                          ...
);
[units,get_featureSet,sampleRate,pfsArgs,tag,stcMode,states,overwrite] =                         ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% TAGS ---------------------------------------------------------------------------------------------
pfsTag = DataHash({tag,mfilename(),func2str(get_featureSet),sampleRate,...
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

% LOAD Spike clusters
spk  = create(copy(Trial.spk),Trial,sampleRate,'',units,'deburst');
% LOAD Feature space
fet  = get_featureSet(Trial,sampleRate);
% STCM : State Collection Matrix
stcm = stc2mat(Trial.load('stc',stcMode),fet,states);
% APER : theta periods without sit or groom
FFLAG = '';
switch tag
  case 'cv1'
    FFLAG = 'last';
  case 'cv2'
    FFLAG = 'first';
end
for s = 2:6
    stcm(find(stcm(:,s),round(sum(nniz(stcm(:,s)))/2),FFLAG),s) = 0;
end
aper = stcm(:,1)==1 & any(stcm(:,2:6),2);


% ASSIGN name and units to MTAApfs computational arguments
pfsArgs.tag   = pfsTag;
pfsArgs.units = units;

% MAYBE overwrite all values
pfsArgs.overwrite    = overwrite;
pfsArgs.units        = units;
pfsArgs.autoSaveFlag = true;
pfsArgs.xyzp         = fet;
pfsArgs.spk          = spk;
pfsArgs.states = MTADepoch([],                                           ...
                           [],                                           ...
                           ThreshCross(aper,0.5,0),                      ...
                           fet.sampleRate,fet.sync.copy(),               ...
                           fet.origin,'TimePeriods','sts',[],'tdrz','d');

pfsArgsArray = struct2varargin(pfsArgs);
pfs = MTAApfs(Trial,pfsArgsArray{:});

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