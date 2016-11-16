function [units,gunits,trajFeatures] = select_units_with_pfs_center_traversal(Trial,varargin);
%function units = select_units_with_pfs_center_traversal(Trial,varargin);
%
% Select units for which the subject traversed some minimum number
% of times and within some distance from each unit's place field center.
%
% varargin:
%   state:     string,  def('theta-sit-groom')
%              The state which forms the basis of the place field 
%
%   minTrajCount:  integer, def(5)
%   minDist:       numeric, def(200)     units(mm)
%   unitType:      string,  def('pyr')   enum('pyr','int')
%   tag:           string,  def('')      
%   overwrite:     logical, def(false)   
%
% assumptions: 
%   Trial.stc.mode is your standard labeling method and must be the same
%   mode used in other analyses when calling this function.
%
%

% CHECK that Trial is an MTATrial object and load if not
Trial= MTATrial.validate(Trial);    

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',         [],                                                            ...
                 'state',         'theta-sit-groom',                                             ...
                 'minTrajCount',  5,                                                             ...
                 'minDist',       200,                                                           ...
                 'minOcc',        round(Trial.xyz.sampleRate/2),                                 ...
                 'unitType',      'pyr',                                                         ...
                 'tag',           '',                                                            ...
                 'overwrite',     false                                                          ...
);%-------------------------------------------------------------------------------------------------

[units,state,minTrajCount,minDist,minOcc,unitType,tag,overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');


% SELECT units
if isempty(units)
    pft = pfs_2d_theta(Trial);
    mrt = pft.maxRate;
    units = select_units(Trial,18,unitType);
    units = units(mrt(pft.data.clu(units))>1);
end


% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
%
%    stcMode 
%    states 
if isempty(tag),
    tag = DataHash(struct('stcMode',Trial.stc.mode,    ...
                          'states', state,             ...
                          'minTrajCount',minTrajCount, ...
                          'minDist',minDist,           ...
                          'minOcc',minOcc,           ...
                          'unitType',unitType          ...
                          )                            ...
                   );
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

filename = fullfile(Trial.spath,[Trial.filebase,'_selUnitsPfcTraj_',tag,'.mat']);

if ~exist(filename,'file') || overwrite,

    % LOAD place fields
    defargs = get_default_args('MjgEdER2016','MTAApfs','struct');
    defargs.units = units;
    defargs.states = state;
    defargs.numIter = 1;
    defargs = struct2varargin(defargs);        
    pfs = MTAApfs(Trial,defargs{:});      
    
    % COLLECT place field features
    clear('pfstats')
    pfstats = {};
    for unit = units,
        pfstats(unit==units) = {PlaceFieldStats(Trial,pfs,unit)};
    end
    pfstats = CatStruct([pfstats{:}],[],1);
        
    % Detect trajectories for each place field center
    % pfstats.
    %               peakFR: [33x1 double]
    %        rateThreshold: [33x1 double]
    %     spatialCoherence: [33x1 double]
    %            patchArea: [33x2 double]
    %             patchCOM: [33x2x2 double]
    %             patchPFR: [33x2 double]
    %             patchMFR: [33x2 double]
    %             patchCnt: [33x2 double]
    %         patchRateInd: [33x2x2x300]
    %         patchRateMap: [33x2x300 double]

    % LOAD marker (rat) positions
    xyz = Trial.load('xyz');
    mid = xyz.model.gmi(Trial.trackingMarker);

    trajParameters.minDist = minDist;

    clear('trajFeatures')
    trajFeatures.pathPer = [];
    trajFeatures.ucDist = [];    
    trajFeatures.pathDist = [];        
    trajFeatures.pathSpeed = [];
    trajFeatures.isFullTrav = [];    
    trajFeatures.isHalfComDist = [];
    for unit = units,

        for p = 1:size(pfstats.patchCOM,2)
            % IGNORE patch if com isnan
            if any(isnan(pfstats.patchCOM(unit==units,p,2))),continue, end
            % COMPUTE distance of xy trajectiories to the place field's COM        
            patchComDist = sqrt(sum(bsxfun(@minus,sq(pfstats.patchCOM(unit==units,p,2)), ...
                                           sq(xyz(:,mid,[1,2]))).^2,2));
            % FIND time periods where rat's trajectory is with some minimum distance
            tper = ThreshCross(-patchComDist,-minDist,minOcc);

            for t = 1:size(tper,1)
                % Point to point distance on Threshold circumference
                trajFeatures(unit==units,p).ucDist(t,1) = sqrt(sum(diff(sq(xyz(tper(t,:)',mid,[1,2]))).^2,2));
                % 
                trajFeatures(unit==units,p).isHalfComDist(t,1) = any(patchComDist(tper(t,1):tper(t,2))<minDist/2);
                trajFeatures(unit==units,p).pathDist(t,1) = sum(sqrt(sum(diff(sq(xyz(tper(t,:),mid,[1,2]))).^2,2)));
                trajFeatures(unit==units,p).pathSpeed(t,1) = trajFeatures(unit==units,p).pathDist(t,1)./diff(tper(t,:)).*xyz.sampleRate/10;
            end
            
            trajFeatures(unit==units,p).pathPer = tper;
            trajFeatures(unit==units,p).isFullTrav = trajFeatures(unit==units,p).ucDist>sqrt(2)*minDist;
        end
    end
    
    gunits = units(any(cell2mat(arrayfun(@(x) sum(x.isFullTrav&x.isHalfComDist)>minTrajCount,...
                                        trajFeatures,...
                                        'UniformOutput',false)...
                               ),...
                      2)...
                  );
    
    save(filename,'units','gunits','trajFeatures');
else
    load(filename);
end


