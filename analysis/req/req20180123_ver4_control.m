function [pfd,LRG,VTG,VSrG,UID,GUNITS] = req20180123_ver4_control(varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Trials',                        {{}},                                          ...
                 'sessionListName',               'MjgER2016',                                   ...
                 'version',                       '7',                                           ...
                 'overwrite',                     false,                                         ...
                 'display',                       false,                                         ...
                 'figDir',  '/storage/gravio/figures/analysis/placefields_nonSpatialFeatures',   ...
                 'maxPfIndex',                    4                                              ...
);
[Trials,sessionListName,version,overwrite,display,figDir,maxPfIndex] =                           ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


create_directory(figDir); 

% LOAD session list
% LOAD Trials
% FILTER units for analysis
sessionList = get_session_list(sessionListName);

if isa(Trials,'MTATrial')
    Trials = {Trials};
end


if isempty(Trials),
    Trials  = af(@(S)  MTATrial.validate(S),   sessionList);
end

units = cf(@(T)  select_placefields(T),  Trials);

numTrials = numel(Trials);

pfd = {};
UID = {};
LRG = {};
VTG = {};
VSrG = {};
%% pfd_BPITCHxHPITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERPCA BPITCHxHPITCH -----------------------------------------------------------------------------
analDir = ['control-v',version];
create_directory(fullfile(figDir,analDir));


pfindex = 1;
for tind = 1:numTrials,
    Trial = Trials{tind}; 

% LOAD vars
    if display,
        pft = pfs_2d_theta(Trial);
        xyz = preproc_xyz(Trial,'trb');
        pch = fet_HB_pitchB(Trial);
        drz = compute_drz(Trial,units{tind},pft);%,pfstats);
        tper =[Trial.stc{'theta-groom-sit'}];
        tper.resample(xyz);
    end

% COMPUTE HPITCH x BPITCH | DRZ [-0.5,0.5]
    pargs.tag            = ['DRZxHBPITCHxBPITCH_v',version];
    if overwrite,  
        pargs = get_default_args('MjgER2016','MTAApfs','struct');        
        pargs.tag            = ['DRZxHBPITCHxBPITCH_v',version];        
        pargs.units   = units{tind};
        pargs.numIter = 1001;
        pargs.halfsample = true;
        pargs.boundaryLimits = [-2,2;-2,2];
        pargs.binDims        = [0.1,0.1];
        pargs.SmoothingWeights = [1.5,1.5];
        pargs.overwrite = true;    
        pargs.xyzp = MTADxyz('data',pch.data,'sampleRate',xyz.sampleRate);
        pargs.autoSaveFlag = false;
        drzState = {};

% COMPUTE first to initialize object
        u = 1;
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u}  = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs      = struct2varargin(pargs);
        pfTemp = MTAApfs(Trial,pfsArgs{:});
        pfTemp.save();
        for u = 2:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        end
        pfTemp.save();
    end
    pfd{tind,pfindex} = MTAApfs(Trial,'tag',pargs.tag);

% VISUALIZE HPITCH x BPITCH | DRZ [-0.5,0.5]
    if display,  req20180123_vis_HPITCHxBPITCH;  end    
end%for tind

if (nargout>1||display)&&numel(sessionList)==numel(Trials)
% COMPUTE erpPCA_HPITCHxBPITCH eigenvectors
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
%rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,pfindex));              % GET rate maps
rmaps = cat(2, rmaps{:});                                            % CONCATENATE rate maps
[smnd,sind] = sort(sum(isnan(rmaps),1));                             % SORT by total occupancy
selectedUnits = sind(smnd>1100&smnd<1450);                           % SELECT rate maps by occupancy
zrmaps = rmaps(:,selectedUnits);                                     % GET rate map subset
zdims = size(zrmaps);                                                % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                                           % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/3;                         % SELECT subspace {1/3 non-zero samples}
zrmaps = zrmaps(validDimsInds,:);                                    % REDUCE to selected subspace
zrmaps = 10.*(zrmaps>0)+randn(size(zrmaps)).*(zrmaps>0);

[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                                  % COMPUTE covariance-based PCA with
V = LR;                                                              %         Varimax rotation

UID{pfindex} = selectedUnits;
LRG{pfindex} = LR;
VTG{pfindex} = VT;
VSrG{pfindex} = FSr;
if display,  req20180123_vis_HPITCHxBPITCH_erpPCA;  end
end

% ERPPCA HPITCHxBPITCH END ---------------------------------------------------------------------------





%% erpPCA_HPITCHxBSPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERPPCA HPITCHxBSPEED START -------------------------------------------------------------------------
analDir = ['HBPITCHxBSPEED-v',version];
create_directory(fullfile(figDir,analDir));

pfindex = 2;

for tind = 1:numTrials,
    Trial = Trials{tind}; %15,16,17,18

    if display || overwrite,
        pft = pfs_2d_theta(Trial);        
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('RectFilter');
        pch = fet_HB_pitchB(Trial);
        vxy = xyz.vel(['spine_lower'],[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);
        fet = pch.copy('empty');
        fet.label = 'fet_VP';       
        fet.data = [pch(:,1),vxy(:)];
        drz = compute_drz(Trial,units{tind},pft);%,pfstats);
        tper =[Trial.stc{'theta-groom-sit-rear'}];
        tper.resample(xyz);
    end

% COMPUTE pfd HBPITCHxBSPEED
    pargs.tag            = ['DRZxHBPITCHxBSPEED_v',version];
    if overwrite,  
        pargs = get_default_args('MjgER2016','MTAApfs','struct');        
        pargs.tag              = ['DRZxHBPITCHxBSPEED_v',version];
        pargs.units            = units{tind};
        pargs.numIter          = 1001;
        pargs.halfsample       = true;
        pargs.overwrite        = true;
        pargs.boundaryLimits   = [-2,2;-2,2];
        pargs.binDims          = [0.1,0.1];
        pargs.SmoothingWeights = [1.5,1.5];
        pargs.xyzp = MTADxyz('data',fet.data,'sampleRate',xyz.sampleRate);
        pargs.autoSaveFlag = false;
        
        drzState = {};
        u = 1;
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(Trial,pfsArgs{:});
        pfTemp.save();
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        end
        pfTemp.save();        
    end
    pfd{tind,pfindex} = MTAApfs(Trial,'tag',pargs.tag);
    
% DISPLAY pfd HBPITCHxBSPEED
    if display,  req20180123_vis_HPITCHxBSPEED;  end

end%for tind



% COMPUTE erpPCA_HPITCHxBSPEED eigenvectors ------------------------------------------------------------------
% CREATE SELCETED VECTOR SUBSPACE
if (nargout>1||display)&&numel(sessionList)==numel(Trials)
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
%rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,pfindex));
rmaps = cat(2, rmaps{:});
[smnd,sind] = sort(sum(isnan(rmaps),1));
selectedUnits = sind(smnd>950&smnd<1450);
zrmaps = rmaps(:,selectedUnits);
zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/2;        % SELECT subspace {1/2 non-zero samples} 
zrmaps = zrmaps(validDimsInds,:);                   % REDUCE to selected subspace
zrmaps = 10.*(zrmaps>0)+randn(size(zrmaps)).*(zrmaps>0);

% DECOMPOSE 
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                 % COMPUTE covariance-based PCA with Varimax rotation
V = LR;

UID{pfindex} = selectedUnits;
LRG{pfindex} = LR;
VTG{pfindex} = VT;
VSrG{pfindex} = FSr;

if display, req20180123_vis_HPITCHxBSPEED_erpPCA;  end
end
% ERPPCA HPITCHxSPEED END -------------------------------------------------------------------------






%% erpPCA_BPITCHxBSPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERPPCA BPITCHxBSPEED START -------------------------------------------------------------------------
analDir = ['BPITCHxBSPEED-v',version];
create_directory(fullfile(figDir,analDir));

pfindex = 3;

for tind = 1:numTrials,
    Trial = Trials{tind}; %15,16,17,18    

    if display || overwrite,
        pft = pfs_2d_theta(Trial);
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('RectFilter');
        pch = fet_HB_pitchB(Trial);
        vxy = xyz.vel(['spine_lower'],[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);
        fet = pch.copy('empty');
        fet.label = 'fet_VP';       
        fet.data = [pch(:,2),vxy(:)];
        drz = compute_drz(Trial,units{tind},pft);%,pfstats);
        tper =[Trial.stc{'theta-groom-sit-rear'}];
        tper.resample(xyz);
    end

% COMPUTE pfd HBPITCHxBSPEED
    pargs.tag            = ['DRZxBPITCHxBSPEED_v',version];
    if overwrite,  
        pargs = get_default_args('MjgER2016','MTAApfs','struct');        
        pargs.tag              = ['DRZxBPITCHxBSPEED_v',version];
        pargs.units            = units{tind};
        pargs.numIter          = 1001;
        pargs.halfsample       = true;
        pargs.overwrite        = true;
        pargs.boundaryLimits   = [-2,2;-2,2];
        pargs.binDims          = [0.1,0.1];
        pargs.SmoothingWeights = [1.5,1.5];
        pargs.xyzp = MTADxyz('data',fet.data,'sampleRate',xyz.sampleRate);
        pargs.autoSaveFlag = false;        
        
        drzState = {};

        u = 1;
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(Trial,pfsArgs{:});
        pfTemp.save();
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        end
        pfTemp.save();        

    end
    pfd{tind,pfindex} = MTAApfs(Trial,'tag',pargs.tag);
    
% DISPLAY pfd BPITCHxBSPEED
    if display,  req20180123_vis_BPITCHxBSPEED;  end

end%for tind



% COMPUTE erpPCA_BPITCHxBSPEED eigenvectors ------------------------------------------------------------------
% CREATE SELCETED VECTOR SUBSPACE
if (nargout>1||display)&&numel(sessionList)==numel(Trials)
%rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,pfindex));
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
rmaps = cat(2, rmaps{:});
[smnd,sind] = sort(sum(isnan(rmaps),1));
selectedUnits = sind(smnd>1000&smnd<1450);
zrmaps = rmaps(:,selectedUnits);
zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/2;        % SELECT subspace {1/3 non-zero samples} 
zrmaps = zrmaps(validDimsInds,:);                   % REDUCE to selected subspace
zrmaps = 10.*(zrmaps>0)+randn(size(zrmaps)).*(zrmaps>0);

% DECOMPOSE 
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                 % COMPUTE covariance-based PCA with Varimax rotation
V = LR;

UID{pfindex} = selectedUnitsMat;
LRG{pfindex} = LR;
VTG{pfindex} = VT;
VSrG{pfindex} = FSr;


if display, req20180123_vis_BPITCHxBSPEED_erpPCA;  end
end
% ERPPCA BPITCHxBSPEED END -------------------------------------------------------------------------





% ERPPCA BPITCHxHSPEED START -------------------------------------------------------------------------

analDir = ['BPITCHxHSPEED-v',version];
create_directory(fullfile(figDir,analDir));

pfindex = 4;

for tind = 1:numTrials,
    Trial = Trials{tind}; %15,16,17,18    

    if display || overwrite,
        pft = pfs_2d_theta(Trial);
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('RectFilter');
        pch = fet_HB_pitchB(Trial);
        vxy = xyz.vel(['hcom'],[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);
        fet = pch.copy('empty');
        fet.label = 'fet_VP';       
        fet.data = [pch(:,2),vxy(:)];
        drz = compute_drz(Trial,units{tind},pft);%,pfstats);
        tper =[Trial.stc{'theta-groom-sit'}];
        tper.resample(xyz);
    end

% COMPUTE pfd BPITCHxHSPEED
    pargs.tag            = ['DRZxBPITCHxHSPEED_v',version];
    if overwrite,  
        pargs = get_default_args('MjgER2016','MTAApfs','struct');        
        pargs.tag              = ['DRZxBPITCHxHSPEED_v',version];
        pargs.units            = units{tind};
        pargs.numIter          = 1001;
        pargs.halfsample       = true;
        pargs.overwrite        = true;
        pargs.boundaryLimits   = [-2,2;-2,2];
        pargs.binDims          = [0.1,0.1];
        pargs.SmoothingWeights = [1.5,1.5];
        pargs.xyzp = MTADxyz('data',fet.data,'sampleRate',xyz.sampleRate);
        pargs.autoSaveFlag = false;
        
        drzState = {};

        u = 1;
        dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                         xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
        drzState{u} = dper&tper;
        pargs.units  = units{tind}(u);
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(Trial,pfsArgs{:});
        pfTemp.save();
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        end
        pfTemp.save();        

    end
    pfd{tind,pfindex} = MTAApfs(Trial,'tag',pargs.tag);
    
% DISPLAY pfd BPITCHxHSPEED
    if display,  req20180123_vis_BPITCHxHSPEED;  end

end%for tind

if (nargout>1||display)&&numel(sessionList)==numel(Trials)
% COMPUTE erpPCA_BPITCHxHSPEED eigenvectors
% CREATE SELCETED VECTOR SUBSPACE
rmaps = cf(@(p,u) mean(p.data.rateMap(:,ismember(p.data.clu,u),:),3,'omitnan'), pfd(:,pfindex),units');
rmaps = cat(2, rmaps{:});
[smnd,sind] = sort(sum(isnan(rmaps),1));
selectedUnits = sind(smnd>1000&smnd<1400);
zrmaps = rmaps(:,selectedUnitsMat);
zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/2;        % SELECT subspace {1/3 non-zero samples} 
zrmaps = zrmaps(validDimsInds,:);                   % REDUCE to selected subspace
zrmaps = 10.*(zrmaps>0)+randn(size(zrmaps)).*(zrmaps>0);

% DECOMPOSE 
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                 % COMPUTE covariance-based PCA with Varimax rotation
V = LR;

UID{pfindex} = selectedUnits;
LRG{pfindex} = LR;
VTG{pfindex} = VT;
VSrG{pfindex} = FSr;


if display, req20180123_vis_BPITCHxHSPEED_erpPCA;  end
end
% ERPPCA BPITCHxHSPEED END -------------------------------------------------------------------------




%GUNITS = intersect(intersect(intersect(UID{1},UID{2}),UID{3}),UID{4});



%% pfd plot parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTS  

if display, req20180123_vis_parts;  end