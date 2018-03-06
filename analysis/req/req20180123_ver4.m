


% SET Analysis Parameters
version = '4';
overwrite = true;
display = true;
sessionListName = 'MjgER2016';
figDir = create_directory('/storage/gravio/figures/analysis/placefields_nonSpatialFeatures'); 

% LOAD session list
% LOAD Trials
% SELECT units for analysis
% LOAD theta placefields
sessionList = get_session_list(sessionListName);

Trials  = af(@(S)  MTATrial.validate(S),   sessionList);
          cf(@(T)  T.load('nq'),           Trials);
units   = cf(@(T)  select_placefields(T),  Trials);

pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);

numTrials = numel(Trials);


%% pfd_BPITCHxHPITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERPCA BPITCHxHPITCH -----------------------------------------------------------------------------
analDir = 'BPITCHxHPITCH';
pfd = {};
pfindex = 1;
for tind = 1:numTrials,
    Trial = Trials{tind}; 

% LOAD vars
    xyz = preproc_xyz(Trial,'trb');
    pch = fet_HB_pitchB(Trial);
    drz = compute_drz(Trial,units{tind},pft{tind});%,pfstats);
    tper =[Trial.stc{'theta-groom-sit'}];
    tper.resample(xyz);

% COMPUTE HPITCH x BPITCH | DRZ [-0.5,0.5]
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units   = units{tind};
    pargs.numIter = 1001;
    pargs.halfsample = true;
    pargs.tag            = 'DRZxHBPITCHxBPITCH_v4';
    pargs.boundaryLimits = [-2,2;-2,2];
    pargs.binDims        = [0.2,0.2];
    pargs.SmoothingWeights = [1.5,1.5];
    if overwrite,  
        pargs.overwrite = true;    
        pargs.xyzp = MTADxyz('data',pch.data,'sampleRate',xyz.sampleRate);
        drzState = {};
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            MTAApfs(Trial,pfsArgs{:});
        end
    end
    pargs.states    = 'tdrz';
    pargs.units     = units{tind};
    pargs.overwrite = false;
    pfsArgs = struct2varargin(pargs);
    pfd{tind,pfindex} = MTAApfs(Trial,pfsArgs{:});

% VISUALIZE HPITCH x BPITCH | DRZ [-0.5,0.5]
    if display,  req20180123_vis_HPITCHxBPITCH;  end    
end%for tind

% COMPUTE erpPCA_HPITCHxBPITCH eigenvectors
rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,1));
rmaps = cat(2, rmaps{:});
si    = cf(@(p) p.data.si(:,:,1), pfd(:,1));
si = cat(2, si{:});
[~,sind] = sort(si);
[smnd,sind] = sort(sum(isnan(rmaps),1));

zrmaps = rmaps(:,sind(smnd>600));
zdims = size(zrmaps);                                                % NOTE the dimensions of the original vector space  
zrmaps(isnan(zrmaps)) = 0;                                           % SET all nan valued elements to zeros              
validDimsInds = sum(zrmaps==0,2)<zdims(2)/3;                         % SELECT subspace {1/3 non-zero samples}            
zrmaps = zrmaps(validDimsInds,:);                                    % REDUCE to selected subspace                       
                                                                                                                         
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                                                                                      
V = LR;                                                              % COMPUTE covariance-based PCA with Varimax rotation

if display,  req20180123_vis_HPITCHxBPITCH_erpPCA;  end

% ERPPCA HPITCHxBPITCH END ---------------------------------------------------------------------------





%% erpPCA_HPITCHxBSPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERPPCA HPITCHxBSPEED START -------------------------------------------------------------------------


pfindex = 2;

for tind = 1:numTrials,
    Trial = Trials{tind}; %15,16,17,18

    if display || overwrite,
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('RectFilter');
        pft = pfs_2d_theta(Trial,'overwrite',false);
        pch = fet_HB_pitchB(Trial);
        vxy = xyz.vel(['spine_lower'],[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);
        fet = pch.copy('empty');
        fet.label = 'fet_VP';       
        fet.data = [pch(:,2),vxy(:)];
        drz = compute_drz(Trial,units{tind},pft{tind});%,pfstats);
        tper =[Trial.stc{'theta-groom-sit-rear'}];
        tper.resample(xyz);
    end

% COMPUTE pfd HBPITCHxBSPEED
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units   = units{tind};
    pargs.numIter = 101;
    pargs.halfsample = true;
    pargs.tag            = 'DRZxHBPITCHxBSPEED_v4';
    pargs.boundaryLimits = [-2,2;-2,2];
    pargs.binDims        = [0.2,0.2];
    pargs.SmoothingWeights = [1.5,1.5];
    if overwrite,  
        pargs.overwrite = true;    
        pargs.xyzp = MTADxyz('data',fet.data,'sampleRate',xyz.sampleRate);
        drzState = {};
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            MTAApfs(Trial,pfsArgs{:});
        end
    end
    pargs.states    = 'tdrz';
    pargs.units     = units{tind};
    pargs.overwrite = false;
    pfsArgs = struct2varargin(pargs);
    pfd{tind,pfindex} = MTAApfs(Trial,pfsArgs{:});

% DISPLAY pfd HBPITCHxBSPEED
    if display,  req20180123_vis_HPITCHxBSPEED;  end

end%for tind



% COMPUTE erpPCA_HPITCHxBSPEED eigenvectors ------------------------------------------------------------------
% CREATE SELCETED VECTOR SUBSPACE
rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,pfindex));
si    = cf(@(p) p.data.si(:,:,1),      pfd(:,pfindex));
rmaps = cat(2, rmaps{:});
si    = cat(2, si{:}   );
[~,sind] = sort(si);
[smnd,sind] = sort(sum(isnan(rmaps),1));
zrmaps = rmaps(:,sind(smnd>800));
zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/2;        % SELECT subspace {1/3 non-zero samples} 
zrmaps = zrmaps(validDimsInds,:);                   % REDUCE to selected subspace

% DECOMPOSE 
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                 % COMPUTE covariance-based PCA with Varimax rotation
V = LR;

if display, req20180123_vis_HPITCHxBSPEED_erpPCA;  end

% ERPPCA HPITCHxSPEED END -------------------------------------------------------------------------




%% pfd plot parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTS  

if display, req20180123_vis_parts;  end