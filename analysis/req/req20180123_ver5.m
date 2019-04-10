function [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(varargin)
%function  req20180123_ver5(varargin)
% 
% T := number of Trials
% S := number of 
% Outputs:
%    pfd,              CellArray of MTAApfs-Object :[Tx
%    tags,             CellArray of String-Char    :
%    eigVec,           CellArray of Matrix-Numeric :
%    eigVar,           CellArray of Matrix-Numeric :
%    eigScore,         CellArray of Matrix-Numeric :
%    validDims,        CellArray of Matrix-Numeric :
%    unitSubsets,      CellArray of Matrix-Numeric :
%    unitIntersection, Matrix-Numeric :
%    zrmMean, 
%    zrmStd
%
% The structure of this analysis pipeline is sub-optimal 

%varargin = {};
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Trials',                        {{}},                                          ...
                 'sessionListName',               'MjgER2016',                                   ...
                 'version',                       '13',                                          ...
                 'overwritePfdFlag',              false,                                         ...
                 'overwriteErpPCAFlag',           false,                                         ...
                 'overwriteBhvContoursFlag',      false,                                         ...
                 'display',                       false,                                         ...
                 'figDir',  '/storage/gravio/figures/analysis/placefields_nonSpatialFeatures',   ...
                 'loadPFDFlag',                   true                                           ...
);
[Trials,sessionListName,version,overwritePfdFlag,overwriteErpPCAFlag,overwriteBhvContoursFlag,   ...
        display,figDir,loadPFDFlag] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

create_directory(figDir); 

% LOAD session list
% LOAD Trials
% FILTER units for analysis
sessionList = get_session_list(sessionListName);
if isa(Trials,'MTATrial'),
    Trials = {Trials};
    units = cf(@(T)  select_placefields(T),  Trials);    
elseif isempty(Trials)||numel(Trials)==numel(sessionList),% don't mess up here
    Trials  = af(@(S)  MTATrial.validate(S),   sessionList);
    units = cf(@(T)  select_placefields(T),  Trials);
    units = req20180123_remove_bad_units(units);
end


numTrials = numel(Trials);
numComp = 10;
numUnits = sum(cellfun(@numel,units));


% SET analysis paremeters
tags={}; fetSets={}; fetInds={}; pfdParam={}; states={}; ranges={};
%pfdParam = { boudaries   binDims   sweights };


%1. HPITCHxBPITCH version 14 running
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2,0.8;-0.8,2],[0.2,0.2],[2,2]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,300];



%1. HPITCHxBPITCH version 13 running
tags{end+1}        = ['HBPITCHxBPITCH'];
fetSets{end+1}     = 'fet_HB_pitchB';
fetInds{end+1}     = [1,2];
pfdParam{end+1}    = {[-1.8,1;-0.8,2],[0.1,0.1],[2,2]};
states{end+1}      = 'theta-groom-sit';
ranges{end+1}      = [0,300];


%1. HPITCHxBPITCH version 12 running
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2,1;-0.5,2],[0.5,0.5],[]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,20];


%1. HPITCHxBPITCH version 12 running
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2,1;-0.5,2],[0.5,0.5],[]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,20];


%1. HPITCHxBPITCH version 11 -- failed: need bigger bins
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2,1;-0.5,2],[0.25,0.25],[]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,80];


%1. HPITCHxBPITCH version 10
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2.25,1;-0.5,2],[0.25,0.25],[0.8,0.8]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,60];


%1. HPITCHxBPITCH version 9
% $$$ tags{end+1}        = ['HBPITCHxBPITCH'];
% $$$ fetSets{end+1}     = 'fet_HB_pitchB';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2.25,1;-0.5,2],[0.2,0.2],[1.1,1.1]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [0,110];%[1100,1550];%1450];

% $$$ %2. HPITCHxBSPEED
% $$$ tags{end+1}        = 'HBPITCHxBSPEED';
% $$$ fetSets{end+1}     = 'fet_HB_HPS';
% $$$ fetInds{end+1}     = [1,3];
% $$$ pfdParam{end+1}    = {[-2.25,1;-3,2],[0.1,0.1],[0.8,0.8]};
% $$$ states{end+1}      = 'theta-groom-sit-rear';
% $$$ ranges{end+1}      = [650,1350];
% $$$ 
% $$$ %2. HPITCHxHSPEED
% $$$ tags{end+1}        = 'HBPITCHxHSPEED';
% $$$ fetSets{end+1}     = 'fet_HB_HPS';
% $$$ fetInds{end+1}     = [1,4];
% $$$ pfdParam{end+1}    = {[-2.25,1;-3,2],[0.1,0.1],[0.8,0.8]};
% $$$ states{end+1}      = 'theta-groom-sit-rear';
% $$$ ranges{end+1}      = [650,1350];

%2. 'BSPEEDxHSPEED'
% $$$ tags{end+1}        = 'BSPEEDxHSPEED';
% $$$ fetSets{end+1}     = 'fet_HB_HPS';
% $$$ fetInds{end+1}     = [3,4];
% $$$ pfdParam{end+1}    = {[-2,2;-2,2],[0.1,0.1],[1.5,1.5]};
% $$$ states{end+1}      = 'theta-groom-sit-rear';
% $$$ ranges{end+1}      = [950,1550];

% $$$ 
% $$$ %3. BPITCHxBSPEED
% $$$ tags{end+1}        = 'BPITCHxBSPEED';
% $$$ fetSets{end+1}     = 'fet_HB_HPS';
% $$$ fetInds{end+1}     = [2,3];
% $$$ pfdParam{end+1}    = {[-2,2;-2,2],[0.1,0.1],[1.5,1.5]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [1000,1450];
% $$$ 
% $$$ %4. BPITCHxHSPEED
% $$$ tags{end+1}        = 'BPITCHxHSPEED';
% $$$ fetSets{end+1}     = 'fet_HB_HPS';
% $$$ fetInds{end+1}     = [2,4];
% $$$ pfdParam{end+1}    = {[-2,2;-2,2],[0.1,0.1],[1.5,1.5]};
% $$$ states{end+1}      = 'theta-groom-sit';
% $$$ ranges{end+1}      = [1000,1400];
% $$$ 
% $$$ %5. HPITCHxRHM   
% $$$ tags{end+1}        = 'HPITCHxRHM';
% $$$ fetSets{end+1}     = 'fet_HB_HPR';
% $$$ fetInds{end+1}     = [1,2];
% $$$ pfdParam{end+1}    = {[-2,2;-9,-3],[0.1,0.2],[1.5,1.5]};
% $$$ states{end+1}      = 'theta-groom-sit-rear';
% $$$ ranges{end+1}      = [700,1100];
% $$$                      


% COMPUTE and DISPLAY drz restricted rate maps
pfd = {};

if loadPFDFlag,
for pfindex = 1:numel(tags),
    create_directory(fullfile(figDir,[tags{pfindex},'_v',version]));  
    analDir = [tags{pfindex},'_v',version];
    for tind = 1:numTrials,
        Trial = Trials{tind};
% LOAD features
        if (display || overwritePfdFlag),        
            pft = pfs_2d_theta(Trial);
            xyz = preproc_xyz(Trial,'trb');
            drz = compute_drz(Trial,units{tind},pft);
            ddz = compute_ddz(Trial,units{tind},pft);
            tper =[Trial.stc{states{pfindex}}];        
            tper.resample(xyz);        
            fet = feval(fetSets{pfindex},Trial);
            fet.data = fet(:,fetInds{pfindex});
        end    
% COMPUTE pfd 
        if overwritePfdFlag,  
            [drzState] = req20180123_pfd_compute(                    ... COMUTE DRZ selected rate Maps
                Trial,                                               ... Trial
                fet,                                                 ... xyzp
                [tags{pfindex},'_v',version],                        ... analysis tag
                units{tind},                                         ... units
                drz,                                                 ... directed rate zones
                ddz,                                                 ... directed distance zones
                tper,                                                ... theta periods
                pfdParam{pfindex}{:}                                 ... 
            );
            [drzStateShuff] = req20180123_pfd_shuffled(              ... COMUTE DRZ selected rate Maps
                Trial,                                               ... Trial
                fet,                                                 ... xyzp
                [tags{pfindex},'_v',version],                        ... analysis tag
                units{tind},                                         ... units
                drz,                                                 ... directed rate zones
                ddz,                                                 ... directed distance zones                
                tper,                                                ... theta periods
                pfdParam{pfindex}{:}                                 ... 
            );
        end
        
        
% LOAD pfds        
        pfd{tind,pfindex} = MTAApfs(Trial,'tag',[tags{pfindex},'_v',version]);    
% PLOT pfds    
        if display, eval(['req20180123_vis_',tags{pfindex}]); end  

    end%for tind
end%for pfindex
end


%% bhv space erpPCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout>2||display) && numel(sessionList)==numel(Trials),
% COMPUTE erpPCA on units firing rates conditioned on behavioral space
    for pfindex = 1:numel(tags),
        [eigVec{pfindex},eigScore{pfindex},eigVar{pfindex},...
         unitSubsets{pfindex},validDims{pfindex},zrmMean{pfindex},zrmStd{pfindex}] = ...
            req20180123_pfd_erpPCA(pfd,units,[tags{pfindex},'_v',version],ranges{pfindex},...
                                   pfindex,numComp,overwriteErpPCAFlag);
% PLOT eigenvectors with behavioral state contours
        if display||overwriteErpPCAFlag, eval(['req20180123_vis_',tags{pfindex},'_erpPCA']); end
    end
    
% GET intersections of units common to each erpPCA analysis
    unitIntersection = unitSubsets{1};
    for pfindex = 2:numel(tags),
        unitIntersection = intersect(unitIntersection,unitSubsets{pfindex});
    end
end



% PLOT pfd parts ----------------------------------------------------------------------------------------
if display, 
    req20180123_vis_parts;  
end