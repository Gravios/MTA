function MTAConfiguration(root_dir,varargin)
% MTAConfiguration(root_dir,flag)
% Sets up the directories, and variable lists required for MTA
%
% Variables:
%   root_dir: string, Path to the data directory
%
%   flag:     string, Path options
%     values: absolute,  assumes the given path starts in the root dir
%             otherwise, assumes the given path starts in the home dir
%
% Example:
%   MTAConfiguration('/data/home/username/data_root_dir','absolute')
%
%   MTAConfiguration('C:\Users\username\data_root_dir','absolute')
%
%   MTAConfiguration('data_root_dir','absolute')
%
% Configuration Files:
%    MTAPaths.mat  
%       The paths where the data and config folders exist
%
%    MTAMarkers.mat
%       The marker list whith valid marker names (May remove in the future)
%
%    MTAMazes.mat 
%       Cell array containing all information over mazes (Will be replaced by a
%       database)
%
%    MTAMarkerConnections.mat
%       Defined connections between markers (Not required if vsk file is
%       present)
% 

global MTA_PROJECT_PATH % use this in the future to set root dir
    
[flag,project_name,host_server,data_server,overwrite] = DefaultArgs(varargin,{'','','','',true});



if ispc, 
    userdir= getenv('USERPROFILE'); 
else
    userdir= getenv('HOME');
end

% SET path to MTA toolbox directory
mtap = fileparts(mfilename('fullpath'));


% SET paths to the data directory
switch flag
    case 'absolute'
        data = root_dir;
    otherwise
        data = fullfile(userdir,root_dir);
end


if ~exist(fullfile(data,'config'),'dir')
    mkdir(fullfile(data,'config'));
end

cfg = fullfile(data,'config','MTA');
if ~exist(cfg,'dir'),
    mkdir(cfg);
end


% SETUP paths
analysisDir = fullfile(MTA_PROJECT_PATH,'analysis');
if ~exist(analysisDir,'dir'),
    mkdir(analysisDir);
end

matlabDir = fullfile(MTA_PROJECT_PATH,'matlab');
if ~exist(matlabDir,'dir'),
    mkdir(matlabDir);
end

classifiersDir = fullfile(MTA_PROJECT_PATH,'matlab','classifiers');
if ~exist(classifiersDir,'dir'),
    mkdir(classifiersDir);
end

analysisModelDir = fullfile(analysisDir,'models');
if ~exist(analysisModelDir,'dir'),
    mkdir(analysisModelDir);
end

arm = fullfile(cfg,'arm');
if ~exist(arm,'dir')||overwrite
    copyfile(fullfile(mtap,'arm'),arm);
end

web = fullfile(cfg,'web');
if ~exist(web,'dir')||overwrite
    copyfile(fullfile(mtap,'web'),web);
end

% final path layout
%
% ../root_dir/config/
% ../root_dir/config/MTA
% ../root_dir/config/MTA/arm
% ../root_dir/config/MTA/web
%
% ../MTA_PROJECT_PATH/analysis
% ../MTA_PROJECT_PATH/analysis/models
% ../MTA_PROJECT_PATH/matlab
% ../MTA_PROJECT_PATH/matlab/classifiers







%% List of the accepted marker names
MTAMarkers ={'hip_right','pelvis_root','hip_left','spine_lower','knee_left',...
             'knee_right','tail_base','tail_end','spine_middle','spine_upper',...
             'head_back','head_left','head_front','head_right','shoulder_left',...
             'elbow_left','shoulder_right','elbow_right','head_top'};

%% List of the accepted mazes  
MTAMazes = {{'cof','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'chr','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'cpm','circle',   [-1000,1000;-1000,1000;-100,500],[-1000,1000;-1000,1000;-100,500]},...
            {'cpf','circle',   [-1000,1000;-1000,1000;-100,500],[-1000,1000;-1000,1000;-100,500]},...
            {'hcf','rectangle',[-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'nor','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'odr','W',        [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'cot','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'ont','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'sof','square',   [-700,700;-700,700;0,1000],[-700,700;-700,700;0,1000]},...
            {'sov','square',   [-700,700;-700,700;0,1000],[-700,700;-700,700;0,1000]},...
            {'spb','square',   [-700,700;-700,700;-100,500],[-700,700;-700,700;-100,500]},...
            {'shr','rectangle',[-300,500;-300,300;0,300],[-300,300;-300,300;0,360]},...
            ... %{'rof','rectangle',[-1100,1100;-600,600;0,400],[-1100,1100;-600,600;0,400]},...
            {'rov','rectangle',[-800,800;-500,500;0,300],[-600,600;-350,350;0,300]},...
            {'tfr','rectangle',[-800,800;-500,500;-100,300],[-800,800;-500,500;-100,300]},...
            {'frr','rectangle',[-700,700;-600,600;-100,300],[-700,700;-600,600;-100,300]},...
            {'ofr','rectangle',[-700,700;-600,600;-100,300],[-700,700;-600,600;-100,300]},...
            {'orr','rectangle',[-700,700;-600,600;-100,300],[-700,700;-600,600;-100,300]},...
            {'rof','rectangle',[-700,700;-600,600;-100,300],[-700,700;-600,600;-100,300]},...
            {'tm','circle',    [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]}};

%% List of the accepted connections between markers
MTAMarkerConnections = ...
   {{'head_front','head_left' },...
    {'head_front','head_right'},...
    {'head_right','head_left' },...
    {'head_left' ,'head_back' },...
    {'head_right','head_back' },...
    {'head_back' ,'spine_upper'},...
    {'spine_upper','shoulder_right'},...
    {'spine_upper','shoulder_left' },...
    {'spine_upper','spine_middle'  },...
    {'spine_middle','pelvis_root'  },...
    {'pelvis_root', 'hip_right'  },...
    {'pelvis_root', 'hip_left'   },...
    {'pelvis_root', 'spine_lower'},...
    {'spine_lower', 'hip_right'  },...
    {'spine_lower', 'hip_left'},...
    {'hip_left','hip_right'},...
    {'hip_left','knee_left'},...
    {'hip_right','knee_right'},...
    {'head_top','head_front'},...
    {'head_top','head_left'},...
    {'head_top','head_right'},...
    {'head_top','head_back'}, ...
    {'tail_proximal','spine_lower'}};

%% Save configuations

if ~exist(fullfile(cfg, 'MTAConf.mat'),'file')||overwrite,
    save(fullfile(cfg, 'MTAConf.mat'  ),'project_name','host_server','data_server');end

if ~exist(fullfile(cfg, 'MTAPaths.mat'),'file')||overwrite,
    save(fullfile(cfg, 'MTAPaths.mat'),'cfg','data','arm','web');end

if ~exist(fullfile(cfg, 'MTAMarkers.mat'),'file')||overwrite,
    save(fullfile(cfg, 'MTAMarkers.mat'),'MTAMarkers');end

if ~exist(fullfile(cfg, 'MTAMazes.mat'),'file')||overwrite,
    save(fullfile(cfg, 'MTAMazes.mat'),'MTAMazes');end

if ~exist(fullfile(cfg, 'MTAMarkerConnections.mat'),'file')||overwrite,
    save(fullfile(cfg, 'MTAMarkerConnections.mat'),'MTAMarkerConnections');end



addpath(genpath(cfg));
% $$$ try
% $$$     savepath
% $$$ end


end