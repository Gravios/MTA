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

arm = fullfile(cfg,'arm');
if ~exist(arm,'dir')||overwrite
    copyfile(fullfile(mtap,'arm'),arm);
end

web = fullfile(cfg,'web');
if ~exist(web,'dir')||overwrite
    copyfile(fullfile(mtap,'web'),web);
end







%% List of the accepted marker names
MTAMarkers ={'hip_right','pelvis_root','hip_left','spine_lower','knee_left',...
             'knee_right','tail_base','tail_end','spine_middle','spine_upper',...
             'head_back','head_left','head_front','head_right','shoulder_left',...
             'elbow_left','shoulder_right','elbow_right','head_top'};

%% List of the accepted mazes
MTAMazes = {{'cof','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'chr','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'nor','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'odr','W',        [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'cot','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'ont','circle',   [-500,500;-500,500;0,300],[-500,500;-500,500;0,360]},...
            {'sof','square',   [-700,700;-700,700;0,1000],[-700,700;-700,700;0,1000]},...
            {'sov','square',   [-700,700;-700,700;0,1000],[-700,700;-700,700;0,1000]},...
            {'shr','rectangle',[-300,500;-300,300;0,300],[-300,300;-300,300;0,360]},...
            {'rof','rectangle',[-1100,1100;-600,600;0,400],[-1100,1100;-600,600;0,400]},...
            {'rov','rectangle',[-800,800;-500,500;0,300],[-600,600;-350,350;0,300]},...
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