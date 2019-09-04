function MTAstartup(varargin)
%function MTAstartup(varargin)
%
% variables:
%   
%   hostServer: string, your personal tag which denotes the working
%                        environment you are currently using.
%
%   dataServer: string, your personal tag which denote in which
%                        mounted, file system your target data is located.
%
%   addBasicPaths: logical, flag which adds a set of paths which
%                    are important to the functionality of MTA,
%                    please ensure they point to the correct locations.
%                    (These paths can be found at the end of the script)
%
%   configure: logical, update mazes, models, ect...
%

clearvars('-except','varargin','MTA_CURRENT_PROJECT','MTA_PROJECT_PATH')

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('projectName',                      '',                                         ...
                 'hostServer',                       'lmu',                                      ...
                 'dataServer',                       'lmu',                                      ...
                 'addBasicPaths',                    true,                                       ...
                 'configure',                        true                                        ...
);
[projectName,hostServer,dataServer,addBasicPaths,configure] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

global MTA_PATH;
global MTA_DATA_PATH;
global MTA_PROJECT_PATH;
global MTA_CURRENT_PROJECT;

MTA_PATH            = getenv('MTA_PATH');
MTA_DATA_PATH       = getenv('MTA_DATA_PATH');
MTA_PROJECT_PATH    = getenv('MTA_PROJECT_PATH');

if ~isempty(projectName),
    MTA_CURRENT_PROJECT = projectName;
else
    configure = false;
    addBasicPaths = false;
end


switch hostServer % The severer where matlab is running. 
                   % Note: the addpath statement must point to 
                   % the "local" version of MTA.
  
  case 'lmu'
% ADD MTA to matlab paths
    addpath(MTA_PATH);

    switch dataServer % Where the data is located
      case 'lmu'
% SET project path        
        projectPath = fullfile(MTA_PROJECT_PATH,projectName);        
        if ~exist(projectPath),
            create_directory(projectPath);
        end
        if configure, 
% CONFIGURE project directories
            MTAConfiguration(projectPath,'absolute',projectName,hostServer,dataServer); 
        end
% SET the project path to global var        
        MTA_PROJECT_PATH = projectPath;
    end

end


% most likely will need corrections when you set
% up MTA for the first time
if addBasicPaths 
    if ispc,
        userpath = getenv('HOMEPATH');
    else
        userpath = getenv('HOME');
    end

    switch dataServer
      case 'lmu'
        addpath(genpath('/storage/share/matlab/Third-Party_Toolboxes/HMM/hmmbox/'));
        addpath(genpath('/storage/share/matlab/Third-Party_Toolboxes/netlab/'));
        addpath(genpath('/storage/share/matlab/MTA/'));
        rmpath(genpath('/storage/share/matlab/MTA/.git'));
        rmpath(genpath('/storage/share/matlab/MTA/config'));
        cd('/storage/share/matlab/MTA/');
        addpath('home/antsiro/matlab/General','-END');
        addpath('home/antsiro/matlab/draft','-END');
      
      otherwise        
        addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/HMM/hmmbox/')));
        addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/netlab/')));
        addpath(genpath(fullfile(userpath,'../../homes/share/matlab/MTA/')));
        rmpath(genpath(fullfile(userpath,'../../homes/share/matlab/MTA/.git')));
        rmpath( '../sirota/homes/share/matlab/MTA/config');
        cd(fullfile(userpath,'../../homes/share/matlab/MTA/'));
        addpath(fullfile(userpath,'../../homes/antsiro/matlab/General'),'-END');
        addpath(fullfile(userpath,'../../homes/antsiro/matlab/draft'),'-END');
    end
end