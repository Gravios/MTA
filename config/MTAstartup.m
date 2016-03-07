function MTAstartup(varargin)
%function MTAstartup(varargin)
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','bach',true});
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'mypc','mycy',true});
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'mypc','mycy',true});
%
% variables:
%   
%   host_server: string, your personal tag which denotes the working
%                        environment you are currently using.
%
%   data_server: string, your personal tag which denote in which
%                        mounted, file system your target data is located.
%
%   add_basic_paths: logical, flag which adds a set of paths which
%                    are important to the functionality of MTA,
%                    please ensure they point to the correct locations.
%                    (These paths can be found at the end of the script)

[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','cin',true});

switch host_server % The severer where matlab is running. 
                   % Note: the addpath statement must point to 
                   % the "local" version of MTA.
    
  case 'lmu'
    addpath('/storage/share/matlab/MTA/');

    switch data_server % Where the data is located
      case 'lmu'
        projPath = fullfile('/storage/gravio/data/project/',project_name);
        if ~exist(projPath),
            mkdir(projPath);
        end
        
        if ~isempty(project_name),            
            MTAConfiguration(projPath,'absolute');
        else
            MTAConfiguration(fullfile(projPath,'general'),'absolute');
        end
    end
  
  case 'cin'
    addpath('/gpfs01/sirota/homes/share/matlab/MTA/');

    switch data_server % Where the data is located
      case 'cin'
        MTAConfiguration('/gpfs01/sirota/home/gravio/data','absolute');
      case 'bach'
        %MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
        MTAConfiguration('/gpfs01/sirota/data/bachdata/data/gravio','absolute');
    end

  case 'bach'
    addpath('/data/homes/share/matlab/MTA/');

    switch data_server % Where the data is located
      case 'cin'
        MTAConfiguration('/mnt/gpfs/home/gravio/data','absolute');  
      case 'bach'
        MTAConfiguration('/gpfs01/sirota/data/bachdata/data/gravio','absolute');
        %MTAConfiguration('/data/homes/gravio/data','absolute');
    end

  case 'mypc'

    switch data_server % Where the data is located
      case 'mypc'
        MTAConfiguration('C:\Users\justi_000\data','absolute');        
      case 'myhd'
        MTAConfiguration('F:\data','absolute');        
      case 'mysd'
        MTAConfiguration('E:\data','absolute');      
      case 'mycy'
        MTAConfiguration('C:\cygwin64\home\justi_000\data','absolute');
    end
    return

end

% most likely will need corrections when you set
% up MTA for the first time
if add_basic_paths 
                   
    if ispc,
        userpath = getenv('HOMEPATH');
    else
        userpath = getenv('HOME');
    end

    switch data_server
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