function MTAstartup(varargin)
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','bach',true});
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','cin',true});
%[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'mypc','mycy',true});


[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','bach',true});

switch host_server
    
    case 'cin'
        switch data_server
          case 'cin'
              addpath('/gpfs01/sirota/homes/share/matlab/MTA/');
              MTAConfiguration('/gpfs01/sirota/home/gravio/data','absolute');
          case 'bach'
              addpath('/gpfs01/sirota/homes/share/matlab/MTA/');
              %MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
              MTAConfiguration('/gpfs01/sirota/data/bachdata/data/gravio','absolute');
        end
        
        
    case 'bach'
        switch data_server
          case 'cin'
            addpath('/data/homes/share/matlab/MTA/');
            MTAConfiguration('/mnt/gpfs/home/gravio/data','absolute');  
          case 'bach'
            MTAConfiguration('/gpfs01/sirota/data/bachdata/data/gravio','absolute');
            %MTAConfiguration('/data/homes/gravio/data','absolute');
        end

    case 'mypc'
        switch data_server
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

if add_basic_paths    
    if ispc,
        userpath = getenv('HOMEPATH');
    else
        userpath = getenv('HOME');
    end
    addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/HMM/hmmbox/')));
    addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/netlab/')));
    addpath(genpath(fullfile(userpath,'../../homes/share/matlab/MTA/')));
    rmpath(genpath(fullfile(userpath,'../../homes/share/matlab/MTA/.git')));
    addpath(fullfile(userpath,'../../homes/antsiro/matlab/General'),'-END');
    addpath(fullfile(userpath,'../../homes/antsiro/matlab/draft'),'-END');
    cd(fullfile(userpath,'../../homes/share/matlab/MTA/'));
end

% $$$ if add_basic_paths    
% $$$     if ispc,
% $$$         userpath = getenv('HOMEPATH');
% $$$     else
% $$$         userpath = getenv('HOME');
% $$$     end
% $$$     cd(fullfile(userpath,'../share/matlab/MTA/'));
% $$$     addpath(genpath(fullfile(userpath,'../share/matlab/Third-Party_Toolboxes/HMM/hmmbox/')),'-END')
% $$$     addpath(genpath(fillfile(userpath,'../share/matlab/Third-Party_Toolboxes/netlab/')),'-END')
% $$$ endC:\cygwin64\home