function MTAstartup(varargin)
[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','bach',true});

switch host_server
    
    case 'cin'
        switch data_server
            case 'cin'
                MTAConfiguration('/gpfs01/sirota/home/gravio/data','absolute');
            case 'bach'
              %MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
                MTAConfiguration('/gpfs01/sirota/data/bachdata/data/gravio','absolute');
        end
        
        
    case 'bach'
        switch data_server
            case 'cin'
                MTAConfiguration('/mnt/gpfs/home/gravio/data','absolute');
            case 'bach'
                MTAConfiguration('/data/homes/gravio/data','absolute');
        end

    case 'mypc'
        switch data_server
            case 'mypc'
                MTAConfiguration('C:\Users\justi_000\data','absolute');        
            case 'myhd'
                MTAConfiguration('F:\data','absolute');        
            case 'mysd'
                MTAConfiguration('E:\data','absolute');        
        end
        return
        
        
end

if add_basic_paths    
    if ispc,
        userpath = getenv('HOMEPATH');
    else
        userpath = getenv('HOME');
    end
    cd(fullfile(userpath,'../../homes/share/matlab/MTA/'));
    addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/HMM/hmmbox/')),'-END')
    addpath(genpath(fullfile(userpath,'../../homes/share/matlab/Third-Party_Toolboxes/netlab/')),'-END')
end
% $$$ 
% $$$ if add_basic_paths    
% $$$     if ispc,
% $$$         userpath = getenv('HOMEPATH');
% $$$     else
% $$$         userpath = getenv('HOME');
% $$$     end
% $$$     cd(fullfile(userpath,'../share/matlab/MTA/'));
% $$$     addpath(genpath(fullfile(userpath,'../share/matlab/Third-Party_Toolboxes/HMM/hmmbox/')),'-END')
% $$$     addpath(genpath(fillfile(userpath,'../share/matlab/Third-Party_Toolboxes/netlab/')),'-END')
% $$$ end