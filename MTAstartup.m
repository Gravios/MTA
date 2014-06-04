function MTAstartup(varargin)
[host_server,data_server,add_basic_paths] = DefaultArgs(varargin,{'cin','bach',true});

switch host_server

  case 'cin'
    switch data_server
      case 'cin'
        MTAConfiguration('/gpfs01/sirota/home/gravio/data','absolute');
      case 'bach'
        MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
    end


  case 'bach'
    switch data_server
      case 'cin'
        MTAConfiguration('/mnt/gpfs/home/gravio/data','absolute');
      case 'bach'
        MTAConfiguration('/data/homes/gravio/data','absolute');
    end


end

if add_basic_paths
    cd('/gpfs01/sirota/bach/homes/gravio/MTA/');
    addpath(genpath('/gpfs01/sirota/bach/homes/share/matlab/Third-Party_Toolboxes/HMM/hmmbox/'))
    addpath(genpath('/gpfs01/sirota/bach/homes/share/matlab/Third-Party_Toolboxes/netlab/'))
end