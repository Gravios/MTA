diary off

%% Detect system and set root directories
if ispc,

else
    [~,hostname] = system('hostname');
    userdir= getenv('HOME');
    [~,username] = system('whoami');
    username = username(1:length(username)-1);

    filesystem(1) = {'sirota_lab'};
    host_list_regex_patterns(1) = {{'bach','liszt'}};

    filesystem(2) = {'cin_cluster'};
    host_list_regex_patterns(2) = {{'hn\d','cn\d'}};

    filesystem(3) = {'bw_cluster'};
    host_list_regex_patterns(2) = {{'bw?'}};

    hlrps = size(host_list_regex_patterns);
    host_list_regex_patterns = cellfun(@strjoin,host_list_regex_patterns,repmat({')|('},hlrps),'UniformOutput',false);
    host_list_regex_patterns = cellfun(@cat,repmat({2},hlrps),repmat({'('},hlrps),host_list_regex_patterns,repmat({')'},hlrps),'uniformoutput',false);
    hlid = ~cellfun(@isempty,regexp(hostname,host_list_regex_patterns));

    switch filesystem{find(hlid)}
      case 'sirota_lab'
        path.home = userdir;
        path.bach_data = fullfile(userdir,'data');
        path.bach_share_matlab = fullfile(userdir,'matlab')
        path.bach_share_matlab = fullfile('/data/homes/share/matlab/')
        path.cin_cluster_data = fullfile('/mnt/gpfs',username,'data');
      case 'cin_cluster'
        path.home = userdir;
        path.bach_data        = fullfile('/gpfs01/sirota/bach/homes/',username,'data')
        path.bach_share_matlab = fullfile('/gpfs01/sirota/bach/homes/share/matlab/')
        path.cin_cluster_data = fullfile('/gpfs01/sirota',username,'data');
        
      case 'bw_cluster'
        %path.home = userdir;
        .path.project = fullfile(userdir,'data');
        %path.cin_cluster_data = fullfile('/mnt/gpfs',username,'data');

    end

end

% Need one of these for the hmm stuf ... can't remember which one
% at the moment
%addpath(genpath('/gpfs01/sirota/bach/homes/gravio/root/data/homes/share/matlab/Third-Party_Toolboxes/circStat2010b'));
%addpath(genpath('/gpfs01/sirota/bach/homes/gravio/root/data/homes/share/matlab/Third-Party_Toolboxes/HMM/hmmbox'));
%addpath(genpath('/gpfs01/sirota/bach/homes/gravio/root/data/homes/share/matlab/Third-Party_Toolboxes/netlab'));

%% Labbox
addpath(genpath(fullfile(path.bach_share_matlab, 'labbox')));
 rmpath(genpath(fullfile(path.bach_share_matlab, 'labbox/.git'));

%% Personal utilities
addpath(fullfile(path.bach_matlab,'utilities'));

%% MTA
addpath /gpfs01/sirota/gravio/data/config/MTA
addpath(genpath('/gpfs01/sirota/bach/homes/gravio/MTA'),'-begin');
rmpath(genpath('/gpfs01/sirota/bach/homes/gravio/MTA/.git'));
addpath(genpath('/gpfs01/sirota/bach/homes/gravio/data/config/MTA'),'-begin');

%% Antons stuff
addpath /gpfs01/sirota/bach/homes/antsiro/matlab/General -END
addpath /gpfs01/sirota/bach/homes/antsiro/matlab/draft -END
