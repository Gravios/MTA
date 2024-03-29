function MTASetup(root_dir,flag)

paths = MTAConfiguration(root_dir);
paths_check = struct2cell(paths);

% Make Directory Tree
for i = 1:length(paths_check)
    if ~exist(paths_check{i},'file'),
        mkdir(paths_check{i});
    end
end
addpath(paths.cfg)
savepath
system(['mv /tmp/MTA{Paths,Markers,Mazes,MarkerConnections}.mat ' paths.cfg ]) 

% Flag for linking data from other users
switch flag
  case 'link_to_gravio'
    session_list = '/data/homes/gravio/data/analysis/session_list';
    paths_to_gravio = load('/data/homes/gravio/data/config/MTA/MTAPaths.mat');
    gtdb('link_session_to_project',session_list,[],[],[],paths_to_gravio,paths);
  otherwise
    error(['Flag: ' flag ', not recognized']);
end

