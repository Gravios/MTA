function linkSessionData(varargin)
% function linkSessionData(sessionName,nlx_path,xyz_path)
%
% linkSessionData takes two paths and links
%
% Work in Progress

Session = MTASession([]);

%% DefaultArgs


% Assumes the nlx data is in the data path under nlx
ses_path = fullfile(Session.path.data,sessionName);


%% Setup session nlx folder in users data root in case the
%  nlx data is stored elsewhere.
ses_nlx_path = fullfile(Session.path.data,'nlx',sessionName);
if ~exist(ses_nlx_path,'dir'), 
    mkdir(sessionName); 
    cd(ses_nlx_path); 
    system(['ln -s ' nlx_path '/* .']);
end



system(['ls ' nlx_path '/*'])

if ~exist(ses_path,'dir'), mkdir(ses_path); end

cd(ses_path)
system(['ln -s '  '/* .']);