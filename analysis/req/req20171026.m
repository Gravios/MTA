% req20171026 ----------------------------------------------------
%  Status: active
%  Type: Test
%  Final_Forms: NA
%  Description: Debug and Characterize the label_bhv_reduced_aux classifier
%  Bugs: NA

% Metrics: Segmentation variance of high and low state
%          Cogruence between original and aux states
%          magnitude of Union between 

%sessionListName = 'ncp';
%sessionListName = 'BHV_S4H5';
sessionListName = 'MjgER2016';


sessionList = get_session_list(sessionListName);
numTrials   = numel(sessionList);
sessionList = mat2cell(sessionList,1,ones([1,numTrials]));
Trials = cf(@(t)   MTATrial.validate(t)                    ,sessionList);
stc    = cf(@(t)   t.load('stc','msnn_ppsvd')              ,Trials);
stc    = cf(@(s,t) label_bhv_reduced(s,t), stc,Trials);


stc    = cf(@(t)   t.load('stc','msnn_ppsvd_raux')              ,Trials);

batch_compute_pfstats_bs(sessionListName);