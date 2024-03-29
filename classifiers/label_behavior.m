function label_behavior(sessionList,varargin)
%function label_behavior(sessionList)
% sessionList: (string) list of sessions found in get_session_list.m

defargs = struct('featureSet','fet_mis','mode','fullbody');
[featureSet,mode] = DefaultArgs(varargin,defargs,'--struct');


switch mode

% REQUIRES marker set: spine_lower, pelvis_root, spine_middel, spine_upper, head_back, head_left,
%                      head_front, head_right
  case 'fullbody'
    Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
    cf(@(t)  label_behavior_with_heuristics(t), Trials);
% LABEL sessions with multi-session patternnet classifier 
    label_bhv_msnn('msnn',...
                   get_session_list(sessionList,[],[],struct('stcMode','default')),...
                   [],...
                   featureSet,...
                   'nIter',10,...
                   'randomizationMethod','WSBNT');
% OPTIMIZE state labels based on heuristics
    optimize_stc_transition('msnn',Trials);
    %optimize_stc_transition('msnn',Trials);    
% LABEL Shakes based on shaking of the body 
    cf(@(t)  label_bhv_shake('msnn_ppsvd',t),  Trials);
% CREATE composite state of walk and turn
    stc = cf(@(t)  reduce_stc_to_loc(t.load('stc','msnn_ppsvd')),  Trials);
% LABEL sniffing periods
    cf(@(s,t)  label_bhv_reduced(s,t),  stc,Trials);
% LABEL homebase behavior
    %cf(@(t)  label_bhv_homebase([],t),  Trials);

  case 'head'
    % todo implement req20190615
    
end
