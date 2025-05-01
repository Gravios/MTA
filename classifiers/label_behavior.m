function label_behavior(Trial,varargin)
%function label_behavior(sessionList)
% sessionList: (string) list of sessions found in get_session_list.m

defargs = struct('featureSet','fet_mis','mode','fullbody');
[featureSet,mode] = DefaultArgs(varargin,defargs,'--struct');


switch mode


    
% REQUIRES marker set: spine_lower, pelvis_root, spine_middel, spine_upper, head_back, head_left,
%                      head_front, head_right
  case 'fullbody'
    label_behavior_with_heuristics(Trial)
% LABEL sessions with multi-session patternnet classifier 
% $$$     label_bhv_msnn('msnn',...
% $$$                    get_session_list_v2(sessionList,[],[],struct('stcMode','default')),...
% $$$                    [],...
% $$$                    featureSet,...
% $$$                    'nIter',10,...
% $$$                    'randomizationMethod','WSBNT');
    label_bhv_msnn('msnn',                    ...
                   {Trial},                   ...
                   [],                        ...
                   featureSet,                ...
                   'nIter', 10,               ...
                   'randomizationMethod','WSBNT');

% OPTIMIZE state labels based on heuristics
    optimize_stc_transition('msnn',{Trial});
    %optimize_stc_transition('msnn',Trials);    
% LABEL Shakes based on shaking of the body 
    label_bhv_shake('msnn_ppsvd',Trial);
% CREATE composite state of walk and turn
% $$$ cf(@(t)  reduce_stc_to_loc('msnn_ppsvd',t),  Trials);
% LABEL sniffing periods
    label_bhv_reduced('msnn_ppsvd',Trial);
% LABEL homebase behavior
    %cf(@(t)  label_bhv_homebase([],t),  Trials);

  case 'head'
    % todo implement req20190615
    
end
