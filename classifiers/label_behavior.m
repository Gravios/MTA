function label_behavior(sessionList)
%function label_behavior(sessionList)
% sessionList: (string) list of sessions found in get_session_list.m

% LABEL sessions with multi-session patternnet classifier 
label_bhv_msnn('msnn',get_session_list(sessionList,[],[],struct('stcMode','default')));

% OPTIMIZE state labels based on heuristics
optimize_stc_transition('msnn',sessionList);

% LABEL Shakes based on shaking of the body 
Trials = af(@(t) MTATrial.validate(t), get_session_list(sessionList));
cf(@(t)  label_bhv_shake('msnn_ppsvd',t),  Trials);

% $$$ for t = 1:numel(Trials),
% $$$     label_bhv_shake('msnn_ppsvd',Trials{t});
% $$$ end

% CREATE composite state of walk and turn
cf(@(t)  reduce_stc_to_loc(t.load('stc','msnn_ppsvd')),  Trials);

% LABEL sniffing periods


% LABEL homebase behavior
cf(@(t)  label_bhv_homebase([],t),  Trials);
