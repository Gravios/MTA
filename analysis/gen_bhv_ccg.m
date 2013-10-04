function  Bccg = gen_bhv_ccg(Trial,varargin)
[state,dur_thresh] = DefaultArgs(varargin,{'rear',180});

sper = Trial.Bhv.getState(state).state;

sdur = diff(sper,1,2);
sper = sper(sdur>dur_thresh,:);

Bccg = MTAccg(Trial,state,['CCG around' state 'and offset'], ...
              {sper(:,1),sper(:,2)},{1,2},{[state ' onset'],[state ' offset']});