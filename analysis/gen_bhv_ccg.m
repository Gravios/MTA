function  Bccg = gen_bhv_ccg(Trial,varargin)
[state,dur_thresh] = DefaultArgs(varargin,{'rear',1.5});

sper = Trial.stc{state,Trial.lfp.sampleRate}.data;

sdur = diff(sper,1,2);
sper = sper(sdur>(dur_thresh.*Trial.lfp.sampleRate),:);

Bccg = MTAccg(Trial,state,['CCG around' state 'and offset'], ...
              {sper(:,1),sper(:,2)},{[state ' onset'],[state ' offset']});