function  Bccg = gen_bhv_ccg(Trial,varargin)
[state,dur_thresh] = DefaultArgs(varargin,{'rear',1.5});

sper = Trial.stc.filter(Trial.lfp.sampleRate,{state,{'exclusion',{state},2},{'duration',dur_thresh}});

if ~iscell(sper), sper = mat2cell(sper,size(sper,1),[1,1]);end

Bccg = MTAccg(Trial,state,['CCG around' state 'and offset'], ...
              sper,{[state ' onset'],[state ' offset']},'overwrite',true,...
              'normalization','count');