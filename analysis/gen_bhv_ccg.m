function  [Bccg,sper] = gen_bhv_ccg(Trial,varargin)
[state,durThresh] = DefaultArgs(varargin,{'rear',1.5});


sper = Trial.stc.filter(Trial.lfp.sampleRate,{state,{'exclusion',{state},2},{'duration',durThresh}});

if ~iscell(sper), sper = mat2cell(sper,size(sper,1),[1,1]);end

try,
Bccg = MTAccg(Trial,state,['CCG around' state 'and offset'], ...
              sper,{[state ' onset'],[state ' offset']},'overwrite',true,...
              'normalization','hz');
catch err,
    disp(err);
    af(@(err) disp(err),err.stack);
    Bccg = {};
end
