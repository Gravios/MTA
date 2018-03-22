function  [Bccg,sper] = gen_bhv_ccg(Trial,varargin)
[state,durThresh,units,numPartitions,pft] = DefaultArgs(varargin,{'rear',1.5,[],4,[]});


sper = Trial.stc.filter(Trial.lfp.sampleRate,{state,{'exclusion',{state},2},{'duration',durThresh}});

if ~iscell(sper), sper = mat2cell(sper,size(sper,1),[1,1]);end

if numPartitions>1,
    if isempty(pft),  pft = pfs_2d_theta(Trial,units);  end
    [~,mrp] = pft.maxRate();
else
    mrp = [];
end

try,
    Bccg = MTAccg(Trial,                                              ... Trial              
                  state,                                              ... name              
                  ['CCG around' state 'and offset'],                  ... Description              
                  sper,                                               ... ResTrain
                  {[state ' onset'],[state ' offset']},               ... CluTags
                  units,                                              ... units              
                  true,                                               ... overwrite
                  [],                                                 ... rand_tag
                  'abs_dist',                                         ... method
                  numPartitions,                                      ... partitions
                  mrp,                                                ... partition_feature
                  [],                                                 ... surrogate_sample_size
                  1,                                                  ... numIterations
                  200,                                                ... binSize
                  20,                                                 ... halfBins
                  'hz'                                                ... normalization
    );
catch err,
    disp(err);
    af(@(err) disp(err),err.stack);
    Bccg = {};
end
