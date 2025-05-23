
% Combined walk and turn into loc
states = {'loc','rear','pause','groom','sit'};

%% convert hand labeled to a reduced state set
tlist = get_session_list('hand_labeled_Ed');
tlist = get_session_list('hand_labeled_jg');
for t = 1:numel(tlist),
    Trial = MTATrial.validate(tlist(t));
    stc = Trial.load('stc',tlist(t).stcMode);
    strprts = regexp(stc.mode,'[_]','split');
    stc.updateMode([strprts{1}(1),strprts{2}(1),'_',...
                    cell2mat(regexp(strprts{3},'\d+','match')),'_',...
                    strprts{4},'_r']);

    mper = stc{'w'};
    mper.data = stc{'w+n'}.data;
    mper.name = 'locomotion';
    mper.label = 'loc';
    mper.key   = 'x';
    mper.updateFilename(Trial);

    stc.states(stc.gsi({'w','n'})) = [];
    stc.states = cat(2,{mper},stc.states(:)');

    stc.save(1);

end


defargs = struct('fetSet',   'fet_mis',                                      ...
                 'mode',     'train',                                        ...
                 'tag_preprocessing',   '+seh+',                             ...
                 'tag_postprocessing',  '',                                  ...
                 ... 'rlist',    'training_hand_labeled_reduced',                ...
                 'rlist',        'hand_labeled_Ed_reduced',                ...                 
                 'slist',     {{'hand_labeled_jg_reduced';'hand_labeled_Ed_reduced'}},       ...
                 'sampleRate',12,                                            ...
                 'nNeurons',  100,                                           ...
                 'nIter',     100,                                           ...
                 'states',    {states},        ...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',      true,                                          ...
                 'mref',      true,                                          ...
                 'prctTrain', []                                             ...
);
defargs = struct2varargin(defargs);
req20160128(defargs{:});

defargs{4} = 'compute';
req20160128(defargs{:});

defargs{4} = 'display';
req20160128(defargs{:});



t = 3;
Trial = MTATrial.validate(tlist(t));
stc = Trial.load('stc',tlist(t).stcMode);
mper = stc{'w'};
mper.data = stc{'w+n'}.data;
mper.name = 'locomotion';
mper.label = 'loc';
mper.key   = 'x';
mper.updateFilename(Trial);
stc.states = cat(2,{mper},stc.states(:)');

ds = load('/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat');

stcn = ds.stc{t}.copy;
mper = stcn{'w'};
mper.data = stcn{'w+n'}.data;
mper.name = 'locomotion';
mper.label = 'loc';
mper.key   = 'x';
mper.updateFilename(Trial);
stcn.states = cat(2,{mper},stcn.states(:)');


labelStats = cmp_stcs(Trial,stc,stcn,states,true);




%% Optimization of Stc
defargs = struct('rlist',    'training_hand_labeled_reduced',                ...
                 'slist',     {{'hand_labeled_jg_reduced';'hand_labeled_Ed_reduced'}},       ...
                 'fetSet', 'fet_mis',                                        ...
                 'tag_preprocessing', '+seh+',                               ...
                 'tag_postprocessing','_PPV2RGS',                              ...
                 'sampleRate', 12,                                           ...
                 'nNeurons',   100,                                          ...
                 'nIter',      100,                                          ...
                 'states',    {{'loc','rear','pause','groom','sit'}},        ...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',       true,                                         ...
                 'mref',       true,                                         ...
                 'prctTrain', []                                             ...                 
);
defargs = struct2varargin(defargs);
optimize_stc_transition_mis_reduced(defargs{:});




defargs = struct('fetSet',   'fet_mis',                                      ...
                 'mode',     'display',                                      ...
                 'tag_preprocessing',   '+seh+',                             ...
                 'tag_postprocessing',  '_PPV2RGS',                            ...
                 'rlist',    'training_hand_labeled_reduced',                ...
                 'slist',     {{'hand_labeled_jg_reduced';'hand_labeled_Ed_reduced'}},       ...
                 'sampleRate',12,                                            ...
                 'nNeurons',  100,                                           ...
                 'nIter',     100,                                           ...
                 'states',    {states},        ...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',      true,                                          ...
                 'mref',      true,                                          ...
                 'prctTrain', []                                             ...
);
defargs = struct2varargin(defargs);
req20160128(defargs{:});



ds = load('/storage/gravio/data/project/general/analysis/hand_labeled_Ed_reduced-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hl_3_jg_r_NN_100_NI_100_NN_multiPN_RAND_WSBNT_PPV2R-map2ref.mat');

ds = load('/storage/gravio/data/project/general/analysis/hand_labeled_jg_reduced-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hl_1_Ed_r_NN_100_NI_100_NN_multiPN_RAND_WSBNT_PPV2R-map2ref.mat');

ds = load('/storage/gravio/data/project/general/analysis/hand_labeled_jg-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat');

df = dir('/storage/gravio/data/project/general/analysis/hand_labeled_jg-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat')


df = dir('/storage/gravio/data/project/general/analysis/hand_labeled_jg-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat')
ds = load('/storage/gravio/data/project/general/analysis/hand_labeled_jg-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat')

fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname);
ds = load(fname);


fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)


fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed_NN_100_NI_100_NN_multiPN_RAND_WSBNT_PPNM-map2ref.mat';
df = dir(fname)
ds = load(fname)




fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT_PPNM-map2ref.mat';
df = dir(fname)
ds = load(fname)

fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)

fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)


fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)

fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hand_labeled_rev1_Ed_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)



fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hl_1_Ed_r_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat';
df = dir(fname)
ds = load(fname)



fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed_reduced-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_Ed03-20140625.cof.all_STC_hl_1_Ed_r_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat'
df = dir(fname)
ds = load(fname)




fname = '/storage/gravio/data/project/general/analysis/hand_labeled_Ed_reduced-MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hl_3_jg_r_NN_100_NI_100_NN_multiPN_RAND_WSBNT-map2ref.mat'
df = dir(fname)
ds = load(fname)






