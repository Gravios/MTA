function [stc,d_state,ls,lsm,mdl,p_state] = gen_stc_mpn(Trial,varargin)

% Default Arguments 
defargs = struct('refFilebase', 'jg05-20120317.cof.all',...
                 'fetSet'     , 'fet_mis',...
                 'states'     , {{'walk','rear','turn','pause','groom','sit'}},...
                 'stcMode'    , '',...
                 'refStcMode' , '',...
                 'sampleRate' , 12,...
                 'norm'       , true,...
                 'map2ref'    , true,...
                 'nNeurons'   , 100,...
                 'nIter'      , 100,...
                 'rndMethod'  , 'WSBNT',...
                 'overwrite'  , false ...
);

                 

% Assign Default Arguments to Workspace
[refFilebase, fetSet, states, stcMode, refStcMode, sampleRate, ...
 norm, map2ref, nNeurons, nIter, rndMethod, overwrite] = ...
DefaultArgs(varargin,defargs,'--struct');

% Load Trial
Trial = MTATrial.validate(Trial);

% Load refStcMode
if isempty(refStcMode),
    refStcMode = MTATrial.validate(refFilebase).stc.mode;
end

% Initialize outputs
stc = Trial.stc.copy;
d_state = [];
ls      = [];
lsm     = [];
mdl     = [];
p_state = [];

% Construct Model Name 
model = ['MTAC_BATCH-' fetSet ...
         '_SR_' num2str(sampleRate) ...
         '_NORM_' num2str(norm) ...
         '_REF_' refFilebase ...
         '_STC_' refStcMode ...
         '_NN_' num2str(nNeurons) ...
         '_NI_' num2str(nIter) ...
         '_NN_multiPN_RAND_' rndMethod];

% Try to find existing stc
stcList = listFiles(Trial.name,model);


if isempty(stcList)||overwrite, % Comp new stc if none exists / overwrite
    mod.states     = states;
    mod.stcMode    = stcMode;
    mod.featureSet = fetSet;
    mod.model      = model;
    mod.sampleRate = sampleRate;
    mod.nNeurons   = nNeurons;
    mod.nIter      = nIter;
    mod.map2reference = map2ref;
    mod.normalize = norm;
    argin = struct2varargin(mod);
    [stc,d_state,ls,lsm,mdl,p_state] = ...
        bhv_nn_multi_patternnet(Trial,argin{:});

    stc.save(1);
    
elseif numel(stcList)==1&&iscell(stcList), % return 
    stc.load(Trial,stcList{1});
    
else
    error('MTA:utilities:gen_stc_mpn:tooManyStcExist');
end


