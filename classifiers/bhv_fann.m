function [Stc,d_state,Model_Information] = bhv_fann(Trial,varargin)
%function [Stc,d_state] = bhv_nn(Trial,varargin)
%
% varargin:
%
%   train:       logical, def - false
%
%   states:    cellArray, def - Trial.stc.list_state_attrib('label')
%
%   fet:          string, def - 'fet_tsne'
%
%   modelName:   string, def - ['MTAC_' Trial.stc.mode '_LGR']
%
%   display:     logical, def - true
%
%   other_state: logical, def - false
%
%   nNeurons:    numeric, def - 100
%
%   subset:    MTADepoch, def - []
%
%   subset:      numeric, def - []


% Constants
MODEL_TYPE = 'NN';
SAMPLE_RATE = 20;

% Load Trial and or check if Trial is of the MTASession class
Trial = MTATrial.validate(Trial);
assert(isa(Trial,'MTASession'),['MTA:classifiers:' mfilename ':Trial not found']);


defArgs = {... 
 ...
 ...           trainModel,
               false,                                  ...
 ...          
 ...           states
               Trial.stc.list_state_attrib('label'),   ...
 ...
 ...           stc
               Trial.stc.copy,   ...
 ...
 ...           feature
               {'fet_tsne',Trial,10},  ...
 ...           
 ...           modelName
               ['MTAC_' Trial.stc.mode '_' MODEL_TYPE],...
 ...
 ...           display
               true,                                   ...
 ...
 ...           other_state
               false,                                  ...
 ... 
 ...           nNeurons
               100,                                    ...
 ...
 ...           subset
               [],                                     ...
 ...
 ...           index
               []                                      ...
};



[trainModel,states,Stc,feature,modelName,...
 display,other_state,nNeurons,subset,index] = DefaultArgs(varargin,defArgs);

Stc = Stc.copy;
%Stc.states = Stc(states{:});

keys = Stc.list_state_attrib('key');    


% LOAD feature if feature is a feature name
% RESAMPLE to conserve memory during model fitting
% feature.resample(30); for now leave it to the input SR
% LOAD feature if feature is a feature name
if iscell(feature),
    feature = feval(feature{:});
elseif ischar(feature);
    try
        feature = feval(feature,Trial,SAMPLE_RATE);
    catch err
        error(err.msg)
    end
end
assert(isa(feature,'MTAData'),['MTA:classifiers:' mfilename ':Feature not found']);


nind = nniz(feature);

% Create model filename based on the model name and the feature name
% default modelName = ['MTAC_' Trial.stc.mode '_' MODEL_TYPE];
modelName = [modelName '-' feature.label '-model.net'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,modelName);


%% Get or Train LGR Model
if trainModel||~exist(model_loc,'file'),
    
    % create NxS array to store states as integers (nomial data) {S=numel(states)}
    [smat] = stc2mat(Stc,feature,states);
    
    % Create struct to store model meta-data
    Model_Information = struct(...
        'filename',              modelName,         ...
        'path',                  model_path,             ...
        'description',           '',                     ...
        'StcMode',               Stc.mode,         ...
        'StcFilename',           Stc.filename,     ...
        'Trial',                 Trial.filebase,         ...
        'state_labels',          {states},               ...
        'state_keys',            {keys});

    % Flag to include unlabeled regions of data within the good
    % periods of the session. 
    % IF TRUE  -> Create classifier model with the specified states and
    %             all other states as a composite state
    if other_state, 
        ind = resample(Stc{'a'}.cast('TimeSeries'),feature);
        ind = logical(ind.data);
        smat = cat(2,smat,ind);
        smat(all(smat,2),end) = 0;
        if sum(smat(:,end))>0,
            Model_Information.state_labels =  [Model_Information.state_labels{:}, {'other'}];
            Model_Information.state_keys   =  [Model_Information.state_keys{:},       {'o'}];
        end
    else % IF FALSE -> Create classifier model with only the specified states
        ind = any(smat,2)&nniz(feature);
    end
    
    
    if ~isempty(subset),
        if isa(subset,'MTADepoch'),
            subset.resample(feature);
            subset.cast('TimeSeries',feature);
        end
        ind = ind&logical(subset.data);
    end
        

    fetMat = feature(ind,:);
    stsMat = ~~smat(ind,:);
    
    tmpfile_trn = fullfile(tempdir,['bhv_fann_train',num2str(randi(1e15,1,1)),'.dat']);
    fid = fopen(tmpfile_trn,'w+');
    % print header values
    fprintf(fid,'%i %i %i \n',size(fetMat,1),size(fetMat,2),size(stsMat,2));
    % print over loop input answer pairs
    for ind = 1:size(fetMat,1),
        fprintf(fid,[repmat('%12.12f ',[1,size(fetMat,2)]),'\n'],fetMat(ind,:));
        fprintf(fid,[repmat('%i ',[1,size(stsMat,2)]),'\n'],stsMat(ind,:));
    end
    fclose(fid);
    system(sprintf('bhv_fann_train %s %s %i %i %i %i',...
                   model_loc,tmpfile_trn,size(fetMat,2),size(stsMat,2),3,nNeurons));
    delete(tmpfile_trn);

    save([model_loc,'.mat'],'model_loc','Model_Information');
    return
else
    % Load the classifier model
    load([model_loc,'.mat']);
end    


% Create new StateCollection ... well copy

Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};

% Used to put labels into xyz sampleRate
xyz = Trial.load('xyz');



% Compute scores for the neural network

tmpfile_dat = fullfile(tempdir,['bft_dat',num2str(randi(1e15,1,1)),'.dat']);
tmpfile_res = fullfile(tempdir,['bft_res',num2str(randi(1e15,1,1)),'.dat']);
fid = fopen(tmpfile_dat,'w+');
% print over loop input
for ind = 1:size(feature,1),
    fprintf(fid,[repmat('%12.12f ',[1,size(feature,2)]),'\n'],feature(ind,:));
end
fclose(fid);
system(sprintf('bhv_fann_test %s %s %s',model_loc,tmpfile_dat,tmpfile_res));
d_state = load(tmpfile_res,'-ASCII');

% Clean up tmp files
delete(tmpfile_dat);
delete(tmpfile_res);


d_state = MTADxyz('data',d_state,'sampleRate',feature.sampleRate);
d_state.resample(xyz);
d_state = d_state.data;

% Separate the winners from the losers
[~,maxState] = max(d_state,[],2);
maxState(~nind,:) = 0;


% Smooth decision boundaries - 200 ms state minimum
bwin = round(.2*xyz.sampleRate)+double(mod(round(.2*xyz.sampleRate),2)==0);
mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
maxState=circshift(sq(mode(mss))',floor(bwin/2));

% Populate Stc object with the new states
for i = 1:numel(Model_Information.state_labels),

    sts = ThreshCross(maxState==i,0.5,1);
    if ~isempty(sts),
        sts = bsxfun(@plus,sts,[1,0]);
    end
    
    Stc.addState(Trial.spath,...
             Trial.filebase,...
             sts,...
             xyz.sampleRate,...
             feature.sync.copy,...
             feature.origin,...
             Model_Information.state_labels{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
end

%Stc.save(1);
