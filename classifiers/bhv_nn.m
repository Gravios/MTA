function [Stc,d_state,Model_Information] = bhv_nn(Trial,varargin)
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
%

% Constants
MODEL_TYPE = 'NN';
SAMPLE_RATE = 20;

% Load Trial and or check if Trial is of the MTASession class
if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end
assert(isa(Trial,'MTASession'),['MTA:classifiers:' mfilename ':Trial not found']);


defArgs = {... 
 ...
 ...           trainModel,
               false,                                  ...
 ...          
 ...           states
               Trial.stc.list_state_attrib('label'),   ...
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
               []                                      ...
};



[trainModel,states,feature,modelName,display,other_state,nNeurons,subset] = DefaultArgs(varargin,defArgs);

if isa(states,'MTAStateCollection'),
    Stc = states.copy;
    states = Stc.list_state_attrib('label');
else
    Stc = Trial.stc.copy;
    Stc.states = Stc(states{:});
end

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
modelName = [modelName '-' feature.label '-model.mat'];
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
        ind = any(smat,2);
    end
    
    
    if ~isempty(subset),
        if isa(subset,'MTADepoch'),
            subset.resample(feature);
            subset.cast('TimeSeries',feature);
        end
        ind = ind&logical(subset.data);
    end
        
    % Train classifie
    net = patternnet(nNeurons);
    %view(net);    
    [net,tr] = train(net,feature(ind,:)',~~smat(ind,:)');

    save(model_loc,'net','tr','Model_Information');
    return
else
    % Load the classifier model
    load(model_loc);
end    


% Create new StateCollection ... well copy

Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};

% Used to put labels into xyz sampleRate
xyz = Trial.load('xyz');

% Compute scores for the neural network
d_state = net(feature.data')';

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