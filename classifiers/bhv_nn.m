function [Stc,d_state,Model_Information,n_state] = bhv_nn(Trial,varargin)
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
%   modelName:    string, def - ['MTAC_' Trial.stc.mode '_LGR']
%
%   display:     logical, def - true
%
%   otherState:  logical, def - false
%
%   nNeurons:    numeric, def - 100
%
%   subset:    MTADepoch, def - []
%
%   subset:      numeric, def - []


% Constants
MODEL_TYPE = 'NN';

% Load Trial and or check if Trial is of the MTASession class
if iscell(Trial)
    
else
    Trial = MTATrial.validate(Trial);
end



if isa(Trial,'MTASession'),
    stcMode = Trial.stc.mode;
else
    stcMode = '';
end


defArgs = struct('trainModel',    false,...
                 'states',        {Trial.stc.list_state_attrib('label')},...
                 'stc',           Trial.stc.copy,   ...
                 'feature',       [],  ...
                 'modelName',     ['MTAC_' stcMode '_' MODEL_TYPE],...
                 'display',       true,...
                 'otherState',   false,...
                 'nNeurons',      100,...
                 'subset',        [],...
                 'index',         []...
);



[trainModel,states,Stc,feature,modelName,...
 display,otherState,nNeurons,subset,index] = DefaultArgs(varargin,defArgs,'--struct');

Stc = Stc.copy;
%Stc.states = Stc(states{:});


assert(isa(feature,'MTAData'),['MTA:classifiers:' mfilename ':Feature not found']);


nind = nniz(feature);

% Create model filename based on the model name and the feature name
% default modelName = ['MTAC_' Trial.stc.mode '_' MODEL_TYPE];
modelName = [modelName '-' feature.label '-model.mat'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,modelName);


%% Get or Train NN Model
if trainModel||~exist(model_loc,'file'),

    keys = Stc.list_state_attrib('key');
    
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
    if otherState, 
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
        
    % Train classifier
    net = patternnet(nNeurons);
    %net.trainParam.showWindow = true;
    net.trainParam.showWindow = false;
    
% SET mapminmax with a Session group
% $$$     net = struct(net);
% $$$     netInputMapminmaxIndex = ~cellfun(@isempty,regexp(net.inputs{1}.processFcns,'mapminmax'));
% $$$     net.inputs{1}.processSettings{netInputMapminmaxIndex}  = ...
% $$$         load_normalization_parameters_mapminmax(feature.label,...
% $$$                                                 [],... need feature.treatmentRecord
% $$$                                                 'hand_labeled');
% $$$     net = network(net);                                            

    psa = load_normalization_parameters_mapminmax(feature.label,...
                                                 [],... need feature.treatmentRecord
                                                 'hand_labeled');
    %view(net);        
    %[net,tr] = train(net,feature(ind,:)',~~smat(ind,:)');
    net.inputs{1}.processFcns(2) = [];
    %mfet = ;
    [net,tr] = train(net,mapminmax('apply',feature(ind,:)',psa),~~smat(ind,:)');
% $$$     [net,tr] = train(net,...
% $$$                      [ones([1,psa.xrows]),-ones([1,psa.xrows]),mfet],...
% $$$                      [zeros([2,size(smat,2)]),~~smat(ind,:)']);

    save(model_loc,'net','tr','Model_Information');
    return
else
    % Load the classifier model
    load(model_loc);
end    


% Create new StateCollection ... well copy

Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys(~isempty(Model_Information.state_keys)))]);
Stc.states = {};

% Used to put labels into xyz sampleRate
xyz = Trial.load('xyz');

psa = load_normalization_parameters_mapminmax(feature.label,...
                                              [],... need feature.treatmentRecord
                                              'hand_labeled');

% Compute scores for the neural network
d_state = net(mapminmax('apply',feature.data',psa))';

d_state = MTADxyz('data',d_state,'sampleRate',feature.sampleRate);
d_state.resample(xyz);
d_state = d_state.data;
n_state = d_state;

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