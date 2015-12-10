function [Stc,d_state] = bhv_lgr(Trial,varargin)
%function [Stc,d_state] = bhv_lgr(Trial,varargin)
%
% varargin:
%
%   train:       logical, def - false
%
%   states:    cellArray, def - Trial.stc.list_state_attrib('label')
%
%   fet:          string, def - 'fet_lgr'
%
%   model_name:   string, def - ['MTAC_' Trial.stc.mode '_LGR']
%
%   display:     logical, def - true
%
%
%
%

MODEL_TYPE = 'LGR';

[train,states,fet,model_name,display,other_state] = DefaultArgs(varargin,...
    {false,Trial.stc.list_state_attrib('label'),'fet_lgr',...
    ['MTAC_' Trial.stc.mode '_' MODEL_TYPE],true,false});



keys = subsref(Trial.stc.list_state_attrib('key'),...
               substruct('()',{Trial.stc.gsi(states)}));

% LOAD fet if fet is a feature name
% RESAMPLE to conserve memory during model fitting
% lrfet.resample(30); for now leave it to the input SR
if ischar(fet),
    lrfet = feval(fet,Trial,30);
else
    lrfet = fet;    
end

nind = nniz(lrfet);

% Create model filename based on the model name and the feature name
% default model_name = ['MTAC_' Trial.stc.mode '_' MODEL_TYPE];
model_name = [model_name '-' lrfet.label '-model.mat'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,model_name);


%% Get or Train LGR Model
if train||~exist(model_loc,'file'),
    
    % create Nx1 array to store states as integers (nomial data)
    [smat] = max(stc2mat(Trial.stc,lrfet,states),[],2);
    
    % Create struct to store model meta-data
    Model_Information = struct(...
        'filename',              model_name,         ...
        'path',                  model_path,             ...
        'description',           '',                     ...
        'StcMode',               Trial.stc.mode,         ...
        'StcFilename',           Trial.stc.filename,     ...
        'Trial',                 Trial.filebase,         ...
        'state_labels',          {states},               ...
        'state_keys',            {keys});

    % Flag to include unlabeled regions of data within the good
    % periods of the session. 
    % IF TRUE  -> Create classifier model with the specified states and
    %             all other states as a composite state
    if other_state, 
        ind = resample(Trial.stc{'a'}.cast('TimeSeries'),lrfet);
        ind = logical(ind.data);
        smat(smat==0) = numel(states)+1;
        if sum(smat(ind)==1)>0,
            Model_Information.state_labels =  [Model_Information.state_labels{:}, {'other'}];
            Model_Information.state_keys   =  [Model_Information.state_keys{:},       {'o'}];
        end
    else % IF FALSE -> Create classifier model with only the specified states
        ind = any(smat,2);
    end
    
    % Train classifier
    [B,dev,stats] = mnrfit(lrfet(ind,:),smat(ind),'model','nominal');
    save(model_loc,'B','dev','stats','Model_Information');
    return
else
    % Load the classifier model
    load(model_loc);
end    


% Create new StateCollection ... well copy
Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode...
                '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};

% Used to put labels into xyz sampleRate
xyz = Trial.load('xyz');

% Compute scores for the logistic regression
% I know this is irresposible code... don't give me that
% look... that was meant for me not you. It's really not as bad as
% it sounds.
d_state = mnrval(B,lrfet.data);
d_state = MTADxyz('data',d_state,'sampleRate',lrfet.sampleRate);
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
Stc.addState(Trial.spath,...
             Trial.filebase,...
             bsxfun(@plus,ThreshCross(maxState==i,0.5,1),[1,0]),...
             xyz.sampleRate,...
             lrfet.sync.copy,...
             lrfet.origin,...
             Model_Information.state_labels{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
%Stc.states{end} = Stc.states{end}+[1/lrfet.sampleRate,0];
end

Stc.save(1);