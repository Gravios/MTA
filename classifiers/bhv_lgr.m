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
if ischar(fet),
    lrfet = feval(fet,Trial);
else
    lrfet = fet;    
end

% RESAMPLE to conserve memory during model fitting
% lrfet.resample(30); for now leave it to the input SR
nind = nniz(lrfet);

model_name = [model_name '-' fet.label '-model.mat'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,model_name);



%% Get or Train LGR Model
if train||~exist(model_loc,'file'),

    [smat] = max(stc2mat(Trial.stc,lrfet,states),[],2);
    
    Model_Information = struct(...
        'filename',              model_name,         ...
        'path',                  model_path,             ...
        'description',           '',                     ...
        'StcMode',               Trial.stc.mode,         ...
        'StcFilename',           Trial.stc.filename,     ...
        'Trial',                 Trial.filebase,         ...
        'state_labels',          {states},               ...
        'state_keys',            {keys});

    if other_state, 
        ind = resample(Trial.stc{'a'}.cast('TimeSeries'),lrfet);
        ind = logical(ind.data);
        smat(smat==0) = numel(states)+1;
        if sum(smat(ind)==1)>0,
            Model_Information.state_labels =  [Model_Information.state_labels{:}, {'other'}];
            Model_Information.state_keys   =  [Model_Information.state_keys{:},       {'o'}];
        end
    else
        ind = any(smat,2);
    end
    
    
    
    [B,dev,stats] = mnrfit(lrfet(ind,:),smat(ind),'model','nominal');
    save(model_loc,'B','dev','stats','Model_Information');
    return
else
    load(model_loc);
end    


% Create new Stc
Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};


% Compute scores for the logistic regression
% You know this is irresposible code... don't give me that look
d_state = mnrval(B,lrfet.data);
d_state = MTADxyz('data',d_state,'sampleRate',lrfet.sampleRate);
d_state.resample(Trial.load('xyz'));
d_state = d_state.data;

% Separate the winners from the losers
[~,maxState] = max(d_state,[],2);
maxState(~nind,:) = 0;

% Smooth decision boundaries - 200 ms state minimum
bwin = round(.2*lrfet.sampleRate)+double(mod(round(.2*lrfet.sampleRate),2)==0);
mss = GetSegs(maxState,1:size(maxState,1),bwin,nan);
maxState=circshift(sq(mode(mss))',floor(bwin/2));

% Populate Stc object with the new states
for i = 1:numel(Model_Information.state_labels),
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,1),...
             lrfet.sampleRate,...
             lrfet.sync.copy,...
             lrfet.origin,...
             Model_Information.state_labels{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
Stc.states{i} = Stc.states{i}+[1/lrfet.sampleRate,0];
end


