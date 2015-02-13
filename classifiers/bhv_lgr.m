function [Stc,d_state] = bhv_lgr(Trial,varargin)
MODEL_TYPE = 'LGR';
[train,states,model_name,display] = DefaultArgs(varargin,...
    {false,Trial.stc.list_state_attrib('label'),...
    ['MTAC_' Trial.stc.mode '_' MODEL_TYPE],true});

keys = subsref(Trial.stc.list_state_attrib('key'),...
               substruct('()',{Trial.stc.gsi(states)}));


model_name = [model_name  '-model.mat'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,model_name);


lrfet = fet_lgr(Trial);
nind = nniz(lrfet);

%% Get or Train LGR Model
if train||~exist(model_loc,'file'),

    Model_Information = struct(...
        'filename',              model_filename,         ...
        'path',                  model_path,             ...
        'description',           '',                     ...
        'StcMode',               Trial.stc.mode,         ...
        'StcFilename',           Trial.stc.filename,     ...
        'Trial',                 Trial.filebase,         ...
        'state_labels',          {states},               ...
        'state_keys',            {keys});

    [smat] = max(stc2mat(Trial.stc,lrfet,states),[],2);
    ind = any(smat,2);
    B = mnrfit(lrfet(ind,:),smat(ind),'model','nominal');
    save(model_loc,'B','Model_Information');
    return
else
    load(model_loc);
end    


Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '-' Model_Information.StcMode '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};


d_state = mnrval(B,lrfet.data);

[~,maxState] = max(d_state,[],2);
maxState(~nind,:) = 0;

for i = 1:numel(Model_Information.state_labels),
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,10),...
             lrfet.sampleRate,...
             lrfet.sync.copy,...
             lrfet.origin,...
             Model_Information.state_labels{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
end


