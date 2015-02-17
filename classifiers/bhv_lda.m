function [Stc,d_state] = bhv_lda(Trial,varargin)
%function [Stc,d_state] = bhv_lda(Trial,varargin)
%
% varargin:
%
%   train:       logical, def - false
%
%   states:    cellArray, def - Trial.stc.list_state_attrib('label')
%
%   fet:          string, def - 'fet_lda'
%
%   model_name:   string, def - ['MTAC_' Trial.stc.mode '_LDA']
%
%   display:     logical, def - true
%
%
%
%

MODEL_TYPE = 'LDA';

%Default Args
[train,states,fet,model_name,display] = DefaultArgs(varargin,...
    {false,Trial.stc.list_state_attrib('label'),'fet_lda',...
     ['MTAC_' Trial.stc.mode '_' MODEL_TYPE],true});


sind = Trial.stc.gsi(states);
keys = subsref(Trial.stc.list_state_attrib('key'),...
               substruct('()',{sind}));
Trial.stc.states = Trial.stc.states(sind);
ns = numel(states);

model_name = [model_name '-' fet '-model.mat'];
model_path = fileparts(mfilename('fullpath'));
model_loc = fullfile(model_path,model_name);

fet = feval(fet,Trial);
nind = nniz(fet);

%% Get or Train QDA Model
if train||~exist(model_loc,'file'),

    Model_Information = struct(...
        'filename',              model_name,         ...
        'path',                  model_path,             ...
        'description',           '',                     ...
        'StcMode',               Trial.stc.mode,         ...
        'StcFilename',           Trial.stc.filename,     ...
        'Trial',                 Trial.filebase,         ...
        'state_labels',          {states},               ...
        'state_keys',            {keys});

    fet_mean_state = zeros([1,size(fet,2),ns]);
    for i = 1:ns
        fet_state = fet(Trial.stc{states{i}},:);
        fet_state(~nniz(fet_state),:) = [];
        fet_mean_state(1,:,i) = nanmean(fet_state);
        fet_cov_state(:,:,i) = cov(fet_state);
    end
    save(model_loc,'fet_mean_state','fet_cov_state','Model_Information');
    return
else
    load(model_name);
end


Stc = Trial.stc.copy;
Stc.updateMode([MODEL_TYPE '_' Model_Information.StcMode '-' cell2mat(Model_Information.state_keys)]);
Stc.states = {};


states = Model_Information.state_labels;
ns = numel(states);

mean_fet_state = repmat(fet_mean_state,[fet.size(1),1,ns]);


%% Transform Features to QDA scores
d_state = zeros(fet.size(1),ns);
for i =  1:ns
    d_state(nind,i) = -.5*log(det(fet_cov_state(:,:,i)))...
        -.5*dot(((fet.data(nind,:)...
        -mean_fet_state(nind,:,i))/fet_cov_state(:,:,i))',(fet.data(nind,:)...
        -mean_fet_state(nind,:,i))');
end


d_state = Filter0(gtwin(0.50,30),d_state);

[~,maxState] = max(d_state,[],2);


% Push the new labels into the MTAStateCollection.
for i = 1:ns,
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,1),...
             fet.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             states{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
end


%keyboard
if display
sts_colors = 'brcmgky';
d_colors   = 'brcmgky';
keys = Model_Information.state_keys;
figure
hold on
for i = 1:ns,
plot(Filter0(gtwin(0.50,30),d_state(:,i)),d_colors(i)),
%Lines(Trial.stc{keys{i},Trial.xyz.sampleRate}(:),[],sts_colors(i));
end
Lines([],0,'k')
end