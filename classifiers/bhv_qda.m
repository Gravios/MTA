function Stc = bhv_qda(Trial,Stc,varargin)
[train,states,model_filename,display] = DefaultArgs(varargin,{false,[],'MTA_manual_mknsrw_QDA_model.mat',true});


MODEL_TYPE = 'MTAC_QDA';


% If Trial is not a MTASession try loading it.
Trial = MTATrial.validate(Trial);

varargout = cell([1,nargout-1]);

% Default Arguments for uninitiated variables
defArgs = struct('states',      {{'walk','rear','turn','pause','groom','sit'}},...
                 'stcMode',     '',...
                 'featureSet',  'fet_tsne_rev3',...
                 'sampleRate',  12,...
                 'model',       [],...
                 'nNeurons',    100,...          
                 'nIter',       100,...
                 'randomizationMethod','',...
                 'map2reference', false,... 
                 'normalize',   false,...
                 'tag',         '',...
                 'targetState', '',...
                 'prctTrain',   []...
);


% LOAD the feature set
if isa(featureSet,'MTADfet'),
    features = featureSet.copy;
    featureSet = features.label;
else,
    if ~strcmp(featureSet,'fet_all'), % For the love of God remove this!
        features = feval(featureSet,Trial,sampleRate,false);
    end
end

% LOAD the hand labels
if ~isempty(stcMode),
    if ischar(stcMode),
        Trial.load('stc',stcMode);
        StcHL = Trial.stc.copy;
    elseif isa(stcMode,'MTAStateCollection'),
        StcHL = stcMode.copy;
        stcMode = stcMode.mode;
    end
    StcHL.states = StcHL(states{:});
    keys = cellfun(@subsref,StcHL(states{:}),repmat({substruct('.','key')},size(states)));
end


% If the model name is empty then create a composite name and 
% toggle train to create new neural network models
if isempty(model),

    if ~isempty(prctTrain), prctTrainTag = ['_PRT_',prctTrain]; else, prctTrainTag = ''; end
    
    model = ['MTAC_BATCH-' tag targetState featureSet ...
             '_SR_'  num2str(sampleRate)              ...
             '_NORM_' num2str(normalize)              ...             
             '_REF_' Trial.filebase                   ...
             '_STC_' stcMode                          ...
             '_NN_'  num2str(nNeurons)                ...
             '_NI_'  num2str(nIter)                   ...
             prctTrainTag                             ...
             '_' MODEL_TYPE];
    train = true;
end



%% Get or Train QDA Model
if train
    
    if isempty(states),   states = Trial.stc.list_state_attrib('label'); end
    ns = numel(states);

    Model_Information.description = '';
    Model_Information.StcMode     = Trial.stc.mode;
    Model_Information.StcFilename = Trial.stc.filename;
    Model_Information.SessionFilebase = Trial.filebase;
    Model_Information.state_labels = states;
    keys = Trial.stc.list_state_attrib('key');
    fet_mean_state = zeros([1,size(fet,2),ns]);
    for i = 1:ns
        fet_state = fet(Trial.stc{states{i}},:);
        fet_state(sum(isnan(fet_state)|...
                      isinf(fet_state),2)>0,:) = [];
        fet_mean_state(1,:,i) = nanmean(fet_state);
        cov_state(:,:,i) = cov(fet_state);
        sti = Trial.stc.gsi(states{i});
        if ~isempty(sti)
            Model_Information.state_keys(i) = keys(sti);
        else
            spare_keys = 'qjpxothrymunslf';
            spare_keys = spare_keys(~ismember(spare_keys,cell2mat(keys)));
            Model_Information.state_keys(i) = spare_keys(i);
        end
    end
    save(fullfile(fileparts(mfilename('fullpath')),model_filename),...
         'fet_state','fet_mean_state','cov_state','Model_Information');
    return
end

load(model_filename);

states = Model_Information.state_labels;
ns = numel(states);

mean_fet_state = repmat(fet_mean_state,[fet.size(1),1,ns]);




%% Transform Features to QDA scores




d_state = zeros(fet.size(1),ns);
for i =  1:ns
    d_state(fet_not_nan,i) = -.5*log(det(cov_state(:,:,i)))...
        -.5*dot(((fet.data(fet_not_nan,:)...
        -mean_fet_state(fet_not_nan,:,i))/cov_state(:,:,i))',(fet.data(fet_not_nan,:)...
        -mean_fet_state(fet_not_nan,:,i))');
end



d_state = Filter0(dwin,d_state);

[~,maxState] = max(d_state,[],2);

for i = 1:ns,
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,10),...
             Trial.xyz.sampleRate,...
             Trial.xyz.sync.copy,...
             Trial.xyz.origin,...
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
plot(Filter0(dwin,d_state(:,i)),d_colors(i)),
%Lines(Trial.stc{keys{i},Trial.xyz.sampleRate}(:),[],sts_colors(i));
end
Lines([],0,'k')
end