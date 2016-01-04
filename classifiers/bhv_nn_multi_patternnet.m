function [stc,varargout] = bhv_nn_multi_patternnet(Trial,varargin)
%function [stc,varargout] = bhv_nn_multi_patternnet(Trial,varargin)
% 
%
% varargin:
%    states,
%    stcMode,
%    featureSet,
%    sampleRate,
%    model,
%    nNeurons,
%    nIter,
%    randomizationMode,
%


% If Trial is not a MTASession try loading it.
if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end


% Default Arguments for uninitiated variables
defArgs = {...
           ... states
               {'walk','rear','turn','pause','groom','sit'}             ...
           ...
           ... stcMode
               'hand_labeled_rev2',                                     ...
           ...
           ... featureSet 
               'fet_tsne_rev3',                                         ...
           ...
           ... sampleRate
               12,                                                      ...
           ...
           ... model
               [],                                                      ...
           ...
           ... nNeurons          
               100,                                                     ...
           ...
           ... nIter
               100,                                                     ...
           ...
           ... randomizationMode
               ''                                                       ...
};

[states,stcMode,featureSet,sampleRate,model,nNeurons,nIter,randomizationMode] = DefaultArgs(varargin,defArgs);


% Load the feature set
features = feval(featureSet,Trial,sampleRate,false);

% Default mode is labeling
train = false;

% If the model name is empty then create a composite name and 
% toggle train to create new neural network models
if isempty(model),
    model = ['MTAC_BATCH-' featureSet ...
             '_SR_'  num2str(sampleRate) ...
             '_REF_' Trial.filebase ...
             '_NN_'  num2str(nNeurons) ];
    %model = 'fet_tsne_REFjg0520120317_NN';
    train = true;
end


% Use xyz as base for labeling output ... quit staring at this part.
xyz = Trial.load('xyz');

% Initialize outputs
tysm = zeros([xyz.size(1),numel(states)]);
labelingStats.confussionMatrix = zeros([nIter,numel(states),numel(states)]);
labelingStats.precision =   zeros([nIter,numel(states)]);
labelingStats.sensitivity = zeros([nIter,numel(states)]);
labelingStats.accuracy =    zeros([nIter,1]);

% Initialize labeling timeperiods
labelingEpochs = Trial.stc{'a'}.cast('TimeSeries');

% Load the hand label
if ~isempty(stcMode),
    Trial.load('stc',stcMode);
    StcHL = Trial.stc.copy;
end

for iter = 1:nIter,
    try,
        if train,
            % Create two partitions for validation
            %blocs  = reshape(1mod(features.size(1), 4.*20 )
            rndInd = randperm(features.size(1))';
            rndInd = rndInd(1:floor(features.size(1)/2));
            trndInd = false([features.size(1),1]);
            trndInd(rndInd) = true;
            rndInd = trndInd;


            trainingEpochs = MTADepoch([],[],              ...
                                       rndInd,             ...
                                       features.sampleRate,... 
                                       features.sync.copy, ... 
                                       features.origin,    ... 
                                       'TimeSeries',       ...
                                       [],[],              ...
                                       'training',         ... label
            't'                ... key
            );

            labelingEpochs = MTADepoch([],[],...
                                       ~rndInd,...
                                       features.sampleRate,... 
                                       features.sync.copy, ... 
                                       features.origin,    ... 
                                       'TimeSeries',       ...
                                       [],[],              ...
                                       'labeling',         ... label
            'l'                ... key
            );



            % Train Model
            bhv_nn (Trial,                                           ... Trial
            true,                                            ... ifTrain
            states,                                          ... States
            features,                                        ... feature set
            [model '_' num2str(iter)],                                           ... model name
            'subset',trainingEpochs);                                        
        

        end
        % Label States
        Stc = bhv_nn (Trial,                                         ... Trial
        false,                                         ... ifTrain
        states,                                        ... States
        features,                                      ... feature set
        [model '_' num2str(iter)]);                                          % model name

        


        ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
        tysm = ysm.data+tysm;
        
        if ~isempty(stcMode),
            shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);

            labelingEpochs.resample(xyz);

            ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;

            tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
            
            labelingStats.precision(iter,:) = round(diag(tcm)'./sum(tcm),4).*100;
            labelingStats.sensitivity(iter,:) = round(diag(tcm)./sum(tcm,2),4).*100;
            labelingStats.accuracy(iter) = sum(diag(tcm))/sum(tcm(:));
            labelingStats.confusionMatrix(iter,:,:) = round(tcm./xyz.sampleRate,2);        
        end
    catch
        continue
    end
    
end



% $$$ 
% $$$ figure,imagesc(tysm')
% $$$ Lines(StcHL{'w'}(:),[],'m');
% $$$ Lines(StcHL{'n'}(:),[],'g');
% $$$ Lines(StcHL{'r'}(:),[],'r');
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ if train,
% $$$     cm_jg05          = cm;
% $$$     precision_jg05   = precision;
% $$$     sensitivity_jg05 = sensitivity;
% $$$     accuracy_jg05    = accuracy;
% $$$ end
% $$$ 
% $$$ save(['/storage/gravio/data/project/general/analysis/req20151216-' ...
% $$$      model '-' Trial.filebase '.mat']);
% $$$ 
% $$$ load(['/storage/gravio/data/project/general/analysis/req20151216-' ...
% $$$      model '-' Trial.filebase '.mat']);
% $$$ 
% $$$ 
% $$$ tysm = ysm.data+tysm;
% $$$ 
% $$$ [~,sind] = max(tysm,[],2);
% $$$ 
% $$$ ysm = MTADxyz('data',zeros([shl.size]),'sampleRate',xyz.sampleRate); 
% $$$ ysm.data = ysm.data';
% $$$ ysm.data([1:6:size(ysm,2).*6]+sind'-1) = 1;
% $$$ ysm.data = ysm.data';
% $$$ 
% $$$ ind = any(shl.data,2)&any(ysm.data,2)&labelingEpochs.data;
% $$$ tcm = confmat(shl(ind&labelingEpochs,:),ysm(ind&labelingEpochs,:)); % DEP: netlab
% $$$ nprecision = round(diag(tcm)'./sum(tcm),4).*100;
% $$$ nsensitivity = round(diag(tcm)./sum(tcm,2),4).*100;
% $$$ ncm = round(tcm./xyz.sampleRate,2);
% $$$ naccuracy = sum(diag(tcm))/sum(tcm(:));
% $$$ 
% $$$ 
% $$$ % Plot the distribution of model precision
% $$$ plot_precision(precision,states,'fet_tsne')
% $$$ 
% $$$ 
