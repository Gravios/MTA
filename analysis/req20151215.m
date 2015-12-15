function req20151215(Trial,varargin)
%req20151215
% Primary goal:
%   Compare feature correction methods base on the quality of
%   labels assigned by a logistic regression or neural network 
%   model using jg05-20120317 as a reference.
%
%   vars: features {fet_tsne,fet20151007}
%
%   Model is trained with feature 'fet_tsne'.   
%
%   Conditions:
%     1. fet_tsne w/o normalization 
%     2. fet_tsne w/ normalization
%     3. fet_fet20151007 w/o normalization
%     4. fet_fet20151007 w/ normalzation
%

defArgs = {...
    ... states 
        {'walk','rear','turn','pause','groom','sit'},...
    ... 
    ... normalize
        false,                                       ...
    ...
    ... sampleRate
        20,                                          ...
    ...
    ... modFet
        'fet_tsne',                                  ...
    ...
    ... labFet
        'fet_tsne',                                  ...
    ...

        
         
normalize = true;

states = 
SAVE_FIG = false;

modelName = 'FET20151007REFJG0520120317';

if normalize, 
    modelName = ['U',modelName];
end

% Net feature set fet_all
% Train a logistic regresion model using features of jg05-20120317
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
bhv_lgr(Trial,                                               ... Trial
        true,                                                ... ifTrain
        states,                                              ... States
        {modFet,Trial,sampleRate,normalize},                     ... feature set
        modelName);                         % model name


% Label Ed01-20140707 using model REFJG0520120317
Trial = MTATrial('Ed01-20140707');
[Stc,dsc] = bhv_lgr(Trial,                                  ... Trial
            false,                                          ... ifTrain
            states,                                         ... States
            {'fet20151007',Trial,20,normalize},             ... feature set
            modelName);                    % model 
StcHL = Trial.load('stc','hand_labeled_rev1');
            


% Confusion matrix
xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,states)),'sampleRate',xyz.sampleRate); 
cm = confmat(shl(:,:),ysm(:,:)); % DEP: netlab                
precision = round(diag(cm)'./sum(cm),4).*100;
sensitivity = round(diag(cm)./sum(cm,2),4).*100;
cm = round(cm./xyz.sampleRate,2);


% Plot confusion matrix
hfig = figure(3);clf
hfig.Units = 'normalized';
hfig.Position = [0.15,0.15,0.66,0.75];
uicontrol(hfig,...
          'Style','text',...
          'String','LGR Labeled',...
          'FontSize',14,...
          'Units','normalized',...
          'Position',[0.05,0.5,0.08,0.05]);
uicontrol(hfig,...
          'Style','text',...
          'String','Hand Labeled',...
          'FontSize',14,...
          'Units','normalized',...
          'Position',[0.5,0.01,0.08,0.05]);
plotmat([[cm,sensitivity];[precision,0]],'b','k',14);
h = gca;
h.XAxisLocation = 'top';
h.XTickLabel = cellfun(@horzcat,repmat({'     '},1,8),{states{:},'sensitivity',''},'uniformoutput',false)
tstates = {states{:},'precision',''};
tstates = cellfun(@horzcat,tstates,cellfun(@repmat,repmat({' '},1,numel(tstates)),repmat({1},1,numel(tstates)),mat2cell(abs(cellfun(@length,tstates)-max( cellfun(@length,tstates))),1,ones(1,numel(tstates))),'uniformoutput',false),'uniformoutput',false);
h.YTickLabel = fliplr(tstates);
suptitle({['Logistic Regression vs hand labels for ',Trial.name],'nn trained jg05-20120317'})


reportfig([], hfig, 'confmats', 'req', false,Trial.filebase,'',[],SAVE_FIG,'png');