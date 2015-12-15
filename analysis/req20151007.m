% req20151007
% logistic regression models 
%
% label other sessions with hand labeled states, using the lgr model
% trained on the data set from jg05-20120317, and comapare with
% hand labels.
%



% Train a logistic regresion 
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
bhv_lgr(Trial,                                               ... Trial
        true,                                                ... ifTrain
        {'walk','rear','turn','pause','groom','sit'},        ... States
        'fet_tsne',                                          ... feature set
        'JG05HLR2');                                           % model name


% Generate a state collection based on a logistic regresion model
% train on the hand labeled data of jg05-20120317
GoThroughTrials('req20151009',                               ... SessionList
                @bhv_lgr,                                    ... Method
                false                                       ,... ifTrain
                {'walk','rear','turn','pause','groom','sit'},... States
                'fet_tsne',                                  ... feature set
                'JG05HLR2');                                   % model name

                
% tSNE dimesionallity reduction visuallization of clustering within
% the feature space and an overlay of behavior labeled with
% logistic regression 
GoThroughTrials('ncpH5B4',@req20151007_A1,[],[],{'walk','rear','turn','pause','groom','sit'},'LGR-hand_labeled_rev2-wrnpms');




% Net feature set fet_all
% Train a logistic regresion 
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
bhv_lgr(Trial,                                               ... Trial
        true,                                                ... ifTrain
        {'walk','rear','turn','pause','groom','sit'},        ... States
        'fet_all',                                           ... feature set
        'JG05HLR2FALL');                                           % model name

        
        

% Check that feature transform based on external reference shares
% similar distribution
fet_jg = fet_tsne(MTATrial('jg05-20120317'),10,false); 
fet_edo = fet_tsne(MTATrial('Ed01-20140707'),10,false);
fet_edr = fet20151007(MTATrial('Ed01-20140707'),10);

f = 2;
eds = linspace(-4,200,500);
figure;hold on;
ho = bar(eds,hist(fet_edo(:,f),eds),'histc');
ho.FaceColor = 'k';
ho.FaceAlpha = .6;
ho.EdgeColor = 'k';
ho.EdgeAlpha = .6;
hs = bar(eds,hist(fet_jg(:,f),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .6;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .6;
hr = bar(eds,hist(fet_edr(:,f),eds),'histc');
hr.FaceColor = 'c';
hr.FaceAlpha = .6;
hr.EdgeColor = 'c';
hr.EdgeAlpha = .6;


STATES = {'walk','rear','turn','pause','groom','sit'};

% Net feature set fet_all
% Train a logistic regresion model using features of jg05-20120317
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
bhv_lgr(Trial,                                               ... Trial
        true,                                                ... ifTrain
        STATES,                                              ... States
        {'fet_tsne',Trial,20,true},                          ... feature set
        'UFET20151007REFJG0520120317');                                    % model name


% Label Ed01-20140707 using model REFJG0520120317
Trial = MTATrial('Ed01-20140707');
[Stc,dsc] = bhv_lgr(Trial,                                  ... Trial
            false,                                          ... ifTrain
            STATES,                                         ... States
            {'fet20151007',Trial,20,true},                       ... feature set
            ...{'fet_tsne',Trial,20},                       ... feature set
            'UFET20151007REFJG0520120317');                                  %
                                                                 % model 
StcHL = Trial.load('stc','hand_labeled_rev1');
            
            
            


% Confusion matrix
xyz = Trial.load('xyz');
shl = MTADxyz('data',double(0<stc2mat(StcHL,xyz,STATES)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc,  xyz,STATES)),'sampleRate',xyz.sampleRate); 
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
h.XTickLabel = cellfun(@horzcat,repmat({'     '},1,8),{STATES{:},'sensitivity',''},'uniformoutput',false)
tstates = {STATES{:},'precision',''};
tstates = cellfun(@horzcat,tstates,cellfun(@repmat,repmat({' '},1,numel(tstates)),repmat({1},1,numel(tstates)),mat2cell(abs(cellfun(@length,tstates)-max( cellfun(@length,tstates))),1,ones(1,numel(tstates))),'uniformoutput',false),'uniformoutput',false);
h.YTickLabel = fliplr(tstates);
suptitle({['neural network vs hand labels for ',Trial.name],'nn trained jg05-20120317'})




