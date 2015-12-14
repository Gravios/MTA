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

        
        


        
        
