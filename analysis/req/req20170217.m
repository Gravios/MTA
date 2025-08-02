

%% Labeling statistics 
% Req 1.0: stats between neural network classifier (NNC) labels
% Req 1.1: stats between expert labeler classifier (ELC) labels 
% Req 1.2: stats between expert labeler classifier (ELC) labels 
%          and neural network classifier (NNC) labels 

% Req 1.0
% Comparison of neural network labels from NN's train on different
% animal with different labelers.
Trial = MTATrial.validate('Ed03-20140624.cof.all');
Trial = MTATrial.validate('Ed03-20140625.cof.all');

label_errors(Trial);

labelBhv_NN(Trial,'NN0317','jg05-20120317.cof.all','hand_labeled_rev3_jg');
labelBhv_NN(Trial,'NN0529','Ed05-20140529.ont.all','hand_labeled_rev1_Ed');

stc_jg = Trial.load('stc','NN0317');
stc_Ed = Trial.load('stc','NN0529');

labelStatsNNC = cmp_stcs(Trial,stc_jg,stc_Ed,[],true);

% Req 1.1 
% Comparision of expert labeler labels within same animal

Trial = MTATrial.validate('Ed01-20140707.cof.all');
stc_jg = Trial.load('stc','hand_labeled_rev2_jg');
stc_Ed = Trial.load('stc','hand_labeled_rev2_Ed');
labelStatsELC = cmp_stcs(Trial,stc_jg,stc_Ed,[],true);

Trial = MTATrial.validate('Ed03-20140624.cof.all');
stc_jg = Trial.load('stc','hand_labeled_rev1_jg');
stc_Ed = Trial.load('stc','hand_labeled_rev1_Ed');
labelStatsELC = cmp_stcs(Trial,stc_jg,stc_Ed,[],true);

% Req 1.2
Trial = MTATrial.validate('Ed03-20140624.cof.all');
stc_jg = Trial.load('stc','hand_labeled_rev1_jg');
stc_NN = Trial.load('stc','NN0317');
labelStatsNE = cmp_stcs(Trial,stc_jg,stc_NN);

% Req 1.3
Trial = MTATrial.validate('Ed03-20140625.cof.all');
Trial = MTATrial.validate('Ed03-20140624.cof.all');
Trial = MTATrial.validate('jg05-20120317.cof.all');
labelBhv_NN(Trial,'frNN0317','jg05-20120317.cof.all','hand_labeled_rev3_jg','fet_raw');
labelBhv_NN(Trial,'frNN0625','Ed03-20140625.cof.all','hand_labeled_rev1_Ed','fet_raw');
stc_NN = Trial.load('stc','NN0317');
stc_NN = Trial.load('stc','frNN0317');
stc_NN = Trial.load('stc','frNN0625');
stc_hl = Trial.load('stc','hand_labeled_rev1_Ed');
labelStatsNE = cmp_stcs(Trial,stc_hl,stc_NN);


