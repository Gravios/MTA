
Trial = MTATrial.validate('Ed01-20140707.cof.all');

stc_Ed = Trial.load('stc','hand_labeled_rev2_Ed');
stc_jg = Trial.load('stc','hand_labeled_rev2_jg');
states = {'walk','rear','turn','pause','groom','sit','shake'};

lstats = compute_inter_stc_stats(Trial,stc_Ed,stc_jg,states,119.880135);

% $$$ confutionMatrix:
% $$$ 316.29       0.2     15.52     39.26         0      0.06      1.39
% $$$   1.03     49.62      0.02      2.63      0.15         0         0
% $$$   7.46      0.33     52.68      7.93         0         0         0
% $$$  27.79      1.27     14.22    378.25      2.09     15.52      1.71
% $$$      0         0         0      6.99    181.34         0         0
% $$$   0.53      0.18         0      8.71         0    207.32         0
% $$$   0.78         0         0       0.1         0         0       1.4
% $$$ 
% $$$       precision: [84.86 92.85 77.02 85.8 96.29 95.65 61.31]
% $$$     sensitivity: [89.38 96.17 63.9 85.22 98.78 93.01 31.11]
% $$$        accuracy: 0.883917699972666

Trial = MTATrial.validate('Ed03-20140624.cof.all');

stc_Ed = Trial.load('stc','hand_labeled_rev1_Ed');
stc_jg = Trial.load('stc','hand_labeled_rev1_jg');
states = {'walk','rear','turn','pause','groom','sit','shake'};

lstats = compute_inter_stc_stats(Trial,stc_Ed,stc_jg,states,119.880135);

% $$$ confusionMatrix:
% $$$  131.29      0.14     10.44     19.77       0.1      0.96      0.26
% $$$    0.21     14.92      0.19      0.63      0.33         0         0
% $$$   15.16      0.03     57.21     13.14         0      0.18      0.01
% $$$   16.27      1.58     34.88    493.53      5.93     77.88      0.08
% $$$    0.63      0.21         0      12.2    123.23         0       0.1
% $$$    0.06         0      0.97     21.26         0    823.54         0
% $$$    0.06         0         0      0.13         0         0      1.26
% $$$ 
% $$$       precision: [80.57 91.65 66.74 78.32 90.37 97.36 87.28]
% $$$     sensitivity: [80.22 88.43 55.18 88.03 95.09 91.25 74.02]
% $$$        accuracy: 0.875588866047127



