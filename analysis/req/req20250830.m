

configure_default_args();
% >>> Project Vars >>> --------------------------------------------------------
%sessionlist = get_session_list_v3('BehaviorPlaceCode');
sessionlist = query_session_list('BehaviorPlaceCode');

sessionlist(14:22) = []; % remove jg04 for this analysis

samplerate = 250;% Hz

Trials = af(@(s) MTATrial.validate(s), sessionlist);

bhv_state = 'theta-groom-sit-rear-hpause';

% Ego field bins
xbins = 36;
ybins = 36;
xminmax = [-500,500]; % unused at the moment;

bhv_state = 'theta-groom-sit-rear-hpause';

swfN = bins.swf.count;
hbaN = bins.hba.count;
hvlN = bins.hvl.count;
phzN = bins.phz.count;

sfob = 1;

% >>> PYRAMIDAL CELLS >>> -----------------------------------------------------
otid = tid;
unit = {};
% >>> PYR ER06 >>> ------------------

% >>> PYR ER06-20130612 >>> ---------

tid = find_trial_index(sessionlist,'ER06-20130612.cof.all');

unit{tid}(1).id = 132;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,12,15];

unit{tid}(end+1).id = 151;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [3,6,10,13,16];

unit{tid}(end+1).id = 158;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [3, 4, 6, 8, 10, 11, 12];

unit{tid}(end+1).id = 173;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,3,4,6,10];

% <<< PYR ER06-20130612 <<< ---------
% >>> PYR ER06-20130613 >>> ---------
tid = find_trial_index(sessionlist,'ER06-20130613.cof.all');

unit{tid}(1).id = 4;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 21;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 59;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 65;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 69;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 123;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 131;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 132;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];


% <<< PYR ER06-20130613 <<< ---------

% <<< PYR ER06 >>> ------------------
% >>> PYR Ed10 >>> ------------------
% >>>   PYR Ed10-20140814 >>> ---------

tid = find_trial_index(sessionlist,'Ed10-20140814.cof.all');

unit{tid}(1).id = 7;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,16];

unit{tid}(end+1).id = 39;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10];

unit{tid}(end+1).id = 42;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10,13];

unit{tid}(end+1).id = 43;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10];

unit{tid}(end+1).id = 44;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,6,7,8,10,12];

unit{tid}(end+1).id = 70;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10];

unit{tid}(end+1).id = 77; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,18];

unit{tid}(end+1).id = 95;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22,23,24];

unit{tid}(end+1).id = 103; % STAR;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,9,10,12,13,15];

unit{tid}(end+1).id = 105;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,15,18];

unit{tid}(end+1).id = 108;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,11];

unit{tid}(end+1).id = 111; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,21,22,24]; % MEH

unit{tid}(end+1).id = 112;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,22,24]; 

unit{tid}(end+1).id = 113;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,21,22,24];

unit{tid}(end+1).id = 121;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,10,13,16,19];

unit{tid}(end+1).id = 124;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,8,10,13,16,19];



% <<< PYR Ed10-20140814 <<< ---------
% >>>   PYR Ed10-20140815 >>> ---------

tid = find_trial_index(sessionlist,'Ed10-20140815.cof.all');

unit{tid}(1).id = 21;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19];

unit{tid}(end+1).id = 37;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4];

unit{tid}(end+1).id = 44;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,3,4,7];

unit{tid}(end+1).id = 47; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,14,16,17];

unit{tid}(end+1).id = 48; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,3];

unit{tid}(end+1).id = 51; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10,18];

unit{tid}(end+1).id = 77;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,18];

unit{tid}(end+1).id = 86;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19];

unit{tid}(end+1).id = 92; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,6,10];

unit{tid}(end+1).id = 94; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,11];

unit{tid}(end+1).id = 96;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,10];


% <<< PYR Ed10-20140815 >>> ---------
% >>>   PYR Ed10-20140816 >>> ---------

tid = find_trial_index(sessionlist,'Ed10-20140816.cof.all');

unit{tid}(1).id = 7; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,18,21];

unit{tid}(end+1).id = 10; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [5,7,10,11,16,21];

unit{tid}(end+1).id = 13;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19,22];

unit{tid}(end+1).id = 14;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19,22]; % ???

unit{tid}(end+1).id = 16;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,22]; % 

unit{tid}(end+1).id = 18;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,16]; % 

unit{tid}(end+1).id = 33;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,16,19]; % 

unit{tid}(end+1).id = 38; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % 

unit{tid}(end+1).id = 44;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % Nope

unit{tid}(end+1).id = 45;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % Nope

unit{tid}(end+1).id = 66; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,16,17,19]; %

unit{tid}(end+1).id = 67; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,2,3]; %

unit{tid}(end+1).id = 74; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,2,3]; % 


% <<< PYR Ed10-20140816 <<< ---------
% >>>   PYR Ed10-20140817 >>> ---------

tid = find_trial_index(sessionlist,'Ed10-20140817.cof.gnd');

unit{tid}(1).id = 7;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,16];

unit{tid}(end+1).id = 9;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22];

unit{tid}(end+1).id = 10;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19,22,24];

unit{tid}(end+1).id = 17;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22,24];

unit{tid}(end+1).id = 33;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,16,19,22,23];

unit{tid}(end+1).id = 38;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19];

unit{tid}(end+1).id = 42;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,17];

unit{tid}(end+1).id = 57;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,9,10,11,13,14,16];

unit{tid}(end+1).id = 58;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,6,7,10];

unit{tid}(end+1).id = 60;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [8,10,11,13,14,16]; % NOPE

unit{tid}(end+1).id = 64;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,3]; % spatial cluster split with 66

unit{tid}(end+1).id = 70;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22];


% <<< PYR Ed10-20140817 <<< ---------
% <<< PYR Ed10 <<< ------------------
% >>> PYR jg05 >>> ------------------
% >>>   PYR jg05-20120309 >>> ---------

tid = find_trial_index(sessionlist,'jg05-20120309.cof.all');

unit{tid}(1).id = 13;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,16,19];

unit{tid}(end+1).id = 16;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,16,23];

unit{tid}(end+1).id = 20;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19];

unit{tid}(end+1).id = 34;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13];

unit{tid}(end+1).id = 53;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,14];

unit{tid}(end+1).id = 55;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [8,13,14];

unit{tid}(end+1).id = 59;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,11,13,16];

unit{tid}(end+1).id = 61;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [8,13,16];

unit{tid}(end+1).id = 62;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,19,21];

unit{tid}(end+1).id = 69;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,19];

unit{tid}(end+1).id = 70;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,21,22];

unit{tid}(end+1).id = 78;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,16,22];

unit{tid}(end+1).id = 80;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,16];

unit{tid}(end+1).id = 81;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19,22];


% <<< PYR jg05-20120309 <<< ---------
% >>>   PYR jg05-20120310 >>> ---------

tid = find_trial_index(sessionlist,'jg05-20120310.cof.all');
unit{tid}(1).id = 11;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,14,16]; % !!!
%  CHOOSE 
unit{tid}(end+1).id = 11;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19];

unit{tid}(end+1).id = 13;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,16,19];

unit{tid}(end+1).id = 18;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,15,19];

unit{tid}(end+1).id = 19;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,22];

unit{tid}(end+1).id = 20;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [2,4,8];

unit{tid}(end+1).id = 23;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,20,22];  

unit{tid}(end+1).id = 24;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16];

unit{tid}(end+1).id = 25;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,19];

unit{tid}(end+1).id = 33;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,12,13];

unit{tid}(end+1).id = 34;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,17];

unit{tid}(end+1).id = 35;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,19];

unit{tid}(end+1).id = 42;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,5,7,10];

unit{tid}(end+1).id = 46;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,22];

unit{tid}(end+1).id = 49;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,13];

unit{tid}(end+1).id = 50;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,15,18];

unit{tid}(end+1).id = 54;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,5,22];

unit{tid}(end+1).id = 56;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,16];
% CHOOSE
unit{tid}(end+1).id = 60;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,19,22];
% CHOOSE
unit{tid}(end+1).id = 60;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,19,22];
% CHOOSE
unit{tid}(end+1).id = 60;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,14,16,19];

unit{tid}(end+1).id = 60;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,16];

unit{tid}(end+1).id = 66;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,13,16];

unit{tid}(end+1).id = 78;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19];



% <<< jg05-20120310 <<< ---------
% >>>   PYR jg05-20120312 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120312.cof.all');

unit{tid}(1).id = 11;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22];

unit{tid}(end+1).id = 11;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,22];

unit{tid}(end+1).id = 12;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22];

unit{tid}(end+1).id = 18;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,19];

unit{tid}(end+1).id = 19;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,16];

unit{tid}(end+1).id = 23;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,18,19,21];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,11,14,16,21];

unit{tid}(end+1).id = 33;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16];

unit{tid}(end+1).id = 33;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,12,13,16,21];

unit{tid}(end+1).id = 57;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7];

unit{tid}(end+1).id = 57;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,15,23];

unit{tid}(end+1).id = 59;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,11,13];

unit{tid}(end+1).id = 61;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = 7;

unit{tid}(end+1).id = 61;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[1,4,7,13,15];

unit{tid}(end+1).id = 70;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[5,7];

unit{tid}(end+1).id = 77;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[7,13];

unit{tid}(end+1).id = 77;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[2,7,9,10];

unit{tid}(end+1).id = 83;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,7,13];

unit{tid}(end+1).id = 101;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[16,22];

unit{tid}(end+1).id = 101;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,7,16,22];

unit{tid}(end+1).id = 102;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,7,16,22];
% CHOOSE
unit{tid}(end+1).id = 105;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[13,16,19,22];
% CHOOSE
unit{tid}(end+1).id = 105;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[9,19,22];

unit{tid}(end+1).id = 106;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[19,22,23];
% CHOOSE
unit{tid}(end+1).id = 106;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[1,4,11,13,14,16,17,19,20,22,23];
% CHOOSE
unit{tid}(end+1).id = 106;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[1,4,7,16,17,19,22,23];

unit{tid}(end+1).id = 111;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,6];

unit{tid}(end+1).id = 114;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,10];
%CHOOSE
unit{tid}(end+1).id = 122;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[1,4,7,9,13];
%CHOOSE
unit{tid}(end+1).id = 122;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,8,13,19,21,24];
%CHOOSE
unit{tid}(end+1).id = 128;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[13,19,22,23];
%CHOOSE
unit{tid}(end+1).id = 128;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[1,9,13,19,21,22];
%CHOOSE
unit{tid}(end+1).id = 136;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,10,13];
%CHOOSE
unit{tid}(end+1).id = 136;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,11,13,17];

unit{tid}(end+1).id = 137;
unit{tid}(end).tid = tid;
unit{tid}(end).fset =[4,6];


% <<< jg05-20120312 <<< ---------
% >>>   PYR jg05-20120315 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120315.cof.all');

unit{tid}(1).id = 6;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,13,15,18,19,21];

unit{tid}(end+1).id = 22;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,12,13,15,16,18,19];

unit{tid}(end+1).id = 25;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,16,19];

unit{tid}(end+1).id = 27;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,11,15,22];

unit{tid}(end+1).id = 37;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [8,10];

unit{tid}(end+1).id = 61;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,17,19,20];



% <<< jg05-20120315 <<< ---------
% >>>   PYR jg05-20120316 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120316.cof.all');

unit{tid}(1).id = 13;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,13,15];

unit{tid}(end+1).id = 19;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,14,15,16,18];

unit{tid}(end+1).id = 30;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10,13,16];

unit{tid}(end+1).id = 41;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,12,16];

unit{tid}(end+1).id = 48;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,9,10,13,16];

unit{tid}(end+1).id = 50;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16];

unit{tid}(end+1).id = 51;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,7,13];

unit{tid}(end+1).id = 58;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [6,10,11,16,19];

unit{tid}(end+1).id = 65;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,23];


% <<< PYR jg05-20120316 <<< ---------
% >>>   PYR jg05-20120317 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120317.cof.all');

unit{tid}(1).id = 10;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,8];

unit{tid}(end+1).id = 18;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,9,10,12];

unit{tid}(end+1).id = 31;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10];

unit{tid}(end+1).id = 36;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19,23];

unit{tid}(end+1).id = 40;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,17,19,20];

unit{tid}(end+1).id = 50;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,9,12,16];

unit{tid}(end+1).id = 54;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10];

unit{tid}(end+1).id = 63;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13];

unit{tid}(end+1).id = 66;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4];

unit{tid}(end+1).id = 69;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,13];

unit{tid}(end+1).id = 72;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,10,13];



% <<< PYR jg05-20120317 <<< ---------
% >>>   PYR jg05-20120323 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120323.cof.all');

%CHOOSE
unit{tid}(1).id = 17; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,19];

%CHOOSE
unit{tid}(1).id = 18; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,20,22];
%CHOOSE
unit{tid}(end+1).id = 18; % STAR
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,16,19];

unit{tid}(end+1).id = 22;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,13,16];

unit{tid}(end+1).id = 26;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,10];

unit{tid}(end+1).id = 48;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,10];

unit{tid}(end+1).id = 54;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,10,11];


% <<< PYR jg05-20120323 <<< ---------
% >>>   PYR jg05-20120324 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120324.cof.all');


unit{tid}(1).id = 10;
unit{tid}(1).tid = tid;
unit{tid}(1).fset = [19,20,22];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,9,10];


% <<< PYR jg05-20120324 <<< ---------
% >>>   E PYR jg05-20120325 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120325.cof.all');


% $$$ unit{tid}(1).id = 10;
% $$$ unit{tid}(1).tid = tid;
% $$$ unit{tid}(1).fset = [19,20,22];
% $$$ 
% $$$ unit{tid}(end+1).id = 29;
% $$$ unit{tid}(end).tid = tid;
% $$$ unit{tid}(end).fset = [4,7,9,10];


% <<< PYR jg05-20120325 <<< ---------
% >>>   E PYR jg05-20120326 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120326.cof.all');


% $$$ unit{tid}(1).id = 10;
% $$$ unit{tid}(1).tid = tid;
% $$$ unit{tid}(1).fset = [19,20,22];
% $$$ 
% $$$ unit{tid}(end+1).id = 29;
% $$$ unit{tid}(end).tid = tid;
% $$$ unit{tid}(end).fset = [4,7,9,10];


% <<< PYR jg05-20120326 <<< ---------
% >>>   PYR jg05-20120329 >>> ---------
tid = find_trial_index(sessionlist,'jg05-20120329.cof.all');


unit{tid}(1).id = 4;
unit{tid}(1).tid = tid;
unit{tid}(1).fset = [4,5,7,9,13,15];

unit{tid}(end+1).id = 20;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [ 4,5,7,9,10,12,13,15,16,19];

unit{tid}(end+1).id = 29;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [ 3, 9, 11];

unit{tid}(end+1).id = 62;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];

unit{tid}(end+1).id = 67;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];

unit{tid}(end+1).id = 71;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [2,4,5,7,9,10,11,12,13,15,16,18,19,21,22,24];

unit{tid}(end+1).id = 76;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [12,18];

unit{tid}(end+1).id = 82;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];

unit{tid}(end+1).id = 83;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,8,9,10];

unit{tid}(end+1).id = 84;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [8,10,14,16,18];

unit{tid}(end+1).id = 97;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];


% <<< PYR jg05-20120329 <<< ---------
% <<< PYR jg05 <<< ------------------
tid = otid;
% <<< PYRAMIDAL CELLS <<< -----------------------------------------------------


avar_ccg_inter_swf = {};
avar_ratemap  = {};

for tid = find(~cellfun(@isempty,unit))
    try

% trial index


tid = find_trial_index(sessionlist,'Ed10-20140816.cof.all');
        
% $$$ tid = find_trial_index(sessionlist,'jg05-20120309.cof.all');
% $$$ tid = find_trial_index(sessionlist,'jg05-20120310.cof.all');
% $$$ tid = find_trial_index(sessionlist,'jg05-20120312.cof.all');
% >>> LOAD Trial >>> ----------------------------------------------------------
Trial = MTATrial.validate(sessionlist(tid));
Trial.par = MTAPar(Trial); % TODO : load in class
Trial.load('par'); % TODO : load in class
Trial.load('nq'); % TODO : load in class
Trial.spk.parent = Trial; % TODO : load in class
par = Trial.par;
% <<< LOAD Trial <<< ----------------------------------------------------------
% >>> LOAD Ephys and Bhv Vars >>> ---------------------------------------------
% >>> VAR composite : stc - Behavioral states >>>
stc = Trial.load('stc','msnn_ppsvd_raux');
% <<< VAR composite : stc - Behavioral states <<<
% >>> VAR ephys : spk - Neuron spikes >>>
spk = Trial.load('spk', ...
                 samplerate, ...
                 '', ...
                 arrayfun(@(u) u.id, unit{tid}), ...
                 '', ...
                 true, ...
                 true);
% <<< VAR : spk - Neuron spikes <<< 
% >>> VAR ephys : Normalized spike wave form features >>>
zfet = spk.fet;
for c = unique(spk.clu)'
    zfet(spk.clu==c,:) = p_zscore(zfet(spk.clu==c,:));
end
% <<< VAR ephys : Normalized spike wave form features <<<
% >>> VAR ephys : phz - Theta Phase >>>
phz = load_theta_phase( Trial, samplerate, Trial.meta.channelGroup.theta);
% <<< VAR ephys : phz - Theta Phase <<<
% >>> VAR composite : pft - Place Fields >>>
pft = pfs_2d_theta(Trial, [], 'theta-groom-sit-rear', true, [], true);
% <<< VAR composite: pft - Place Fields <<<
% >>> VAR kinematic : xyz, hba, hvfl, fvxy, roll >>>
xyz = preproc_xyz(Trial,'trb',samplerate);
hba = fet_hba(Trial,samplerate);
hvfl = fet_hvfl(Trial,samplerate);
fvxy = vel(filter(copy(xyz),'ButFilter',4,2.5,'low'),'hcom',[1,2]);
roll = fet_roll(Trial,samplerate); % -:right, +:left
% <<< VAR kinematic : xyz, hba, hvfl, fvxy, roll <<<
% >>> VAR kinematic : har - Head Body angle/roll feature >>>
theta = pi/3; %Ed10
theta = pi/2.8; %jg05
Vr = [cos(theta),-sin(theta); sin(theta),cos(theta)];
har = copy(hba);
har.data = [hba.data, roll.data];
har.data = har.data * Vr;
har.data = har(:,1);
% <<< hba roll feature <<<
% <<< LOAD Ephys and Bhv Vars <<< ---------------------------------------------

for u = 1:numel(unit{tid})
try
% >>> LOAD unit data >>> ------------------------------------------------------
uid = unit{tid}(u).id
[field_mrate, field_center] = pft.maxRate(uid);
ego = fet_ego( Trial, xyz, {'head_right','head_left'}, field_center, 0);
[res, fet, swf] = spk(unit{tid}(u).id, [stc{bhv_state, samplerate}]);
cind = spk.clu==uid & ismember(spk.res,res);
mswf=sq(mean(swf));
% <<< LOAD unit data <<< ------------------------------------------------------
% >>> SETUP figure >>> --------------------------------------------------------

figDir = create_directory('/storage/gravio/figures/analysis/ego-swf-har');
trlDir = create_directory(fullfile(figDir, Trial.filebase));
% phz, swf, hba
hfig = setup_figure(figure(23849843), ...
                    'A4', 'portrait', ...
                    'centimeters', ...
                    1.8, 1.8, 0.05, 0.15);

% <<< SETUP figure <<< --------------------------------------------------------
% >>> COMPUTE EGO postion >>> -------------------------------------------------

[field_mrate, field_center] = pft.maxRate(uid);
unit{tid}(u).mrate = field_mrate;
unit{tid}(u).mpos = field_center;
ego = fet_ego( Trial, xyz, [], field_center, 0);
fd = sqrt(sum(ego(res,:).^2,2));


% <<< COMPUTE EGO postion <<< -------------------------------------------------
% >>> COMPUTE allo occupancy map >>> ------------------------------------------
[ocnt,xb,yb] = hist2(sq(xyz([Trial.stc{bhv_state}],'hcom',[1,2])),...
                     linspace(-500,500,xbins),...
                     linspace(-500,500,ybins));
xc = mean([xb(1:end-1),xb(2:end)]');
yc = mean([xb(1:end-1),xb(2:end)]');
[W, H] = meshgrid( xb, yb);
mazeMaskAllo = sqrt(W.^2 + H.^2) > 480;
% <<< COMPUTE allo occupancy map <<< ------------------------------------------7
% >>> COMPUTE ego occupancy map >>> -------------------------------------------

[ecnt,fb,lb] = hist2(ego([Trial.stc{bhv_state}],:),...
                     linspace(-500,500,xbins),...
                     linspace(-500,500,ybins));
fc = mean([fb(1:end-1),fb(2:end)]');
lc = mean([lb(1:end-1),lb(2:end)]');
[W, H] = meshgrid( lb, lb);
mazeMaskEgo = sqrt(W.^2 + H.^2) > 350;

hba_full = hba([Trial.stc{bhv_state}],:);
hvl_full = hvfl([Trial.stc{bhv_state}],2);
ego_full = ego([Trial.stc{bhv_state}],:);
har_full = har([Trial.stc{bhv_state}],1);

ecnt_hba = {};
ecnt_hvl = {};
ecnt_har = {};
for hbaI = 1:hbaN
    ego_hvl = ego_full(within_ranges(hvl_full, bins.hvl.edges( hbaI+[0,1])),:);
    ego_hba = ego_full(within_ranges(hba_full, bins.hba.edges( hbaI+[0,1])),:);
    ego_har = ego_full(within_ranges(har_full, bins.har.edges( hbaI+[0,1])),:);
    [ecnt_hba{hbaI}] = hist2(ego_hba, linspace(-500,500,xbins), linspace(-500,500,ybins));
    [ecnt_hvl{hbaI}] = hist2(ego_hba, linspace(-500,500,xbins), linspace(-500,500,ybins));
    [ecnt_har{hbaI}] = hist2(ego_har, linspace(-500,500,xbins), linspace(-500,500,ybins));
end

% <<< COMPUTE ego occupancy map <<< -------------------------------------------
% >>> COMPUTE ffet >>> --------------------------------------------------------

sfet = p_zscore( median_filter(fet(:,unit{tid}(u).fset),250));
%rfet = p_zscore( fet(:, swf_fet));
fd = sqrt(sum(ego(res,:).^2,2));
fdist = sqrt(sum(ego(:,:).^2,2));

xswf=sq(swf(:,Trial.nq.maxAmpChan(uid),:))';
mispk_ind=[];
mispk_max=[];
mispk_half=[];
mispk_dhalf=[];
for s=1:size(xswf,2)
    ixswf=interp1(1:size(xswf,1),-xswf(:,s),linspace(1,size(xswf,1),size(xswf,1)*10));
    [mispk_max(s),mispk_ind(s)]=max(ixswf(100:end-100));
    mispk_ind(s) = mispk_ind(s)+100;
    mfh = find(ixswf(mispk_ind(s):end)<mispk_max(s)/2,1,'first');
    if ~isempty(mfh)
        mispk_half(s) = mfh;
    else
        mispk_half(s) = 0;
    end
    mdh = find(ixswf([mispk_ind(s)+mispk_half(s)]:end)<ixswf(mispk_ind(s)+mispk_half(s))/2,1,'first');
    if ~isempty(mdh)
        mispk_dhalf(s) = mdh;
    else
        mispk_dhalf(s) = 0;
    end
end

sfet = p_zscore( median_filter(fet(:,unit{tid}(u).fset),250));
tfet = p_zscore( ...
    [ ...
      (mispk_half+randn([1,size(xswf,2)])/5)',...
      (mispk_dhalf+randn([1,size(xswf,2)])/5)',...
      (mispk_ind+randn([1,size(xswf,2)])*2)' - mean(mispk_ind)' ...
    ] ...
);

[S,U,V] = svd([sfet,tfet],0);
ffet = [sfet,tfet] * V;
[~,pscrI] = max(abs(p_corr(fd(fd<300),ffet(fd<300,:))));

ffet = ffet(:,pscrI);


% <<< COMPUTE ffet <<< --------------------------------------------------------
% >>> Unit Information >>> ----------------------------------------------------
    FigInfo = uicontrol('Parent',hfig,...
                        'Style','text',...
                        'String',{['Unit: ',num2str(uid)],...
                        Trial.filebase,...
                        ['stcMode: ',Trial.stc.mode],...
                        ['eDist:   ',num2str(Trial.nq.eDist(uid))],...
                        ['Refrac:  ',num2str(log10(Trial.nq.Refrac(uid)))],...
                        ['SNR:     ',num2str(Trial.nq.SNR(uid))],...
                        ['AmpSym:  ',num2str(Trial.nq.AmpSym(uid))],...
                        ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(uid))]...
                   },...
                        'Units','centimeters',...
                        'Position',[2, 23.25,5,4]);
% <<< Unit Information <<< ----------------------------------------------------
% >>> ACCG >>> ----------------------------------------------------------------

tpnts = [res];
gpnts = [ones(size(res))];
grps =  [1];
[mccg,tbins,pairs] = CCG(tpnts, gpnts,                                      ...
                         1, 35,                                             ...
                         samplerate, grps,                                  ...
                         'hz2'                                              ...
);
sfob; % because emacs is wierd

% <<< ACCG <<< ----------------------------------------------------------------
% >>> PLOT ACCG >>> -----------------------------------------------------------
sax = setup_axes(                                                   ...
    hfig,                                                           ...
    3, 0,                                                           ...
    1, 0,                                                           ...
    0, 0,                                                           ...
    1, 2);
bar(tbins,mccg(:,1));
ylim(sax, [0,max(mccg(:))]);
grid(sax, 'on');
% <<< PLOT ACCG <<< -----------------------------------------------------------
% >>> COMPUTE ALLO ratemap >>> ------------------------------------------------

scnt = hist2(sq(xyz(res,'hcom',[1,2])), ...
             linspace(-500,500,xbins),...
             linspace(-500,500,ybins));
rmap = (p_imgaussfilt( scnt, 1.2) ./ p_imgaussfilt( ocnt, 1.2));
rmap = rmap(:);
rmap(mazeMaskAllo(:)) = nan;
rmap(ocnt(:)<(samplerate.*0.02)) = nan;
rmap = reshape(rmap,[numel(xb),numel(yb)]);

% <<< COMPUTE ALLO ratemap <<< ------------------------------------------------
% >>> IMAGESC ALLO ratemap >>> ------------------------------------------------

sax = setup_axes(hfig,4,-1,1,0);
imagescnan({xb, yb, rmap'.*samplerate },...
           'colorbarIsRequired', false, ...
           'colorMap', @jet);
text(275,400, num2str(max(rmap(:).*samplerate),'%2.0f'), 'Color', 'w');
sax.XTickLabel = {};
sax.YTickLabel = {};
axis('xy');

% <<< IMAGESC ALLO ratemap <<< ------------------------------------------------
% >>> COMPUTE ego ratemap >>> -------------------------------------------------

scnt = hist2(ego(res,:), ...
             linspace(-500,500,xbins),...
             linspace(-500,500,ybins));
rmap = p_imgaussfilt( scnt, 1.2) ./ p_imgaussfilt( ecnt, 1.2);
rmap = rmap(:);
rmap(mazeMaskEgo(:)) = nan;
rmap(ocnt(:)<(samplerate.*0.02)) = nan;
rmap = reshape(rmap,[numel(xb),numel(yb)]);

% <<< COMPUTE ego ratemap <<<  ------------------------------------------------
% >>> IMAGESC ego ratemap >>> -------------------------------------------------

sax = setup_axes(hfig,4,-1,2,0);
imagescnan({xb, yb, fliplr(rot90(rmap'.*samplerate,-1)) },...
           'colorbarIsRequired', false,  ...
           'colorMap',           @jet,  ...
           'colorMapFlipFlag',   true   ...
);
xlim([-350, 350]);
ylim([-350, 350]);
sax.XTickLabel = {};
sax.YTickLabel = {};
text(170,280, num2str(max(rmap(:).*samplerate),'%2.0f'),'Color','w');
Lines([],0,'w');
Lines(0,[],'w');
axis('xy');

% <<< IMAGESC ALLO ratemap <<< ------------------------------------------------
% >>> SCATTER (x=dist, y=phz, c=sfet) >>> -------------------------------------
sax = setup_axes(hfig,9,0,1,0,0,0,4.5,2);
scatter([fd;fd], [phz(res);phz(res)+2*pi], 3, [ffet;ffet], 'Filled');
colormap(sax, 'jet');
caxis   (sax, [-2,2]);
xlim    (sax, [0,350]);
ylim    (sax, [0,4*pi]);
Lines([],bins.phz.edges,'k');
Lines([],bins.phz.edges+2*pi,'k');

% <<< SCATTER (x=dist, y=phz, c=sfet) <<< -------------------------------------
% >>> CCG >>> -----------------------------------------------------------------
tpnts = [];
gpnts = [];
grps =  [];
dthresh = 600;
for s = 1:3
    pid =  within_ranges(ffet, bins.swf.edges([s,s+1])) ...
           & abs(ego(res,2)) < dthresh ...
           ;% = pid
    tres = res(pid);
    tpnts = cat(1, tpnts, tres);
    gpnts = cat(1, gpnts, s * ones([numel(tres),1]));
    grps = [grps,s];
end
bin_size = 2;
bin_halfwidth = 35;
[mccg,tbins,pairs] = CCG(tpnts,      gpnts,                             ...
                         bin_size,   bin_halfwidth,                     ...
                         samplerate, grps,                              ...
                         'hz2'                                          ...
                         );
avar_ccg_inter_swf{tid}(:,:,:,u) = mccg;
sfob;
% <<< CCG <<< -----------------------------------------------------------------
% >>> PLOT CCG >>> ------------------------------------------------------------
for s = 1:swfN
    for s2 = s:swfN
        sax = setup_axes(                                                   ...
            hfig,                                                           ...
            (s-1)+11, 1,                                                    ...
            s2-1+1, +3*(s2-1),                                              ...
            0, 0,                                                           ...
            1, 2.5);
        bar(tbins,mccg(:,s,s2));
        ylim(sax, [0,max(mccg(:))]);
        grid(sax, 'on');
    end
end
% <<< PLOT CCG <<< ------------------------------------------------------------
% >>> spike waveforms by phase >>> --------------------------------------------

clrs = summer(3);
for phzI = 1:phzN
    % >>> SETUP subplot >>>
    sax = setup_axes(                                                       ...
        hfig,                                                               ...
        (phzI-1)*3+3, -0.05*(phzI-1),                                       ...
        2, 3.2,                                                             ...
        0, 0,                                                               ...
        3, 1.5);
    hold(gca(), 'on');
    % <<< SETUP subplot <<<
    spk_waveforms = diff(swf(:, 1:nchan, :), 1, 3);
    sax.Color = [0.75,0.75,0.75];
    for s = 1:swfN
        ind = within_ranges(ffet, bins.swf.edges([s,s+1])) & ...
              within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]));
        plot(bsxfun(                                                        ...
                 @minus,                                                    ...
                 sq(mean(spk_waveforms(ind,:,5:end-5)))'.*3,1:2500:2500*nchan),           ...
             '-','Color',clrs(s,:)                                                  ...
         );
    end
    sax.YTickLabel = {};
    ylim(sax, [-2.6e4,6e3]);
    grid(sax,'on');
end

% <<< spike waveforms by phase <<< --------------------------------------------
cmax = 10;
% >>> EGO ratemaps ( hba X phz) >>> -------------------------------------------

% $$$ for swfI = 1:swfN
% $$$ for phzI = 1:phzN
% $$$     for hbaI = 1:hbaN
% $$$         dm = res(  ...
% $$$               within_ranges(hba(res),  bins.hba.edges([hbaI, hbaI+1]))  ...
% $$$             & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
% $$$             & within_ranges(sfet,      bins.swf.edges([swfI,swfI+1]))                ...
% $$$         );
% $$$         if ~isempty(dm)
% $$$             % >>> compute ego ratemap >>>
% $$$             scnt = hist2(sq(ego(dm,:)),                                     ...
% $$$                          linspace(-500,500,xbins),                          ...
% $$$                          linspace(-500,500,ybins));
% $$$             rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ecnt_hba{hbaI},1.2));
% $$$             rmap = rmap(:);
% $$$             rmap(mazeMaskEgo(:)) = nan;
% $$$             rmap(ecnt(:)<(samplerate.*0.02)) = nan;
% $$$             rmap = reshape(rmap,[numel(fb),numel(lb)]);
% $$$             % <<< compute ego ratemap <<<
% $$$             sax = setup_axes(hfig, ...
% $$$                              (phzN-phzI)*3+(swfN-swfI+1), -0.05*(swfI-1), ...
% $$$                              hbaI+5, -0.5, ...
% $$$                              0, 0, ...
% $$$                              1, 1);
% $$$             imagescnan({fb, lb, fliplr(rot90(rmap'.*samplerate .*3,-1))}, ...
% $$$                        [0,cmax], ...
% $$$                        'colorbarIsRequired', false, ...
% $$$                        'colorMap', @jet, ...
% $$$                        'colorMapFlipFlag',true);
% $$$             axis('xy');
% $$$             xlim([-300,300]);
% $$$             ylim([-300,300]);
% $$$             text(100,-250, num2str(max(rmap(:).*samplerate .*3),'%2.0f'),'Color','w');
% $$$             sax.XTickLabel = {};
% $$$             sax.YTickLabel = {};
% $$$             Lines([],0,'w');
% $$$             Lines(0,[],'w');
% $$$         end
% $$$     end
% $$$ end
% $$$ end

% <<< EGO ratemaps ( hba X phz) <<< -------------------------------------------
% >>> EGO ratemaps ( har X phz) >>> -------------------------------------------
for swfI = 1:swfN
    for phzI = 1:phzN
        for harI = 1:harN
            % >>> QUERY partitioned spikes >>>

            dm = res(  ...
                  within_ranges(har(res),  bins.har.edges([harI, harI+1]))  ...
                & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
                & within_ranges(ffet,      bins.swf.edges([swfI, swfI+1]))  ...
                & within_ranges(hba(res),  bins.hba.edges([1,end]))  ...
                );

            % <<< QUERY partition indicies <<<
            if ~isempty(dm)
                % >>> COMPUTE ego ratemap >>>
                scnt = hist2(sq(ego(dm,:)),                                 ...
                             linspace(-500,500,xbins),                      ...
                             linspace(-500,500,ybins));
                
                rmap =(p_imgaussfilt(scnt,1.4)                                ...
                       ./p_imgaussfilt(ecnt_har{harI},1.4));
                rmap = rmap(:);

                occMaskEgo = ecnt_har{harI}(:) < (samplerate.*0.02)/3 ;
                
                rmap( mazeMaskEgo(:) ) = nan;
                rmap( occMaskEgo(:)  ) = nan;
                
                rmap = reshape(rmap,[numel(fb),numel(lb)]);
                avar_ratemap{tid}(:,:,:,u,phzI,harI,swfI) = rmap;
                % <<< COMPUTE ego ratemap <<<
                % >>> SETUP subplot >>>
                sax = setup_axes(                                           ...
                    hfig,                                                   ...
                    (phzN-phzI)*3+(swfN-swfI+1), -0.05*(swfI-1),            ...
                    harI+5, -0.5,                                           ...
                    0, 0,                                                   ...
                    1, 1);
                % <<< SETUP subplot <<<
                % >>> IMAGESCNAN ratemap >>>
                imagescnan(                                                 ...
                    {                                                       ...
                        fb,                                                 ...
                        lb,                                                 ...
                        fliplr(rot90(rmap'.*samplerate .*3,-1))             ...
                    },                                                      ...
                    [ 0, cmax],                                             ...
                    'colorbarIsRequired', false,                            ...
                    'colorMap', @jet,                                       ...
                    'colorMapFlipFlag',true);
                % <<< IMAGESCNAN ratemap <<<                
                % >>> FORMAT subplot >>>
                axis('xy');
                xlim([-300,300]);
                ylim([-300,300]);
                text( 100, -250,                                            ...
                      num2str(max(rmap(:).*samplerate .*3),'%2.0f'),        ...
                      'Color','w');
                sax.XTickLabel = {};
                sax.YTickLabel = {};
                Lines([],0,'w');
                Lines(0,[],'w');
                % <<< FORMAT subplot <<<
            end
        end
    end
end


% <<< EGO ratemaps ( har X phz) <<< -------------------------------------------

% >>> PRINT figure >>> --------------------------------------------------------

pause(0.01);
FigName = ['ego-swf-har','_',Trial.filebase,'_unit-',num2str(uid,'%04.f')];
print(hfig, '-dpng', fullfile(trlDir, [FigName,'.png']));

% <<< PRINT figure <<< --------------------------------------------------------c
catch err
continue
end
end

catch err
continue
end
end
    
figure();
hold('on');
plot(hba(:),roll(:),'.k');
scatter(hba(res),roll(res),15,ego(res,2),'Filled');
colormap('jet');
