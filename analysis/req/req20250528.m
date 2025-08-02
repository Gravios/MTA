% req20250528
%

% >>> Project Vars >>> --------------------------------------------------------
sessionlist = get_session_list_v3('BehaviorPlaceCode');
samplerate = 250;% Hz

swfN = bins.swf.count;
hbaN = bins.hba.count;
hvlN = bins.hvl.count;
phzN = bins.phz.count;

sind = {};
sedges = bins.swf.edges;
scntrs = bins.swf.centers;
% <<< Project Vars <<< --------------------------------------------------------
% >>> Trial Ids >>>
tid = find_trial_index(sessionlist,'jg05-20120309.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120310.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120312.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120315.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120316.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120317.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120323.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120324.cof.all');
tid = find_trial_index(sessionlist,'jg05-20120329.cof.all');

tid = find_trial_index(sessionlist,'Ed10-20140814.cof.all');
tid = find_trial_index(sessionlist,'Ed10-20140815.cof.all');
tid = find_trial_index(sessionlist,'Ed10-20140816.cof.all');
tid = find_trial_index(sessionlist,'Ed10-20140817.cof.gnd');

tid = find_trial_index(sessionlist,'ER06-20130612.cof.all');
% <<< Trial Ids <<<
% >>> LOAD Trial >>>

% >>> LOAD Trial >>> ----------------------------------------------------------

Trial = MTATrial.validate(sessionlist(tid));
Trial.par = MTAPar(Trial); % TODO : load in class
Trial.load('par'); % TODO : load in class
Trial.load('nq'); % TODO : load in class
Trial.spk.parent = Trial; % TODO : load in class

% <<< LOAD Trial <<< ----------------------------------------------------------
% >>> LOAD Ephys and Bhv Vars >>> ---------------------------------------------

stc = Trial.load('stc','msnn_ppsvd_raux');
spk = Trial.load('spk', samplerate, '', [], '', true);
phz = load_theta_phase(Trial,samplerate,Trial.meta.channelGroup.theta);
pft = pfs_2d_theta(Trial);
xyz = preproc_xyz(Trial,'trb',samplerate);
hba = fet_hba(Trial,samplerate);
hvfl = fet_hvfl(Trial,samplerate);
fvxy = vel(filter(copy(xyz),'ButFilter',4,2.5,'low'),'hcom',[1,2]);
roll = fet_roll(Trial,samplerate); % -:right, +:left

% VAR : zfet 
% >>> Normalized spike wave form features >>> ---------------------------------

zfet = spk.fet;
for c = unique(spk.clu)'
    zfet(spk.clu==c,:) = zscore(zfet(spk.clu==c,:));
end

% <<< Normalized spike wave form features <<< ---------------------------------

% <<< LOAD Ephys and Bhv Vars <<< ---------------------------------------------
% <<< LOAD Trial <<<

% >>> PYRAMIDAL CELLS >>> -----------------------------------------------------
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

unit{tid}(end+1).id = 77;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,16,18];

unit{tid}(end+1).id = 95;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,22,23,24];

unit{tid}(end+1).id = 103;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,9,10,12,13,15];

unit{tid}(end+1).id = 105;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,15,18];

unit{tid}(end+1).id = 108;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,11];

unit{tid}(end+1).id = 111;
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

unit{tid}(end+1).id = 47;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [10,13,14,16,17];

unit{tid}(end+1).id = 48;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,3];

unit{tid}(end+1).id = 92;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,6,10];

unit{tid}(end+1).id = 94;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,11];

unit{tid}(end+1).id = 96;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [4,7,10];

unit{tid}(end+1).id = 86;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,19];

unit{tid}(end+1).id = 77;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,18];

unit{tid}(end+1).id = 51;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,4,7,10,18];


% <<< PYR Ed10-20140815 >>> ---------
% >>>   PYR Ed10-20140816 >>> ---------

tid = find_trial_index(sessionlist,'Ed10-20140816.cof.all');

unit{tid}(1).id = 7;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [16,18,21];

unit{tid}(end+1).id = 10;
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

unit{tid}(end+1).id = 38;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % 

unit{tid}(end+1).id = 44;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % Nope

unit{tid}(end+1).id = 45;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [13,14,16,17,19,22]; % Nope

unit{tid}(end+1).id = 66;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [7,10,16,17,19]; %

unit{tid}(end+1).id = 67;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [1,2,3]; %

unit{tid}(end+1).id = 74;
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
unit{tid}(end).fset = [10,13,16,19]

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
unit{tid}(1).id = 18;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [19,20,22];
%CHOOSE
unit{tid}(end+1).id = 18;
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
unit{tid}(end).fset = [ 3,9,11];

unit{tid}(end+1).id = 62;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];

unit{tid}(end+1).id = 67;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [];

unit{tid}(end+1).id = 71;
unit{tid}(end).tid = tid;
unit{tid}(end).fset = [2,4,5,7,9,10,11,12,13,15,16,18,19,21,22,24]

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
% <<< PYRAMIDAL CELLS <<< -----------------------------------------------------


% >>> GENERATE automatic spike feature >>> ------------------------------------
normalization = '';
uid = unit{tid}(end).id;
res = spk(uid);
cind = spk.clu==uid&ismember(spk.res,res);
pz = phz(res);
k = 1;
sfet_rho =[];
sfet_pval = [];
for chn = 1:3:22
    for nfet = 0:2
        sfet = spk.fet(cind,chn+nfet);
        [sfet_rho(chn+nfet), sfet_pval(chn+nfet)] = circ_corrcl(pz, sfet);
    end
    k = k + 1;
end
unit{tid}(end).fset = find(sfet_rho>0.15);
unit{tid}(end)
% <<< GENERATE automatic spike feature <<< ------------------------------------
% >>> PLOT SPK PCA FEATURES >>> -----------------------------------------------

normalization = '';
uid = unit{tid}(end).id;
sfet = unit{tid}(end).fset;
res = spk(uid);
cind = spk.clu==uid&ismember(spk.res,res);
pz = phz(res);
hfig = figure();
k = 1;
for chn = 1:3:22
    for nfet = 0:2
    subplot2(3,8,1+nfet,k);
        sfet = spk.fet(cind,chn+nfet);
        hist2([pz, sfet],12,12,normalization);
        Lines([],1500,'k');
    end
    k = k + 1;
end
hfig.Position(3) = hfig.Position(3)*2;



% <<< PLOT SPK PCA FEATURES <<< -----------------------------------------------


% >>> INTERNEURONS >>> --------------------------------------------------------
% >>> ED10 >>> ----------------------
% >>> INT Ed10-20140814 >>> ---------

int = spk(10, [stc{statesSel{sts}, samplerate}]);
int = spk(11, [stc{statesSel{sts}, samplerate}]);
int = spk(16, [stc{statesSel{sts}, samplerate}]);
%int = spk(29, [stc{statesSel{sts}, samplerate}]);
int = spk(30, [stc{statesSel{sts}, samplerate}]);
int = spk(50, [stc{statesSel{sts}, samplerate}]);
int = spk(63, [stc{statesSel{sts}, samplerate}]);
int = spk(85, [stc{statesSel{sts}, samplerate}]);
int = spk(94, [stc{statesSel{sts}, samplerate}]);
int = spk(99, [stc{statesSel{sts}, samplerate}]);

% <<< INT Ed10-20140814 <<< ---------
% >>> INT Ed10-20140815 >>> ---------
int =       ...
    [       ...
        20; ...
        25; ...
        26; ...
        29; ...
        38; ...
        60; ...
        63; ...
        67; ...
    ];

% <<< INT Ed10-20140815 <<< ---------
% >>> INT Ed10-20140816 >>> ---------
int = ...
    [ ...
        29; ... DSC_TRH_---
        30; ... DSC_---_---
        35; ... DSC_---_ASC
        41; ... ---_TRH_ASC
        47; ... ---_TRH_---
        48; ... DSC_TRH_---
        50; ... DSC_TRH_---
        52; ... ---_---_ASC
        53; ... DSC_---_---
    ];
% <<< INT Ed10-20140816 <<< ---------
% >>> INT Ed10-20140817 >>> ---------
int = ...
    [ ...
        11; ... DSC_TRH_ASC 
        31; ... DSC_TRH_---
        39; ... DSC_---_ASC
        47; ... DSC_TRH_---
        48; ... DSC_---_---
        51; ... ---_TRH_ASC
        52; ... ---_TRH_---
        53; ... DSC_TRH_---
        56; ... ---_---_ASC
        72; ... NA
        74; ...
        81; ...
    ];
% <<< INT Ed10-20140817 <<< ---------
% <<< ED10 <<< ----------------------
% >>> jg05 >>> ----------------------
% >>> INT jg05-20120309 >>> ---------
int = ...
    [ ...
        10; ... DSC_TRH_---
        11; ... DSC_TRH_---
        15; ... ---_---_ASC
        26; ... DSC_TRH_---
        27; ... ---_---_---
        33; ... ---_---_ASC
        64; ... ---_---_ASC
        66; ... ---_TRH_---
    ];
% <<< INT jg05-20120309 <<< ---------
% >>> INT jg05-20120310 >>> ---------

int = spk( 7, [stc{statesSel{sts}, samplerate}]);
int = spk( 27, [stc{statesSel{sts}, samplerate}]);
int = spk( 28, [stc{statesSel{sts}, samplerate}]);
int = spk( 43, [stc{statesSel{sts}, samplerate}]); %!!!
int = spk( 59, [stc{statesSel{sts}, samplerate}]); %!!!
int = spk( 71, [stc{statesSel{sts}, samplerate}]);

% <<< INT jg05-20120310 <<< ---------
% >>> INT jg05-20120312 >>> ---------
int = spk( 14, [stc{statesSel{sts}, samplerate}]);
int = spk( 41, [stc{statesSel{sts}, samplerate}]);
int = spk( 43, [stc{statesSel{sts}, samplerate}]);
int = spk( 48, [stc{statesSel{sts}, samplerate}]);
int = spk( 74, [stc{statesSel{sts}, samplerate}]);
int = spk( 75, [stc{statesSel{sts}, samplerate}]);
int = spk( 90, [stc{statesSel{sts}, samplerate}]);
int = spk(104, [stc{statesSel{sts}, samplerate}]);
int = spk(117, [stc{statesSel{sts}, samplerate}]);
% <<< INT jg05-20120312 <<< ---------
% >>> INT jg05-20120316 >>> ---------
int = ...
    [ ...
        27; ... DSC_TRH_---
        28; ... DSC_TRH_---
        29; ... DSC_TRH_---
        30; ... DSC_TRH_---                
    ];
% <<< INT jg05-20120316 <<< ---------
% >>> INT jg05-20120317 >>> ---------
int = ...
    [ ...
        44; ... DSC_TRH_---
        49; ... DSC_TRH_---
        52; ... DSC_TRH_---
    ];
% <<< INT jg05-20120317 <<< ---------
% >>> INT jg05-20120323 >>> ---------
int = ...
    [ ...
        2;  ... DSC_TRH_---
        4;  ... ---_TRH_ASC
        20; ... ---_---_---
        30; ... ---_---_ASC
        39; ... DSC_TRH_---
        40; ... ---_---_ASC
        51; ... ---_---_---        
    ];
% <<< INT jg05-20120323 <<< ---------
% >>> INT jg05-20120324 >>> ---------
int = ...
    [ ...
         1; ... ---_---_---        
         2; ... ---_---_---
         4; ... ---_---_---
         6; ... ---_---_---
         7; ... ---_---_---
         9; ... ---_---_---
        14; ... ---_---_---
        16; ... ---_---_---
        26; ... ---_---_---
        34; ... ---_---_---
        36; ... ---_---_---
        37; ... ---_---_---
        49; ... ---_---_---
        54; ... ---_---_---
        60; ... ---_---_---
        61; ... ---_---_---
    ];
% <<< INT jg05-20120324 <<< ---------
% >>> INT jg05-20120329 >>> ---------
int = ...
    [ ...
         2; ... ---_---_---
         3; ... ---_---_---
         7; ... ---_---_---
        16; ... ---_---_---
        31; ... ---_---_---
        33; ... ---_---_---
        36; ... ---_---_---
        37; ... ---_---_---
        38; ... ---_---_---
        47; ... ---_---_---
        48; ... ---_---_---
        69; ... ---_---_---
        74; ... ---_---_---
        90; ... ---_---_---
        96; ... ---_---_---
        99; ... ---_---_---
    ];
% <<< INT jg05-20120329 <<< ---------
% <<< jg05 <<< ----------------------
% <<< INTERNEURONS <<< --------------------------------------------------------

int1 = spk(33, [stc{statesSel{sts}, samplerate}]);
uint = 36;
int = spk(uint, [stc{statesSel{sts}, samplerate}]);

uid = unit{tid}(end).id;
% >>> FIGURE: ccg, PC X Interneuron >>> ---------------------------------------
figure
sts = 1;
fetI = unit{tid}(u).fset;

% >>> RECOMPUTE vars for theta RUN state >>> ----------------------------------

ufr = Trial.ufr.create(Trial, xyz, spk, uid, 0.1, 'boxcar',true);
[field_mrate, field_center] = pft.maxRate(uid);
ego = fet_ego( Trial, xyz, [], field_center);

% SELECT spikes by clusterId: uid
% SELECT periods by bhv/ephys state: statesSel
% RESAMPLE to new sampling rate: samplerate
res = spk(uid, [stc{ statesSel{sts}, samplerate}]);

% GENERATE Indexing for spike waveforms
cind = spk.clu==uid&ismember(spk.res,res);
% COMPUTE distance to field center
fd = sqrt(sum(ego(res,:).^2,2));
% >>> COMPUTE spike feature >>>

if numel(fetI)<2
    sfet = zfet(cind,fetI);
else
    sfet = zfet(cind,fetI);
    if sts == 1
        %[pcomps,pscrs] = pca(sfet);
        [S,U,V] = svd(sfet);
        sc = sfet * V;
        %[~,pscrI] = max(abs(corr(fd(fd<400),sc(fd<400,:))));
        sfet = sc(:,1);
    else
        sfet = MedianFilter(sfet * V(:,pscrI),50);
    end
end

% <<< COMPUTE spike feature <<<
isiR = log10(diff([1;res]./samplerate+0.0001));
isiL = log10(diff([res;1]./samplerate+0.0001));
pz  = phz(res);
fwd = ego(res,1);
lat = ego(res,2);    
hb  = hba(res);
fv  = fvxy(res);
ur  = ufr(res);

% $$$ figure();
% $$$ for s = 1:size(sc,2)
% $$$     subplot(1,size(sc,2),s);
% $$$     scatter(fd,pz,20,sc(:,s),'Filled');
% $$$ end
% $$$ colormap('jet')

% <<< RECOMPUTE vars for theta RUN state <<< ----------------------------------
% >>> INIT: ARGUMENTSS, ccg >>> -----------------------------------------------
accg = [];
bin_size = 2;
bin_halfwidth = 100;
grps = [1,2];
normalization = 'hz';
% <<< INIT: ARGUMENTSS, ccg <<< -----------------------------------------------
% >>> SP: IMAGESC, INT ratemap >>> --------------------------------------------
subplot2(4,4,1,1);
    plot(pft, uint, 1, 'colorbar');
% <<< SP: IMAGESC, INT ratemap <<< --------------------------------------------
% >>> SP: IMAGESC, PC ratemap >>> ---------------------------------------------

subplot2(4,4,1,2);
    plot(pft, uid, 1, 'colorbar');

% <<< SP: IMAGESC, PC ratemap <<< ---------------------------------------------
% >>> SP: BAR, auto ccg of interneuron >>> ------------------------------------
subplot2(4,4,1,3);
    histcirc(phz(int));
    Lines(rad2deg(bins.phz.edges),[],'k');
% <<< SP: BAR, auto ccg of interneuron <<< ------------------------------------
% >>> SP: BAR, ccg between PC and interneuron >>> -----------------------------

tpnts = ...
    [ ...
        res; ...
        int; ...
        int1;...
    ];
ipnts = ...
    [ ...
        1 * ones([numel(res),1]); ...
        2 * ones([numel(int),1]); ...
        3 * ones([numel(int1),1]); ...
    ];

[mccg,tbins] = CCG( ...
    tpnts, ipnts, bin_size, bin_halfwidth, samplerate, [1,2,3], 'hz');
sfob;
subplot2(4,4,2,4);
    bar(tbins,mccg(:,1,2));
subplot2(4,4,3,4);
    bar(tbins,mccg(:,2,3));
subplot2(4,4,4,4);
    bar(tbins,mccg(:,3,2));

% <<< SP: BAR, ccg between PC and interneuron <<< -----------------------------
for phzI = 1 : phzN
    % >>> SELECT spikes in theta phase bin >>> --------------------------------
    pind = within_ranges( phz(res), bins.phz.edges([phzI,phzI+1]));
    % <<< SELECT spikes in theta phase bin <<< --------------------------------
    for s = 1 : swfN
        % >>> CCG >>>

        sind{s} = pind & within_ranges(sfet,sedges([s,s+1]));
        tpnts = [res(sind{s}); int];
        ipnts = [ones([sum(sind{s}),1]); 2*ones([numel(int),1])];
        [mccg,tbins] = CCG(tpnts,      ipnts,                                ...
                           bin_size,   bin_halfwidth,                        ...
                           samplerate, grps,                                 ...
                           normalization                                     ...
                           );
        accg(s,:,phzI) = mccg(:,1,2);

        % <<< CCG <<<
    end
    % >>> SP: IMAGESC ccg >>> -------------------------------------------------

    subplot2(4,4,phzN - phzI+2,1:3);
        imagesc(tbins,scntrs,(accg(:,:,phzI)));
        %imagesc(tbins,scntrs,log10(accg(:,:,phzI)));
        Lines(0,[],'w');
        axis('xy');
        colorbar();
        caxis(prctile(nonzeros(accg),[0.1,99.9]));
        %caxis([log10(prctile(nonzeros(accg),[0.1,99.9]))]);
    % <<< SP: IMAGESC ccg <<< -------------------------------------------------
end
sfob;
% >>> SP: SCATTER of theta phase VS distance VS spike feature >>> -------------
subplot2(4,4,1,4);
    scatter(pz,fd,20,sfet,'Filled');
    colormap('jet');
    caxis([-3,3]);
    Lines(bins.phz.edges,[],'k');
% <<< SP: SCATTER of theta phase VS distance VS spike feature <<< -------------    

% <<< FIGURE: ccg, PC X Interneuron <<< ---------------------------------------
% >>> DIAGNOSTIC >>>
% $$$ figure();
% $$$ scatter(fd,pz,20,sc(:,1),'Filled');
% $$$ caxis([-2.5,2.5]);
% $$$ colormap('jet');
% <<< DIAGNOSTIC <<<
% >>> FIGURE: isi vs sfet >>> -------------------------------------------------

ufr = Trial.ufr.create(Trial,xyz,spk,uid,0.1,'boxcar',true);
ylm = [-4, 4];
figure();
xlm = [-3,4];
clm = [0,300];
states = {'theta-groom-sit-rear','sit-theta','groom'};

[field_mrate, field_center] = pft.maxRate(uid);
ego = fet_ego( Trial, xyz, [], field_center);

for sts = 1:numel(states)
    res = spk(uid, [stc{states{sts},samplerate}]);
    cind = spk.clu==uid & ismember(spk.res,res);
    if numel(fetI)<2
        sfet = zfet(cind,fetI);
    else
        sfet = zfet(cind,fetI);
        [S,U,V] = svd(sfet);
        sc = sfet * V;
        sfet = sc(:,1);
    end
    isiR = log10(diff([1;res]./samplerate+0.0001));
    isiL = log10(diff([res;1]./samplerate+0.0001));
    pz = phz(res);
    fd = sqrt(sum(ego(res,:).^2,2));
    fwd = ego(res,1);
    lat = ego(res,2);    
    hb  = hba(res);
    fv  = fvxy(res);
    ur  = ufr(res);
    subplot2(4,numel(states),1,sts);
    histcirc(pz)
    subplot2(4,numel(states),4,sts);
    pind = phz(res)<2& phz(res)>0.5;
    scatter(randn([sum(pind),1])/10 + isiL(pind), sfet(pind),10,fd(pind),'Filled');
    caxis(clm);
    ylim(ylm);
    xlim(xlm);
    Lines([],mean(ylm));        
    set(gca,'Color',[0.85,0.85,0.85]);    
    subplot2(4,numel(states),3,sts);
    pind = phz(res)<4 &phz(res)>2;
    scatter(randn([sum(pind),1])/10 + isiL(pind), sfet(pind),10,fd(pind),'Filled');
    caxis(clm);
    ylim(ylm);
    xlim(xlm);
    Lines([],mean(ylm));        
    set(gca,'Color',[0.85,0.85,0.85]);    
    subplot2(4,numel(states),2,sts);
    pind = phz(res)>4;
    scatter(randn([sum(pind),1])/10 + isiL(pind), sfet(pind),10,fd(pind),'Filled');
    caxis(clm);
    ylim(ylm);
    xlim(xlm);
    Lines([],mean(ylm));
    set(gca,'Color',[0.85,0.85,0.85]);        
end
colormap('jet');

% <<< PLOT isi vs sfet <<< ----------------------------------------------------
% >>> FIGURE: CCG, PC X Sfet >>> ----------------------------------------------

figure();
% >>> RECOMPUTE vars for theta RUN state >>> ----------------------------------

ufr = Trial.ufr.create(Trial,xyz,spk,uid,0.1,'boxcar',true);

% >>> COMP: EGO position w.r.t. PC >>> ----------------------------------------
[field_mrate, field_center] = pft.maxRate(uid);
ego = fet_ego( Trial, xyz, [], field_center);
% <<< COMP: EGO position w.r.t. PC <<< ----------------------------------------

res = spk(uid, [stc{statesSel{sts},samplerate}]);
cind = spk.clu==uid&ismember(spk.res,res);
fd = sqrt(sum(ego(res,:).^2,2));
if numel(fetI)<2
    sfet = zfet(cind,fetI);
else
    sfet = zfet(cind,fetI);
    if sts == 1
        %[pcomps,pscrs] = pca(sfet);
        [S,U,V] = svd(sfet);
        sc = sfet * V;
        [~,pscrI] = max(abs(corr(fd(fd<400),sc(fd<400,:))));
        sfet = sc(:,pscrI);
    else
        sfet = sfet * V(:,pscrI);
    end
end
isiR = log10(diff([1;res]./samplerate+0.0001));
isiL = log10(diff([res;1]./samplerate+0.0001));
pz = phz(res);
fwd = ego(res,1);
lat = ego(res,2);    
hb  = hba(res);
fv  = fvxy(res);
ur  = ufr(res);

% <<< RECOMPUTE vars for theta RUN state <<< ----------------------------------
% >>> INIT: ARGUMENTSS, ccg >>> -----------------------------------------------
bin_size = 1;
bin_halfwidth = 75;
grps = [1,2,3];
% <<< INIT: ARGUMENTSS, ccg <<< -----------------------------------------------
for phzI = 1 : phzN
    % >>> SELECT spikes in theta phase bin >>> --------------------------------
    pind = within_ranges( phz(res), bins.phz.edges([phzI,phzI+1]));
    accg = [];
    % <<< SELECT spikes in theta phase bin <<< --------------------------------    
    % >>> CCG >>>

    sind{s} = pind & within_ranges(sfet,sedges([1,2]));
    lres = res(sind{s});
    sind{s} = pind & within_ranges(sfet,sedges([2,3]));
    mres = res(sind{s});
    sind{s} = pind & within_ranges(sfet,sedges([3,4]));
    hres = res(sind{s});
    tpnts = [lres; mres; hres];
    ipnts = [1 * ones([numel(lres),1]);                                     ...
             2 * ones([numel(mres),1]);                                     ...
             3 * ones([numel(hres),1])];                                    ...
    [mccg,tbins,pairs] = CCG(tpnts,      ipnts,                             ...
                             bin_size,   bin_halfwidth,                     ...
                             samplerate, grps,                              ...
                             'hz2'                                          ...
                             );
    accg(:,phzI,1) = mccg(:,1,2);
    accg(:,phzI,2) = mccg(:,2,3);
    accg(:,phzI,3) = mccg(:,1,3);
    accg(:,phzI,4) = mccg(:,1,1);
    accg(:,phzI,5) = mccg(:,2,2);
    accg(:,phzI,6) = mccg(:,3,3);        

    % <<< CCG <<<
    sfob;
    % >>> SP: IMAGESC ccg >>> -------------------------------------------------
    subplot2(3,6,phzN - phzI+1,1);        bar(tbins,accg(:,phzI,1));
    subplot2(3,6,phzN - phzI+1,2);        bar(tbins,accg(:,phzI,2));
    subplot2(3,6,phzN - phzI+1,3);        bar(tbins,accg(:,phzI,3));
    subplot2(3,6,phzN - phzI+1,4);        bar(tbins,accg(:,phzI,4));
    subplot2(3,6,phzN - phzI+1,5);        bar(tbins,accg(:,phzI,5));
    subplot2(3,6,phzN - phzI+1,6);        bar(tbins,accg(:,phzI,6));
    % <<< SP: IMAGESC ccg <<< -------------------------------------------------
end
sfob;
% >>> Formating >>> -----------------------------------------------------------
ForAllSubplots( ...
    [ ...
        'set(gca,''YScale'',''linear'');' ...
        'ylim([0,1]);' ...
        'Lines(0,[],''r'');' ...
        'grid(''on'');' ...
    ] ...
    );
% <<< Formating <<< -----------------------------------------------------------

% <<< FIGURE: CCG, PC X Sfet <<< ----------------------------------------------
% >>> FIGURE: EGO, sfet, hba, speed >>> ---------------------------------------
zm = 1;
statesSel = {'theta-groom-sit-rear','lbhv&theta','hbhv&theta'};
ny = bins.phz.count+1;
nx = 6;
figure();
for sts = 1:3;
    % >>> RECOMPUTE vars for theta RUN state >>> ------------------------------

    ufr = Trial.ufr.create(Trial,xyz,spk,uid,0.1,'boxcar',true);
    [field_mrate, field_center] = pft.maxRate(uid);
    ego = fet_ego( Trial, xyz, [], field_center);
    
    res = spk(uid, [stc{statesSel{sts},samplerate}]);
    cind = spk.clu==uid&ismember(spk.res,res);
    % >>> Spike Waveforme Features >>> ----------------------------------------

    fd = sqrt(sum(ego(res,:).^2,2));
    if numel(fetI)<2
        sfet = zfet(cind,fetI);
    else
        sfet = zfet(cind,fetI);
        if sts == 1
            %[pcomps,pscrs] = pca(sfet);
            [S,U,V] = svd(sfet);
            sc = sfet * V;
            [~,pscrI] = max(abs(corr(fd(fd<400),sc(fd<400,:))));
            sfet = sc(:,1);
        else
            sfet = sfet * V(:,pscrI);
        end
    end

    % <<< Spike Waveforme Features <<< ----------------------------------------    

    isiR = log10(diff([1;res]./samplerate+0.0001));
    isiL = log10(diff([res;1]./samplerate+0.0001));

    pz  = phz(res);
    fwd = ego(res,1);
    lat = ego(res,2);    
    hb  = hba(res);
    fv  = fvxy(res);
    ur  = ufr(res);
    hl  = hvfl(res,2);
    hf  = hvfl(res,1);
    rl  = roll(res);
    %fet = hf;    clim = [-20,30];
    %fet = hb;    clim = [-1,1];
    fet = rl;    clim = [-0.5,0.5]; %clim = [-0.9,0.9];
    %fet = hl;    clim = [-20,20];

    % <<< RECOMPUTE vars for theta RUN state <<< ------------------------------
    % >>> PLOT ego vs caxis sfet >>> ------------------------------------------

    msize = 11;% px
    sax = gobjects([0,1]);
    for phzI = 1:bins.phz.count
        pind = within_ranges(phz(res),bins.phz.edges([phzI,phzI+1])) ...
               ;%& sfet > -1;
        switch sts
          case 1
            % >>> theta >>> 

        subplot2(ny, nx, 1, 1);
            histcirc(phz(res));
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 1);
            scatter( lat(pind), fwd(pind), msize, sfet(pind), 'Filled');
            caxis([-2.5,2.5]);
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 2);
            scatter( lat(pind), fwd(pind), msize, fet(pind), 'Filled');
            caxis(clim);
% $$$         sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 3);
% $$$             scatter(lat(pind),fwd(pind),msize,fv(pind),'Filled');
% $$$             caxis([0,60]);

        % <<< theta <<<
          case 2
            % >>> lbhv >>> 
        subplot2(ny, nx, 1, 3);
            histcirc(phz(res));
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 3);
            scatter( lat(pind), fwd(pind), msize, sfet(pind), 'Filled');
            caxis([-2.5,2.5]);
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 4);
            scatter( lat(pind), fwd(pind), msize, fet(pind), 'Filled');
            caxis(clim);
% $$$         sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 6);
% $$$             scatter(lat(pind),fwd(pind),msize,fv(pind),'Filled');
% $$$             caxis([0,60]);
        % <<< lbhv <<<
          case 3
            % >>> hbhv >>>
        subplot2(ny, nx, 1, 5);
            histcirc(phz(res));
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 5);
            scatter( lat(pind), fwd(pind), msize, sfet(pind), 'Filled');
            caxis([-2.5,2.5]);
        sax(end+1) = subplot2(ny, nx, bins.phz.count+2-phzI, 6);
            scatter( lat(pind), fwd(pind), msize, fet(pind), 'Filled');
            caxis(clim);
        % <<< hbhv <<<
        end
    end
    colormap('jet');

    % <<< PLOT ego vs caxis sfet <<< ------------------------------------------
    % >>> FORMATING >>> -------------------------------------------------------
switch zm
  case 1
af(@(s) daspect(s,[1,1,1]), sax);
af(@(s) ylim(s,[-250,350]), sax);
af(@(s) xlim(s,[-300,300]), sax);    
ForAllSubplots( ...
    ['Lines([],0);'...     
     'Lines(0,[]);']);
  case 2
af(@(s) daspect(s,[1,1,1]), sax);
af(@(s) ylim(s,[-600,600]), sax);
af(@(s) xlim(s,[-600,600]), sax);
ForAllSubplots( ...
    ['Lines([],0);'...     
     'Lines(0,[]);']);
end
% <<< FORMATING <<< -----------------------------------------------------------
end
% <<< FIGURE: EGO, sfet, hba, speed <<<<---------------------------------------


hang = copy(xyz);
hang.data = atan2(xyz(:,'nose',2)-xyz(:,'hcom',2),...
                  xyz(:,'nose',1)-xyz(:,'hcom',1));
mang = copy(xyz);
mang.data = atan2(xyz(:,'hcom',2),xyz(:,'hcom',1));
dang = copy(xyz);
dang.data = circ_dist(mang.data,hang.data);

cstate = 'lloc+lpause&theta';


er = roll([stc{cstate}], 1);
dr = dang([stc{cstate}], 1);
ha = hba([stc{cstate}], 1);
pz = phz([stc{cstate}], 1);

ar = [hba([stc{cstate}], 1),roll([stc{cstate}], 1)];
[Sr,Ur,Vr] = svd(ar,0);

theta = pi/2.8;
theta = pi/2.5;
theta = pi/2.3;
Vr = [cos(theta),-sin(theta); sin(theta),cos(theta)]

theta = pi/3;
Vr = [cos(theta),-sin(theta); sin(theta),cos(theta)]
har = copy(hba);
har.data = [hba.data, roll.data] 
har.data = har.data *Vr;
ar = har([stc{cstate}],:);
har.data = har(:,1);


%u =13
uid = unit{tid}(u).id
[field_mrate, field_center] = pft.maxRate(uid);
ego = fet_ego( Trial, xyz, [], field_center);
ex =  ego([stc{cstate}], 2);
ey =  ego([stc{cstate}], 1);
fd = sqrt(sum(ego([stc{cstate}],:).^2,2));
res = spk(uid, [stc{cstate,samplerate}]);
exs =  ego(res,2);
eys =  ego(res,1);
ers = roll(res,1);
drs = dang(res,1);
has =  hba(res,1);
ars = [ hba(res), roll(res)];
ars = ars*Vr;



figure();
% >>> SUBPLOT hba v roll >>>
subplot(121);hold('on');
plot(ha,er,'.','Color',[0.7,0.7,0.7]);
scatter(has,ers,20,exs,'Filled');
colormap('jet');
caxis([-200,200]);
grid('on');
% <<< SUBPLOT hba v roll <<<
% >>> SUBPLOT pca fet hba v roll >>>
subplot(122); hold('on');
plot( ar(:,1), ar(:,2), '.', 'Color',[0.7,0.7,0.7]);
scatter( ars(:,1), ars(:,2), 20, exs, 'Filled');
colormap('jet');
caxis([-200,200]);
daspect([1,1,1]);
grid('on');
% <<< SUBPLOT pca fet hba v roll <<<


figure(); 
for harI = 1 : harN
    subplot2( 2, harN, 1, harI);
        hold('on');
        eind = within_ranges(ar(:,1), bins.har.edges(harI+[0,1]));
        einds = within_ranges(ars(:,1), bins.har.edges(harI+[0,1]));
        plot( ex(eind), ey(eind), '.' , 'Color',[0.7,0.7,0.7]);
        plot(exs(einds),eys(einds),'.',...
             'Color',bins.hba.color(harI,:),...
             'MarkerSize',10);
        xlim([-250,250]);
        ylim([-250,250]);
end
 
for hbaI = 1 : hbaN
    subplot2( 2, harN, 2, hbaI);
        hold('on');
        eind  = within_ranges(ha,  bins.hba.edges(hbaI+[0,1]));
        einds = within_ranges(has, bins.hba.edges(hbaI+[0,1]));
        plot( ex(eind), ey(eind), '.' , 'Color',[0.7,0.7,0.7]);
        plot(exs(einds),eys(einds),'.',...
             'Color',bins.hba.color(hbaI,:),...
             'MarkerSize',10);
        xlim([-250,250]);
        ylim([-250,250]);
end




figure(); hold('on');
plot(ha,er,'.','Color',[0.7,0.7,0.7]);
scatter(has,ers,20,exs,'Filled');
xlim([-pi/2,pi/2]);
ylim([-pi/2,pi/2]);
colormap('jet');
caxis([-200,200]);


figure();
subplot(221); hold('on');
    ind = dr>0 & er<0;
    inds = drs>0 & ers<0; 
    scatter(ex(ind),ey(ind),11,er(ind),'Filled');
    plot(exs(inds),eys(inds),'.m','MarkerSize',15);
    xlim([-250,250]);ylim([-250,250]);
    colormap('jet');
    caxis([-0.5,0.5]);
subplot(222); hold('on');
    ind = dr>0 & er>0;
    inds = drs>0 & ers>0; 
    scatter(ex(ind),ey(ind),11,er(ind),'Filled');
    plot(exs(inds),eys(inds),'.m','MarkerSize',15);
    xlim([-250,250]);ylim([-250,250]);
    colormap('jet');
    caxis([-0.5,0.5]);
subplot(223); hold('on');
    ind = dr<0 & er<0;
    inds = drs<0 & ers<0; 
    scatter(ex(ind),ey(ind),11,er(ind),'Filled');
    plot(exs(inds),eys(inds),'.m','MarkerSize',15);
    xlim([-250,250]);ylim([-250,250]);
    colormap('jet');
    caxis([-0.5,0.5]);
subplot(224); hold('on');
    ind = dr<0 & er>0;
    inds = drs<0 & ers>0; 
    scatter(ex(ind),ey(ind),11,er(ind),'Filled');
    plot(exs(inds),eys(inds),'.m','MarkerSize',15);
    xlim([-250,250]);ylim([-250,250]);
    colormap('jet');
    caxis([-0.5,0.5]);


% FIGURE roll vs direction through field
figure(); hold('on');
plot(er,dr,'.','MarkerSize',10,'Color',[0.7].*[1,1,1])
scatter(ers,drs,10,has,'Filled');
colormap('jet');
caxis([-1,1]);

figure(); hold('on');
plot(ha,er,'.','MarkerSize',10,'Color',[0.7].*[1,1,1])
scatter(has,ers,10,drs,'Filled');
colormap('hsv');
caxis([-pi,pi]);

figure(); 
dst = 300;
subplot(2,2,1); ind = fd < dst & within_ranges(dr,[-pi,-pi/2]); plot(ha(ind),er(ind),'.','MarkerSize', 10); 
subplot(2,2,2); ind = fd < dst & within_ranges(dr,[-pi/2,0]); plot(ha(ind),er(ind),'.','MarkerSize', 10);
subplot(2,2,3); ind = fd < dst & within_ranges(dr,[0,pi/2]); plot(ha(ind),er(ind),'.','MarkerSize', 10);
subplot(2,2,4); ind = fd < dst & within_ranges(dr,[pi/2,pi]); plot(ha(ind),er(ind),'.','MarkerSize', 10);
ForAllSubplots('xlim([-3,3]); ylim([-2,2]);grid(''on'');')


%11,23,70,77,105,106
% >>> RECOMPUTE vars for theta RUN state >>> ----------------------------------

sts = 1;
uset = [11,60];
fset = { [10,13,14,16], [10,16]};
% $$$ fset = { [10,13,14,16], [7,8,16,19]};
uset = [33, 49];
fset = {[7,10,12,13],[4,10,13]};
uset = [33, 49];
fset = {[7,10,12,13],[4,10,13]};
% >>> Ed10-20140816 >>>

% $$$ uset = [10,66];
% $$$ fset = {[7,10,16,17,19], [5,7,10,11,16,21]};

% <<< Ed10-20140816 <<<

% >>> jg05-20120312 >>> 

% $$$ uset = [29,105];
% $$$ fset = {[10,11,14,16,21],[9,19,22]};
% $$$ uset = [23,106];
% $$$ fset = {[16,18,19,21],[19,22,23]};
% $$$ uset = [23,77];
% $$$ fset = {[16,18,19,21],[2,7,9,10]};
% $$$ uset = [23,12];
% $$$ fset = {[16,18,19,21],[19,22]};
% $$$ uset = [12,77];
% $$$ fset = {[19,22],[2,7,9,10]};
% $$$ uset = [12,106];
% $$$ fset = {[19,22],[19,22,23]};
% $$$ uset = [11,128];
% $$$ fset = {[16,19,22],[10,16,19,22]};

% <<< jg05-20120312 <<<
clear('unit');
unit.id = 0;
for u = 1:numel(uset)
    unit(u).id = uset(u);
    res = spk(unit(u).id, [stc{'lloc+pause&theta', samplerate}]);
    unit(u).fset = fset{u};
    unit(u).res = res;
    unit(u).phz = phz(res);
    unit(u).hba = hba(res);    
    % >>> COMP: EGO position w.r.t. PC >>> ----------------------------------------

    [field_mrate, field_center] = pft.maxRate(unit(u).id);
    unit(u).mrate = field_mrate;
    unit(u).mpos  = field_center;
    unit(u).ego   = fet_ego( Trial, xyz, [], field_center);

    % <<< COMP: EGO position w.r.t. PC <<< ----------------------------------------
    % >>> COMP: SPK feature >>> ---------------------------------------------------
    sfet = zfet( spk.clu==unit(u).id & ismember(spk.res,res), fset{u});
    fd = sqrt(sum(unit(u).ego(:,:).^2,2));
    [S,U,V] = svd(sfet);
    sc = sfet * V;
    [~,pscrI] = max(abs(corr(fd(fd<300),sc(fd<300,:))));
    unit(u).sfet = zscore(MedianFilter(sc(:,pscrI),50));
    % <<< COMP: SPK feature <<< ---------------------------------------------------    
end

%% DOUBLE CHECK PUROPOSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% $$$ % >>> COMP: EGO position w.r.t. Arena >>> -------------------------------------
% $$$ %[mxr,mxp] = pft.maxRate(unit(2).id);
% $$$ mxp = [400,400];
% $$$ %mxp = [150,500];
% $$$ unit(u).mrate = mxr;
% $$$ unit(u).mpos = mxp;
% $$$ ego = MTADfet.encapsulate(                ... 3. Put into MTA object
% $$$     Trial,                                ...
% $$$     multiprod(                            ... 2. Rotate trajs to head
% $$$         bsxfun(                           ... 1. Set MXP center as Origin
% $$$             @minus,                       ...
% $$$             unit(u).mpos,                 ...
% $$$             sq(xyz(:,'hcom',[1,2]))),     ...
% $$$         hvec(:,:,:),2,[2,3]),             ...
% $$$     samplerate,                           ...
% $$$     'egocentric_placefield',              ...
% $$$     'egopfs',                             ...
% $$$     'p'                                   ...
% $$$     );
% $$$ % <<< COMP: EGO position w.r.t. Arena <<< -------------------------------------

% <<< RECOMPUTE vars for theta RUN state <<< ----------------------------------

direction = {@lt,@gt};

dthresh = 300;
figure();
for u = 1:2    
    subplot2(3,3,u,1);
    hold('on');
    plot(pft,unit(u).id,1,'colorbar');
    title(['Unit: ',num2str(unit(u).id)]);
    plot(mxp(:,1),mxp(:,2),'*m');
end
for d = 1:2
    % >>> INIT: ARGUMENTSS, ccg >>> -----------------------------------------------
    bin_size = 2;
    bin_halfwidth = 25;
    % <<< INIT: ARGUMENTSS, ccg <<< -----------------------------------------------
    for phzI = 1 : phzN
        % >>> CCG >>>

        tpnts = [];
        gpnts = [];
        grps = [];
        for u = 1:numel(unit)
            if u == 1
                pid =   ... within_ranges(unit(u).phz,  bins.phz.edges([phzI,phzI+1]))    ...
                ...within_ranges(unit(u).sfet, sedges([1,4])) ...
                within_ranges(unit(u).hba, [-0.4,0.4]) ...
                    & abs(ego(unit(u).res,2)) < dthresh ...
                    & direction{d}(ego(unit(u).res,1), 0)  ...
                    ;% = pid
            elseif u == 2
                pid =  within_ranges(unit(u).phz,  bins.phz.edges([phzI,phzI+1]))    ...
                       & within_ranges(unit(u).hba, [-0.4,0.4]) ...
                       & within_ranges(unit(u).sfet, sedges([1,4])) ...
                       & abs(ego(unit(u).res,2)) < dthresh ...
                       & direction{d}(ego(unit(u).res,1), 0)  ...
                       ;% = pid
            end
            tres = unit(u).res(pid);
            tpnts = cat(1, tpnts, tres);
            gpnts = cat(1, gpnts, u * ones([numel(tres),1]));
            grps = [grps,u];
        end
        [mccg,tbins,pairs] = CCG(tpnts,      gpnts,                             ...
                                 bin_size,   bin_halfwidth,                     ...
                                 samplerate, grps,                              ...
                                 'count'                                          ...
                                 );

        % <<< CCG <<<
        % >>> SP: IMAGESC ccg >>> -------------------------------------------------
% $$$     subplot2(3,3,phzN - phzI+1,1);
% $$$         bar(tbins,mccg(:,1,2));
% $$$         Lines(0,[],'r');
% $$$     subplot2(3,3,phzN - phzI+1,2);
% $$$         bar(tbins,mccg(:,1,1));
% $$$         Lines(0,[],'r');
% $$$     subplot2(3,3,phzN - phzI+1,3);
% $$$         bar(tbins,mccg(:,2,2));
% $$$         Lines(0,[],'r');
        
        subplot2(3,3,phzN - phzI+1,d+1);
        bar(tbins,mccg(:,1,2));
        Lines(0,[],'r');
        title(['Dir: ',num2str(d),...
               ' Phase: ',num2str(phzI) ,...
               ' Unit: ',num2str(unit(1).id), ...
               ' -> ','Unit: ',num2str(unit(2).id)]);
        % <<< SP: IMAGESC ccg <<< -------------------------------------------------
    end
end



hfig = figure();
dthresh = 200;
% >>> INIT: ARGUMENTSS, ccg >>> -----------------------------------------------
bin_size = 1;
bin_halfwidth = 50;
direction = {@gt,@lt};
% <<< INIT: ARGUMENTSS, ccg <<< -----------------------------------------------
for s = 1:3
    for d = 1:2
        % >>> CCG >>>
        tpnts = [];  gpnts = [];  grps = [];
        for u = 1:numel(unit)
            if u == 1
                pid =    within_ranges(unit{tind}(u).sfet, sedges([1,4])) ...
                         & within_ranges(unit(u).hba, [-0.6,0.6]) ...  
                         & abs(ego(unit(u).res,2)) < dthresh ...
                         & direction{d}(ego(unit(u).res,1), 0)  ...
                         ;% = pid
            elseif u == 2
                pid =    within_ranges(unit(u).sfet, sedges([s,s+1])) ...
                         & within_ranges(unit(u).hba, [-0.6,0.6]) ...                       
                         & abs(ego(unit(u).res,2)) < dthresh ...
                         & direction{d}(ego(unit(u).res,1), 0)  ...
                         ;% = pid
            end
            tres = unit(u).res(pid);
            tpnts = cat(1, tpnts, tres);
            gpnts = cat(1, gpnts, u * ones([numel(tres),1]));
            grps = [grps,u];
        end

        [mccg,tbins,pairs] = CCG(tpnts,      gpnts,                             ...
                                 bin_size,   bin_halfwidth,                     ...
                                 samplerate, grps,                              ...
                                 'count'                                          ...
                                 );
        sfob;

        % <<< CCG <<<
        % >>> SP: IMAGESC ccg >>> ---------------------------------------------

        subplot2(6,3,d+2*(s-1),1);
        bar(tbins,mccg(:,1,2));
        Lines(0,[],'r');
        subplot2(6,3,d+2*(s-1),2);
        bar(tbins,mccg(:,1,1));
        Lines(0,[],'r');
        subplot2(6,3,d+2*(s-1),3);
        bar(tbins,mccg(:,2,2));
        Lines(0,[],'r');
        % <<< SP: IMAGESC ccg <<< ---------------------------------------------
    end
end
sfob;
hfig.Position = hfig.Position .* [0,0.25,3,2];



figure
for d = 1:2
tic
dm = unit(2).res(direction{d}(ego(unit(2).res,1), 0) ...
                 & unit(2).phz<2.5 ...
                 & within_ranges(unit(2).hba, [-0.6,0.6]) ...
                 & abs(ego(unit(2).res,2)) < dthresh ...
                 & within_ranges(unit(2).sfet, sedges([1,4])) ...                 
                 ) ...  
     - unit(1).res(direction{d}(ego(unit(1).res,1), 0) ...
                 & unit(1).phz>3.5 ...
                 & within_ranges(unit(1).hba, [-0.6,0.6]) ...
                 & abs(ego(unit(1).res,2)) < dthresh ...
                 & within_ranges(unit(1).sfet, sedges([1,2])) ...
                   )';
dm = dm(:);
dm = dm(dm<50&dm>-50);
toc
subplot(2,1,d);
hist(dm/500,linspace(-0.1, 0.1, 100));
title(num2str(d));
xlim([-0.1,0.1]);
end


% >>> Figure spk feature summary >>> ------------------------------------------

tid = 23;
u = 6;
uid = unit{tid}(u).id
xbins = 36;
ybins = 36;
bhv_state = 'lloc+pause&theta';
res = spk(unit{tid}(u).id, [stc{bhv_state, samplerate}]);

% >>> COMPUTE EGO postion >>> -------------------------------------------------

[unit{tid}(u).mrate, unit{tid}(u).mpos ] = pft.maxRate(unit{tid}(u).id);
ego = fet_ego( Trial, xyz, [], unit{tid}(u).mpos);

% <<< COMPUTE EGO postion <<< -------------------------------------------------
% >>> COMPUTE allo occupancy map >>> ------------------------------------------
[ocnt,xb,yb] = hist2(sq(xyz([Trial.stc{bhv_state}],'hcom',[1,2])),...
                     linspace(-500,500,xbins),...
                     linspace(-500,500,ybins));
xc = mean([xb(1:end-1),xb(2:end)]');
yc = mean([xb(1:end-1),xb(2:end)]');
[W, H] = meshgrid( xb, yb);
mazeMaskAllo = sqrt(W.^2 + H.^2) > 480;
% <<< COMPUTE allo occupancy map <<< ------------------------------------------
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

ecnt_hba = {};
ecnt_hvl = {};
for hbaI = 1:hbaN
    ego_hvl = ego_full(within_ranges(hvl_full, bins.hvl.edges([hbaI,hbaI+1])),:);
    ego_hba = ego_full(within_ranges(hba_full, bins.hba.edges([hbaI,hbaI+1])),:);
    [ecnt_hba{hbaI}] = hist2(ego_hba, linspace(-500,500,xbins), linspace(-500,500,ybins));
    [ecnt_hvl{hbaI}] = hist2(ego_hba, linspace(-500,500,xbins), linspace(-500,500,ybins));
end

% <<< COMPUTE ego occupancy map <<< -------------------------------------------
% >>> COMPUTE sfet >>> --------------------------------------------------------
sfet = zfet( spk.clu==unit{tid}(u).id & ismember(spk.res,res), unit{tid}(u).fset);
fd = sqrt(sum(ego(res,:).^2,2));
[S,U,V] = svd(sfet);
sc = sfet * V;
[~,pscrI] = max(abs(corr(fd(fd<300),sc(fd<300,:))));
sfet = zscore(MedianFilter(sc(:,pscrI),50));
% <<< COMPUTE sfet <<< --------------------------------------------------------


hfig = setup_figure(figure(23849839), ...
                    'A2', 'landscape', ...
                    'centimeters', ...
                    4.5, 4, 0.75, 0.75);


hfig = setup_figure(figure(23849842), ...
                    'A4', 'landscape', ...
                    'centimeters', ...
                    2.5, 2.25, 0.45, 0.45);

% >>> COMPUTE ALLO ratemap >>> ------------------------------------------------
scnt = hist2(sq(xyz(res,'hcom',[1,2])), ...
             linspace(-500,500,xbins),...
             linspace(-500,500,ybins));
rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ocnt,1.2));
rmap = rmap(:);
rmap(mazeMask(:)) = nan;
rmap(ocnt(:)<(samplerate.*0.02)) = nan;
rmap = reshape(rmap,[numel(xb),numel(yb)]);
% <<< COMPUTE ALLO ratemap <<< ------------------------------------------------
% >>> IMAGESC ALLO ratemap >>> ------------------------------------------------
sax = setup_axes(hfig,1,0,1,0);
imagescnan({xb,yb,rmap'.*samplerate }, 'colorbarIsRequired', true, ...
           'colorMap', @jet);
axis('xy');
% <<< IMAGESC ALLO ratemap <<< ------------------------------------------------
% >>> COMPUTE ego ratemap >>> -------------------------------------------------
% $$$ scnt = hist2(ego(res,:), ...
% $$$              linspace(-500,500,xbins),...
% $$$              linspace(-500,500,ybins));
% $$$ rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ecnt,1.2));
% $$$ rmap = rmap(:);
% $$$ rmap(mazeMask(:)) = nan;
% $$$ rmap(ocnt(:)<(samplerate.*0.02)) = nan;
% $$$ rmap = reshape(rmap,[numel(xb),numel(yb)]);
% <<< COMPUTE ALLO ratemap <<< ------------------------------------------------% >>> COMPUTE ego ratemap >>> -------------------------------------------------
% >>> SCATTER (x=dist, y=phz, c=sfet) >>> -------------------------------------
sax = setup_axes(hfig,3,0,1,0,0,0,2,2);
scatter(fd, phz(res), 6, sfet, 'Filled');
colormap(gca(), 'jet');
caxis   (gca(), [-2.5,2.5]);
xlim    (gca(), [0,350]);
% <<< SCATTER (x=dist, y=phz, c=sfet) <<< -------------------------------------
% >>> spike waveforms by phase >>> --------------------------------------------
clrs = 'bgr';
for phzI = 1:phzN
    sax = setup_axes(                                                       ...
        hfig,                                                               ...
        4, 0,                                                               ...
        1, (hfig.UserData.subplot.width/3+0.25)*(phzI-1)+(phzI-1),          ...
        0, 0,                                                               ...
        1, 0.6);
    hold(gca(), 'on');
    spk_waveforms = mspk.spk(  mspk.clu==unit{tid}(u).id                    ...
                               & ismember(mspk.res,res),                    ...
                              :,:);
    for s = 1:numel(scntrs)
        ind = within_ranges(sfet, sedges([s,s+1])) & ...
              within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]));
        plot(bsxfun(                                                        ...
                 @minus,                                                    ...
                 sq(mean(spk_waveforms(ind,:,:)))',1:2500:20000),           ...
             ['-',clrs(s)]                                                  ...
         );
    end
    if phzI > 1
        sax.YTickLabel = {};
    end
    ylim(sax, [-2e4,1e3]);
end
% <<< spike waveforms by phase <<< --------------------------------------------
% >>> spike waveforms by spk feature >>> --------------------------------------
clrs = cool(3);
for s = 1:numel(scntrs)
    sax = setup_axes(hfig, ...
                     5, 0, ...
                     1, (hfig.UserData.subplot.width/3+0.25)*(s-1)+s-1, ...
                     0, 0, ...
                     1, 0.6);
    hold(gca(),'on');
    spk_waveforms = mspk.spk(mspk.clu==unit{tid}(u).id & ismember(mspk.res,res),:,:) ;
    for phzI = 1:phzN
        ind = within_ranges(sfet, sedges([s,s+1])) & ...
              within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]));
        sum(ind)
        plot(bsxfun(@minus,sq(mean(spk_waveforms(ind,:,:)))',1:2500:20000),'-','Color',clrs(phzI,:));
    end
    if s > 1
        sax.YTickLabel = {};
    end
    ylim(sax,[-2e4,1e3]);
end
% <<< spike waveforms by spk feature <<< --------------------------------------
% >>> ALLO ratemaps ( swf X phz) >>> ------------------------------------------
for phzI = 1:phzN
    for s = 1:numel(scntrs)
        dm = res(  within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]))  ...
                 & within_ranges(sfet,      sedges([s,s+1]))                ...
        );
        if ~isempty(dm)
            % >>> compute ego ratemap >>>
            scnt = hist2(sq(xyz(dm,'hcom',[1,2])),                          ...
                         linspace(-500,500,xbins),                          ...
                         linspace(-500,500,ybins));
            rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ocnt,1.2));
            rmap = rmap(:);
            rmap(mazeMaskAllo(:)) = nan;
            rmap(ocnt(:)<(samplerate.*0.02)) = nan;
            rmap = reshape(rmap,[numel(xb),numel(yb)]);
            % <<< compute ego ratemap <<<            
            sax = setup_axes(hfig, (phzN-phzI+1), 0, s+1, 0, 0, 4, 1, 1);
            imagescnan({xb,yb,rmap'.*samplerate .*3},                       ...
                       'colorbarIsRequired', true,                          ...
                       'colorMap', @jet);
            axis('xy');
            sax.XTickLabel = {};
            sax.YTickLabel = {};
        end
    end
end
% <<< ALLO ratemaps ( swf X phz) <<< ------------------------------------------
% >>> EGO ratemaps ( swf X phz) >>> -------------------------------------------
for phzI = 1:phzN
    for s = 1:numel(scntrs)
        dm = res(  within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]))  ...
                 & within_ranges(sfet,      sedges([s,s+1]))                ...
        );
        if ~isempty(dm)
            % >>> compute ego ratemap >>>
            scnt = hist2(sq(ego(dm,:)),                                     ...
                         linspace(-500,500,xbins),                          ...
                         linspace(-500,500,ybins));
            rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ecnt,1.2));
            rmap = rmap(:);
            rmap(mazeMaskEgo(:)) = nan;
            rmap(ecnt(:)<(samplerate.*0.02)) = nan;
            rmap = reshape(rmap,[numel(fb),numel(lb)]);
            % <<< compute ego ratemap <<<
            sax = setup_axes(hfig, (phzN-phzI+1), 0, s+4, 0, 0, 4.25, 1, 1);
            imagescnan({fb, lb, fliplr(rot90(rmap'.*samplerate .*3',-1))}, ...
                       'colorbarIsRequired', true, ...
                       'colorMap', @jet, ...
                       'colorMapFlipFlag',false);
            axis('xy');
            xlim([-300,300]);
            ylim([-300,300]);
            sax.XTickLabel = {};
            sax.YTickLabel = {};
            Lines([],0,'w');
            Lines(0,[],'w');
        end
    end
end
% <<< EGO ratemaps ( swf X phz) <<< -------------------------------------------
s = 3;
% >>> EGO ratemaps ( hvl X phz) >>> -------------------------------------------
for phzI = 1:phzN
    for hvlI = 1:hvlN
        dm = res(  ...
              within_ranges(hvfl(res,2),  bins.phz.edges([hvlI, hvlI+1]))  ...
            & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
            & within_ranges(sfet,      sedges([s,s+1]))                ...
        );
        if ~isempty(dm)
            % >>> compute ego ratemap >>>
            scnt = hist2(sq(ego(dm,:)),                                     ...
                         linspace(-500,500,xbins),                          ...
                         linspace(-500,500,ybins));
            rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ecnt_hvl{hvlI},1.2));
            rmap = rmap(:);
            rmap(mazeMaskEgo(:)) = nan;
            rmap(ecnt(:)<(samplerate.*0.02)) = nan;
            rmap = reshape(rmap,[numel(fb),numel(lb)]);
            % <<< compute ego ratemap <<<
            sax = setup_axes(hfig, (phzN-phzI+1+3), 0, hvlI+4, 0, 0, 4.25, 1, 1);
            imagescnan({fb, lb, fliplr(rot90(rmap'.*samplerate .*3',-1))}, ...
                       'colorbarIsRequired', true, ...
                       'colorMap', @jet, ...
                       'colorMapFlipFlag',false);
            axis('xy');
            xlim([-300,300]);
            ylim([-300,300]);
            sax.XTickLabel = {};
            sax.YTickLabel = {};
            Lines([],0,'w');
            Lines(0,[],'w');
        end
    end
end
% <<< EGO ratemaps ( hvl X phz) <<< -------------------------------------------
% >>> EGO ratemaps ( hba X phz) >>> -------------------------------------------
for phzI = 1:phzN
    for hbaI = 1:hbaN
        dm = res(  ...
              within_ranges(hba(res),  bins.hba.edges([hbaI, hbaI+1]))  ...
            & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
            & within_ranges(sfet,      sedges([s,s+1]))                ...
        );
        if ~isempty(dm)
            % >>> compute ego ratemap >>>

            scnt = hist2(sq(ego(dm,:)),                                     ...
                         linspace(-500,500,xbins),                          ...
                         linspace(-500,500,ybins));
            rmap =(imgaussfilt(scnt,1.2)./imgaussfilt(ecnt_hba{hbaI},1.2));
            rmap = rmap(:);
            rmap(mazeMaskEgo(:)) = nan;
            rmap(ecnt(:)<(samplerate.*0.02)) = nan;
            rmap = reshape(rmap,[numel(fb),numel(lb)]);

            % <<< compute ego ratemap <<<
            sax = setup_axes(hfig, (phzN-phzI+1+3), 0, hbaI+1, 0, 0, 4.25, 1, 1);
            imagescnan({fb, lb, fliplr(rot90(rmap'.*samplerate .*3',-1))}, ...
                       'colorbarIsRequired', true, ...
                       'colorMap', @jet, ...
                       'colorMapFlipFlag',false);
            axis('xy');
            xlim([-300,300]);
            ylim([-300,300]);
            sax.XTickLabel = {};
            sax.YTickLabel = {};
            Lines([],0,'w');
            Lines(0,[],'w');
        end
    end
end
% <<< EGO ratemaps ( hba X phz) <<< -------------------------------------------

% <<< Figure spk feature summary <<< ------------------------------------------


% >>> MAIN FIGURE >>> -------------------------------------------------------------
tid = find_trial_index(sessionlist,'jg05-20120310.cof.all');
tid = find_trial_index(sessionlist,'ER06-20130613.cof.all');
tid = find_trial_index(sessionlist,'Ed10-20140816.cof.all');
tid = find_trial_index(sessionlist,'Ed10-20140817.cof.gnd');

% LOAD Trial
% >>> LOAD Trial >>> ----------------------------------------------------------

Trial = MTATrial.validate(sessionlist(tid));
Trial.par = MTAPar(Trial); % TODO : load in class
Trial.load('par'); % TODO : load in class
Trial.load('nq'); % TODO : load in class
Trial.spk.parent = Trial; % TODO : load in class

% <<< LOAD Trial <<< ----------------------------------------------------------
% >>> LOAD Ephys and Bhv Vars >>> ---------------------------------------------

stc = Trial.load('stc','msnn_ppsvd_raux');
spk = Trial.load('spk', samplerate, '', [], '', true);
phz = load_theta_phase(Trial,samplerate,Trial.meta.channelGroup.theta);
pft = pfs_2d_theta(Trial);
xyz = preproc_xyz(Trial,'trb',samplerate);
hba = fet_hba(Trial,samplerate);
hvfl = fet_hvfl(Trial,samplerate);
fvxy = vel(filter(copy(xyz),'ButFilter',4,2.5,'low'),'hcom',[1,2]);
roll = fet_roll(Trial,samplerate); % -:right, +:left
% VAR - har
% >>> hba roll feature >>> ---------------------------------------------------
har = copy(hba);
har.data = [hba.data, roll.data];
har.data = har.data * Vr;
har.data = har(:,1);
% <<< hba roll feature <<< ---------------------------------------------------
% VAR : zfet 
% >>> Normalized spike wave form features >>> ---------------------------------
zfet = spk.fet;
for c = unique(spk.clu)'
    zfet(spk.clu==c,:) = zscore(zfet(spk.clu==c,:));
end
% <<< Normalized spike wave form features <<< ---------------------------------

% <<< LOAD Ephys and Bhv Vars <<< ---------------------------------------------

mspk = Trial.load('spk', samplerate,'',arrayfun(@(u) u.id, unit{tid}),'', true, true);


u = 3;
unit{tid}(u)
figure,plot(pft,unit{tid}(u).id);

nchan = 8;
uid = unit{tid}(u).id
xbins = 36;
ybins = 36;
bhv_state = 'theta-groom-sit-rear-hpause';
res = spk(unit{tid}(u).id, [stc{bhv_state, samplerate}]);
cind = mspk.clu==uid&ismember(mspk.res,res);

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

% <<< COMPUTE EGO postion <<< -------------------------------------------------
% >>> COMPUTE allo occupancy map >>> ------------------------------------------
[ocnt,xb,yb] = hist2(sq(xyz([Trial.stc{bhv_state}],'hcom',[1,2])),...
                     linspace(-500,500,xbins),...
                     linspace(-500,500,ybins));
xc = mean([xb(1:end-1),xb(2:end)]');
yc = mean([xb(1:end-1),xb(2:end)]');
[W, H] = meshgrid( xb, yb);
mazeMaskAllo = sqrt(W.^2 + H.^2) > 480;
% <<< COMPUTE allo occupancy map <<< ------------------------------------------
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
% >>> COMPUTE sfet >>> --------------------------------------------------------

sfet = zfet( spk.clu==unit{tid}(u).id & ismember(spk.res,res), unit{tid}(u).fset);
fd = sqrt(sum(ego(res,:).^2,2));
[S,U,V] = svd(sfet);
sc = sfet * V;
[~,pscrI] = max(abs(corr(fd(fd<300),sc(fd<300,:))));
sfet = zscore(MedianFilter(sc(:,1),100));

% <<< COMPUTE sfet <<< --------------------------------------------------------
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
                         1, 25,                                             ...
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
rmap = (imgaussfilt( scnt, 1.2) ./ imgaussfilt( ocnt, 1.2));
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
rmap = imgaussfilt( scnt, 1.2) ./ imgaussfilt( ecnt, 1.2);
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
scatter([fd;fd], [phz(res);phz(res)+2*pi], 3, [sfet;sfet], 'Filled');
colormap(sax, 'jet');
caxis   (sax, [-2,2]);
xlim    (sax, [0,350]);
ylim    (sax, [0,4*pi]);
% <<< SCATTER (x=dist, y=phz, c=sfet) <<< -------------------------------------
% >>> CCG >>> -----------------------------------------------------------------
tpnts = [];
gpnts = [];
grps =  [];
for s = 1:3
    pid =  within_ranges(sfet, sedges([s,s+1])) ...
           & abs(ego(res,2)) < dthresh ...
           ;% = pid
    tres = res(pid);
    tpnts = cat(1, tpnts, tres);
    gpnts = cat(1, gpnts, s * ones([numel(tres),1]));
    grps = [grps,s];
end
bin_size = 2;
bin_halfwidth = 15;
[mccg,tbins,pairs] = CCG(tpnts,      gpnts,                             ...
                         bin_size,   bin_halfwidth,                     ...
                         samplerate, grps,                              ...
                         'hz2'                                          ...
                         );
sfob;
% <<< CCG <<< -----------------------------------------------------------------
% >>> PLOT CCG >>> ------------------------------------------------------------
for s = 1:numel(scntrs)
    for s2 = s:numel(scntrs)
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
    spk_waveforms = diff(mspk.spk(cind, 1:nchan, :), 1, 3);
    sax.Color = [0.75,0.75,0.75];
    for s = 1:numel(scntrs)
        ind = within_ranges(sfet, sedges([s,s+1])) & ...
              within_ranges(phz(res),  bins.phz.edges([phzI,phzI+1]));
        plot(bsxfun(                                                        ...
                 @minus,                                                    ...
                 sq(mean(spk_waveforms(ind,:,5:end-5)))'.*3,1:2500:20000),           ...
             '-','Color',clrs(s,:)                                                  ...
         );
    end
    sax.YTickLabel = {};
    ylim(sax, [-2.6e4,6e3]);
    grid(sax,'on');
end

% <<< spike waveforms by phase <<< --------------------------------------------
cmax = 5;
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
                & within_ranges(sfet,      bins.swf.edges([swfI, swfI+1]))  ...
                & within_ranges(hba(res),  bins.hba.edges([1,end]))  ...
                );

            % <<< QUERY partition indicies <<<
            if ~isempty(dm)
                % >>> COMPUTE ego ratemap >>>
                scnt = hist2(sq(ego(dm,:)),                                 ...
                             linspace(-500,500,xbins),                      ...
                             linspace(-500,500,ybins));
                
                rmap =(imgaussfilt(scnt,1.4)                                ...
                       ./imgaussfilt(ecnt_har{harI},1.4));
                rmap = rmap(:);

                occMaskEgo = ecnt_har{harI}(:) < (samplerate.*0.02)/3 ;
                
                rmap( mazeMaskEgo(:) ) = nan;
                rmap( occMaskEgo(:)  ) = nan;
                
                rmap = reshape(rmap,[numel(fb),numel(lb)]);
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

% <<< PRINT figure <<< --------------------------------------------------------


cmax = 25; 
% >>> FIGURE ego har >>>
figure();
for harI = 1:harN
    for phzI = 1:phzI
        dm = res( within_ranges(har(res),  bins.har.edges([harI, harI+ 1])) ...
                & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1])));
        subplot2(3,3,phzN-phzI+1,harI);
        % >>> COMPUTE ratemap >>>
        scnt = hist2(sq(ego(dm,:)),                                 ...
                     linspace(-500,500,xbins),                      ...
                     linspace(-500,500,ybins));
        rmap =(imgaussfilt(scnt,1.2)                                ...
               ./imgaussfilt(ecnt_har{harI},1.2));
        rmap = rmap(:);
        rmap(mazeMaskEgo(:)) = nan;
        rmap(ecnt_har{harI}(:)<(samplerate.*0.02)) = nan;
        rmap = reshape(rmap,[numel(fb),numel(lb)]);
        % <<< COMPUTE ratemap <<<
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
% <<< FIGURE ego har <<<

% >>> FIGURE ego hba >>>
figure();
for hbaI = 1:hbaN
    for phzI = 1:phzI
        dm = res( within_ranges(hba(res),  bins.hba.edges([hbaI, hbaI+ 1])) ...
                & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1])));
        subplot2(3,3,phzN-phzI+1,hbaI);
        % >>> COMPUTE ratemap >>>
        scnt = hist2(sq(ego(dm,:)),                                 ...
                     linspace(-500,500,xbins),                      ...
                     linspace(-500,500,ybins));
        rmap =(imgaussfilt(scnt,1.4)                                ...
               ./imgaussfilt(ecnt_hba{hbaI},1.4));
        rmap = rmap(:);
        rmap(mazeMaskEgo(:)) = nan;
        rmap(ecnt_hba{hbaI}(:)<(samplerate.*0.02)) = nan;
        rmap = reshape(rmap,[numel(fb),numel(lb)]);
        % <<< COMPUTE ratemap <<<
        % >>> IMAGESCNAN ratemap >>>
        imagescnan(                                                 ...
            {                                                       ...
                fb,                                                 ...
                lb,                                                 ...
                fliplr(rot90(rmap'.*samplerate .*3,-1))             ...
            },                                                      ...
            [ 0, cmax],                                               ...
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
% <<< FIGURE ego hba <<<

figure,plot3(phz(res),fd,sfet,'.')

figure
plot( sq(xyz([Trial.stc{'loc+pause&theta'}],'hcom',[1])),sq(xyz([Trial.stc{'loc+pause&theta'}],'hcom',[2])),'.');

% >>> SFET vs other vars >>> --------------------------------------------------

figure();
scatter3(randn([sum(pind),1])/10 + isi(pind),sfet(pind),fd(pind),20,ur(pind).*samplerate,'Filled');

pind = phz(res)<4.5;
figure();
scatter3(randn([sum(pind),1])/10 + isi(pind),sfet(pind),fd(pind),20,ur(pind).*samplerate,'Filled');
zlim([0,500]);
ylim([-4,4]);
colormap('jet');
caxis([0,50]);

figure();
for phzI = 1:bins.phz.count
    pind = within_ranges(phz(res),bins.phz.edges([phzI,phzI+1]));
    subplot2(bins.phz.count, 1, bins.phz.count+1-phzI, 1);
    scatter(sfet(pind),fd(pind),20,ur(pind).*samplerate,'Filled');
    ylim([0,500]);
    colormap('jet');
    caxis([0,30]);
end


pind = phz(res)<9;
figure();
scatter(isi(pind),sfet(pind),20,ur(pind).*samplerate,'Filled');
colormap('jet');
caxis([0,30]);

corr(sfet(pind&fd<550),fd(pind&fd<550))

% <<< SFET vs other vars <<< --------------------------------------------------



% TODO 
% >>> ANALYSIS - compute ratemaps >>> -----------------------------------------

% >>> ANALYSIS  VARS >>> ------------------------------------------------------
bhv_state = 'theta-groom-sit-rear-hpause';
% <<< ANALYSIS  VARS <<< ------------------------------------------------------
% LOAD Trial
% >>> LOAD Trial >>> ----------------------------------------------------------

Trial = MTATrial.validate(sessionlist(tid));
Trial.par = MTAPar(Trial); % TODO : load in class
Trial.load('par'); % TODO : load in class
Trial.load('nq'); % TODO : load in class
Trial.spk.parent = Trial; % TODO : load in class

% <<< LOAD Trial <<< ----------------------------------------------------------
% >>> LOAD Ephys and Bhv Vars >>> ---------------------------------------------

stc = Trial.load('stc','msnn_ppsvd_raux');
spk = Trial.load('spk', samplerate, '', [], '', true);
phz = load_theta_phase(Trial,samplerate,Trial.meta.channelGroup.theta);
pft = pfs_2d_theta(Trial);
xyz = preproc_xyz(Trial,'trb',samplerate);
hba = fet_hba(Trial,samplerate);
hvfl = fet_hvfl(Trial,samplerate);
fvxy = vel(filter(copy(xyz),'ButFilter',4,2.5,'low'),'hcom',[1,2]);
roll = fet_roll(Trial,samplerate); % -:right, +:left

% VAR - har
% >>> hba roll feature >>> ---------------------------------------------------

theta = pi/3;
Vr = [cos(theta),-sin(theta); sin(theta),cos(theta)];

har = copy(hba);
har.data = [hba.data, roll.data];
har.data = har.data * Vr;
har.data = har(:,1);
% <<< hba roll feature <<< ---------------------------------------------------
% VAR : zfet 
% >>> Normalized spike wave form features >>> ---------------------------------
zfet = spk.fet;
for c = unique(spk.clu)'
    zfet(spk.clu==c,:) = zscore(zfet(spk.clu==c,:));
end
% <<< Normalized spike wave form features <<< ---------------------------------

% <<< LOAD Ephys and Bhv Vars <<< ---------------------------------------------
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
% >>> COMPUTE phz, swf, hba - ratemaps >>> ------------------------------------
for swfI = 1:swfN
    for phzI = 1:phzN
        for hbaI = 1:harN
            % >>> QUERY partitioned spikes >>>

            dm = res(  ...
                  within_ranges(hba(res),  bins.hba.edges([hbaI, hbaI+1]))  ...
                & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
                & within_ranges(sfet,      bins.swf.edges([swfI, swfI+1]))  ...
                );

            % <<< QUERY partition indicies <<<
            if ~isempty(dm)
                % >>> COMPUTE ego ratemap >>>
                scnt = hist2(sq(ego(dm,:)),                                 ...
                             linspace(-500,500,xbins),                      ...
                             linspace(-500,500,ybins));
                
                rmap =(imgaussfilt(scnt,1.4)                                ...
                       ./imgaussfilt(ecnt_hba{hbaI},1.4));
                rmap = rmap(:);
                
                occMaskEgo = ecnt_hba{hbaI}(:) < (samplerate.*0.02)/3 ;
                
                rmap( mazeMaskEgo(:) ) = nan;
                rmap( occMaskEgo(:)  ) = nan;
                
                ratemaps{tid}(:, u, phzI, hbaI, swfI) = rmap;
                rmap = reshape(rmap,[numel(fb),numel(lb)]);
                % <<< COMPUTE ego ratemap <<<
            end
        end
    end
end
% <<< COMPUTE phz, swf, hba - ratemaps <<< ------------------------------------
% >>> COMPUTE phz, swf, har - ratemaps >>> ------------------------------------
for swfI = 1:swfN
    for phzI = 1:phzN
        for harI = 1:harN
            % >>> QUERY partitioned spikes >>>

            dm = res(  ...
                  within_ranges(har(res),  bins.har.edges([harI, harI+1]))  ...
                & within_ranges(phz(res),  bins.phz.edges([phzI, phzI+1]))  ...
                & within_ranges(sfet,      bins.swf.edges([swfI, swfI+1]))  ...
                & within_ranges(hba(res),  bins.hba.edges([1,end]))  ...
                );

            % <<< QUERY partition indicies <<<
            if ~isempty(dm)
                % >>> COMPUTE ego ratemap >>>

                scnt = hist2(sq(ego(dm,:)),                                 ...
                             linspace(-500,500,xbins),                      ...
                             linspace(-500,500,ybins));
                
                rmap =(imgaussfilt(scnt,1.4)                                ...
                       ./imgaussfilt(ecnt_har{harI},1.4));
                rmap = rmap(:);
                
                occMaskEgo = ecnt_har{harI}(:) < (samplerate.*0.02)/3 ;
                
                rmap( mazeMaskEgo(:) ) = nan;
                rmap( occMaskEgo(:)  ) = nan;
                
                ratemaps{tid}(:, u, phzI, harI, swfI) = rmap;
                rmap = reshape(rmap,[numel(fb),numel(lb)]);

                % <<< COMPUTE ego ratemap <<<
            end
        end
    end
end
% <<< COMPUTE phz, swf, har - ratemaps <<< ------------------------------------

% <<< ANALYSIS - compute ratemaps <<< -----------------------------------------



% >>> TFET >>> ----------------------------------------------------------------

u = 8;
uid = unit{tid}(u).id
figure();
plot(pft,uid,[],'colorbar',[],false);

[field_mrate, field_center] = pft.maxRate(uid);
ego=fet_ego(Trial,xyz,{'head_right','head_left'},field_center,theta);

mres=mspk(uid);
sind = mspk.clu==uid & ismember( mspk.res, mres);
cind =  spk.clu==uid & ismember(  spk.res, mres);

figure,histcirc(phz(mres));


swf=mspk.spk(sind,:,:);
mswf=sq(mean(swf));
figure,
plot((mswf-repmat([1:1000:8000]',[1,52]))')


sfet = zfet( cind, unit{tid}(u).fset);
fd = sqrt(sum(ego(mres,:).^2,2));
[S,U,V] = svd(sfet,0);
sc = sfet * V;
[~,pscrI] = max(abs(corr(fd(fd<300),sc(fd<300,:))));

figure,
for sid = 1:size(sc,2);
    subplot(1, size(sc,2), sid);
    plot(fd,sc(:,sid),'.');
end
sid = [1,3,4];

[Ss,Us,Vs] = svd(sc(:,sid),0);
sc = sc(:,sid) * Vs;
[~,pscrI] = max(abs(corr(fd(fd<300),sc(fd<300,:))));



figure,
for sid = 1:size(sc,2);
    subplot(1, size(sc,2), sid);
    plot(fd,sc(:,sid),'.');
end
sid = [1,3,4];


sid = 1;
sfet = zscore(MedianFilter(sc(:,sid),100));

figure
scatter(fd,phz(mres),10,sfet,'filled');
colormap('jet');
caxis([-2.5,2.5]);


xswf=sq(swf(:,3,:))';
mispk_ind=[];
mispk_max=[];
mispk_half=[];
mispk_dhalf=[];
for s=1:size(xswf,2)
    ixswf=interp1(1:52,-xswf(:,s),linspace(1,52,520));
    [mispk_max(s),mispk_ind(s)]=max(ixswf(100:end-100));
    mispk_ind(s) = mispk_ind(s)+100;
    mispk_half(s) = find(ixswf(mispk_ind(s):end)<mispk_max(s)/2,1,'first');
    mdh = find(ixswf([mispk_ind(s)+mispk_half(s)]:end)<ixswf(mispk_ind(s)+mispk_half(s))/2,1,'first');
    if ~isempty(mdh)
        mispk_dhalf(s) = mdh;
    else
        mispk_dhalf(s) = 0;
    end
end

tfet=[(mispk_half+randn([1,size(xswf,2)])/5)',...
      (mispk_dhalf+randn([1,size(xswf,2)])/5)'];
[St,Ut,Vt]=svd(tfet);
tfet=tfet*Vt;
edist=sqrt(sum(ego(mres,:).^2,2));
dind=edist<600;


isiR = log10(diff([1;mres]./samplerate+0.0001));
isiL = log10(diff([mres;1]./samplerate+0.0001));


normalization='xprob';
normalization='';

figure();
subplot(121);
hist2([repmat(tfet(:,1),[3,1]),[phz(mres)-2*pi;phz(mres);2*pi+phz(mres)]],linspace(-80,-30,15),linspace(-2*pi,4*pi,30),normalization);
subplot(122);
hist2([repmat(tfet(:,2),[3,1]),[phz(mres)-2*pi;phz(mres);2*pi+phz(mres)]],linspace(-20,20,15),linspace(-2*pi,4*pi,30),normalization);
colormap('jet');
axis('tight');

dind=edist<400;
figure();
scatter( repmat(tfet(dind,1),[3,1]),...
         [phz(mres(dind))-2*pi;phz(mres(dind));2*pi+phz(mres(dind))],...
         20,repmat(sfet(dind),[3,1]),'Filled');
colormap('jet');
ylim([-2*pi,4*pi]);
caxis([-3,3]);
xlim([-100,-20]);


% distance 

figure();
subplot(131); scatter( tfet(:,1), sfet, 15, fd, 'Filled'); caxis([5,250]); colormap(gca(), 'jet');
subplot(132); scatter( tfet(:,2), sfet, 15, fd, 'Filled'); caxis([5,250]); colormap(gca(), 'jet');
subplot(133); scatter( tfet(:,1), tfet(:,2), 15, fd, 'Filled'); caxis([5,250]); colormap(gca(), 'jet');

figure();
subplot(121); scatter( tfet(:,1), tfet(:,2), 15, fd, 'Filled'); caxis([5,250]); colormap(gca(), 'jet');
subplot(122); scatter( tfet(:,1), tfet(:,2), 15, sfet, 'Filled'); caxis([-2.5,2.5]); colormap(gca(), 'jet');

nind = nniz(fd);
[Ub,Sb,Vb] = svd(nunity([tfet(nind,1), tfet(nind,2), sfet(nind)]),0);
lfet = nunity([tfet(:,1), tfet(:,2), sfet(:)])*Vb;


figure
subplot(131); plot(lfet(:,1),log10(fd(:,1)),'.');
subplot(132); plot(lfet(:,2),log10(fd(:,1)),'.');
subplot(133); plot(lfet(:,3),log10(fd(:,1)),'.');
ForAllSubplots('xlim([-3,3]);ylim([1,3])');


lid = [1];
nind = nniz(fd);
[Uf,Sf,Vf] = svd([lfet(nind,lid), sfet(nind)],0);
ffet = [lfet(:,lid), sfet(:)] * Vf;

figure
subplot(131); plot(ffet(:,1),log10(fd(:,1)),'.');
subplot(132); plot(ffet(:,2),log10(fd(:,1)),'.');
subplot(133); plot(ffet(:,3),log10(fd(:,1)),'.');
ForAllSubplots('xlim([-3,3]);ylim([1,3])');


fid = 1;
    
figure();
subplot(221);
    scatter( ego(mres,2), ego(mres,1), 20, phz(mres),'Filled'); 
    colormap(gca(),'hsv'); caxis([0,2*pi]); xlim([-300,300]); ylim([-300,300]);
subplot(222);
    scatter( ego(mres,2), ego(mres,1), 20, ffet(:,fid),'Filled');
    colormap(gca(),'jet'); caxis(gca(),[-2.5,2.5]); xlim([-300,300]); ylim([-300,300]);
subplot(223);
    scatter( ego(mres,2), ego(mres,1), 20, sfet(:,1),'Filled');
    colormap(gca(),'jet'); caxis(gca(),[-2.5,2.5]); xlim([-300,300]); ylim([-300,300]);
subplot(224);
    scatter( ego(mres,2), ego(mres,1), 20, lfet(:,1),'Filled');
    colormap(gca(),'jet'); caxis(gca(),[-2.5,2.5]); xlim([-300,300]); ylim([-300,300]);

%EOF

figure();
nind = nniz(fd);
scatter3(lfet(:,1),lfet(:,2),lfet(:,3),10,fd(:),'filled');
caxis([0,300]);colormap('jet');

pz = phz(mres);
nind = nniz(fd);
figure,scatter(lfet(nind,lid), pz(nind,1), 10, fd(nind),'filled');
caxis([0,300]);colormap('jet');

pz = phz(mres);
nind = nniz(fd);
figure(); hold('on');
scatter(sfet(nind,1), pz(nind,1)-2*pi, 10, fd(nind),'filled');
scatter(sfet(nind,1), pz(nind,1), 10, fd(nind),'filled');
scatter(sfet(nind,1), pz(nind,1)+2*pi, 10, fd(nind),'filled');
caxis([0,300]);colormap('jet');


nind = nniz(fd);
figure();
hold('on');
scatter( lfet(nind,1), pz(nind)-2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( lfet(nind,1), pz(nind)+0*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( lfet(nind,1), pz(nind)+2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
caxis([0,300]);colormap('jet');
caxis([-3,1.5]);colormap('jet');

nind = nniz(fd);
figure();
hold('on');
scatter( ffet(nind,fid), pz(nind)-2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( ffet(nind,fid), pz(nind)+0*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( ffet(nind,fid), pz(nind)+2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
caxis([0,300]);colormap('jet');
caxis([-3,1.5]);colormap('jet');


figure();
nind = nniz(fd) & ffet(:,fid) < 0;
subplot(121);
hold('on');
scatter( log10(fd(nind)), pz(nind)-2*pi, 10, sfet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+0*pi, 10, sfet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+2*pi, 10, sfet(nind,1), 'Filled');
caxis([-2.5,2.5]); colormap('jet');
ylim([-2*pi,4*pi]);
xlim([1.4,2.8]);
nind = nniz(fd) & ffet(:,fid) > 0;
subplot(122);
hold('on');
scatter( log10(fd(nind)), pz(nind)-2*pi, 10, sfet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+0*pi, 10, sfet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+2*pi, 10, sfet(nind,1), 'Filled');
caxis([-2.5,2.5]); colormap('jet');
ylim([-2*pi,4*pi]);
xlim([1.4,2.8]);



figure();
nind = nniz(fd) & sfet(:,1) < 0;
subplot(121);
hold('on');
scatter( log10(fd(nind)), pz(nind)-2*pi, 10, ffet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+0*pi, 10, ffet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+2*pi, 10, ffet(nind,1), 'Filled');
caxis([-2.5,2.5]); colormap('jet');
ylim([-2*pi,4*pi]);
xlim([1.4,2.8]);
nind = nniz(fd) & sfet(:,1) > 0;
subplot(122);
hold('on');
scatter( log10(fd(nind)), pz(nind)-2*pi, 10, ffet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+0*pi, 10, ffet(nind,1), 'Filled');
scatter( log10(fd(nind)), pz(nind)+2*pi, 10, ffet(nind,1), 'Filled');
caxis([-2.5,2.5]); colormap('jet');
ylim([-2*pi,4*pi]);
xlim([1.4,2.8]);

figure();plot(sfet(:,1),ffet(:,fid),'.');

figure();
subplot(121);
    hold('on');
    nind = nniz(fd) & sfet(:,1) < 0;
    hist2([repmat(log10(fd(nind)),[3,1]), ...
           reshape(bsxfun(@plus,repmat(pz(nind),[1,3]),[-2*pi,0,2*pi]),[],1)], ...
           linspace( [1.4, 2.6 ,15] ), ...
           linspace( [-2*pi, 4*pi, 25]));
    axis('tight');
    Lines(2,[],'w');
subplot(122);
    hold('on');
    nind = nniz(fd) & sfet(:,1) > 0;
    hist2([repmat(log10(fd(nind)),[3,1]), ...
           reshape(bsxfun(@plus,repmat(pz(nind),[1,3]),[-2*pi,0,2*pi]),[],1)], ...
           linspace( [1.4, 2.6 ,15]), ...
          linspace( [-2*pi, 4*pi, 25]));
    axis('tight');
    Lines(2,[],'w');
    colormap('jet');
 
    
figure();
subplot(121);
    hold('on');
    nind = nniz(fd) & ffet(:,fid) < 0;
    hist2([repmat(log10(fd(nind)),[3,1]), ...
           reshape(bsxfun(@plus,repmat(pz(nind),[1,3]),[-2*pi,0,2*pi]),[],1)], ...
           linspace( [1.4, 2.6 ,15] ), ...
           linspace( [-2*pi, 4*pi, 25]));
    axis('tight');
    Lines(2,[],'w');
subplot(122);
    hold('on');
    nind = nniz(fd) & ffet(:,fid) > 0;
    hist2([repmat(log10(fd(nind)),[3,1]), ...
           reshape(bsxfun(@plus,repmat(pz(nind),[1,3]),[-2*pi,0,2*pi]),[],1)], ...
           linspace( [1.4, 2.6 ,15]), ...
          linspace( [-2*pi, 4*pi, 25]));
    axis('tight');
    Lines(2,[],'w');
    colormap('jet');

    

figure();
hold('on');
nind = nniz(fd);
scatter( log10(fd(nind)), pz(nind)-2*pi, 10, ffet(nind,fid), 'Filled');
scatter( log10(fd(nind)), pz(nind)+0*pi, 10, ffet(nind,fid), 'Filled');
scatter( log10(fd(nind)), pz(nind)+2*pi, 10, ffet(nind,fid), 'Filled');
caxis([-2.5,2.5]); colormap('jet');
ylim([-2*pi,4*pi]);
xlim([1.4,2.8]);



figure();
hold('on');
scatter( lfet(nind,1), pz(nind)-2*pi, 10, fd(nind), 'Filled');
scatter( lfet(nind,1), pz(nind)+0*pi, 10, fd(nind), 'Filled');
scatter( lfet(nind,1), pz(nind)+2*pi, 10, fd(nind), 'Filled');
caxis([0,200]);colormap('jet');


figure();
scatter3(tfet(nind,1),tfet(nind,2),sfet(nind,1),10,fd(nind),'filled');
caxis([0,300]);colormap('jet');

figure();
scatter3(tfet(nind,1),tfet(nind,2),sfet(nind,1),10,isiR(nind),'filled');
caxis([0,300]);colormap('jet');

figure();scatter( isiR(nind) + randn([sum(nind),1])/5, ffet(nind,fid), 10, fd(nind), 'Filled');
caxis([0,300]);colormap('jet');

figure();scatter( isiR(nind) + randn([sum(nind),1])/5, lfet(:,1), 10, fd(nind), 'Filled');
caxis([0,300]);colormap('jet');

figure();scatter( isiR(nind) + randn([sum(nind),1])/5, sfet(nind,1), 10, fd(nind), 'Filled');
caxis([0,300]);colormap('jet');


% >>> FIGURE far/close isiR, sfet, fd >>> -------------------------------------
figure();
zind = nind & fd>100;
bind = fd(nind)>100;
subplot(221);scatter( isiR(zind) + randn([sum(zind),1])/5, sfet(zind,1), 10, fd(zind), 'Filled');
caxis([0,300]);colormap(gca(),'jet');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd<100;
bind = fd(nind)<100;
subplot(222);scatter( isiR(zind) + randn([sum(zind),1])/5, sfet(zind,1), 10, fd(zind), 'Filled');
caxis([0,300]);colormap(gca(),'jet');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd>100;
bind = fd(nind)>100;
subplot(223);scatter( isiR(zind) + randn([sum(zind),1])/5, sfet(zind,1), 10, pz(zind), 'Filled');
caxis([0,2*pi]);colormap(gca(),'hsv');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd<100;
bind = fd(nind)<100;
subplot(224);scatter( isiR(zind) + randn([sum(zind),1])/5, sfet(zind,1), 10, pz(zind), 'Filled');
caxis([0,2*pi]);colormap(gca(),'hsv');
xlim([-3,2]);ylim([-4,4]);
grid('on')
% <<< FIGURE far/close isiR, sfet, fd <<< -------------------------------------

% >>> FIGURE far/close isiR, lfet, fd >>> -------------------------------------
figure();
nind = nniz(fd);
zind = nind & fd>100;
bind = fd(nind)>100;
subplot(221);scatter( isiR(zind) + randn([sum(zind),1])/5, lfet(bind,1), 10, fd(zind), 'Filled');
caxis([0,200]);colormap(gca(),'jet');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd<100;
bind = fd(nind)<100;
subplot(222);scatter( isiR(zind) + randn([sum(zind),1])/5, lfet(bind,1), 10, fd(zind), 'Filled');
caxis([0,200]);colormap(gca(),'jet');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd>100;
bind = fd(nind)>100;
subplot(223);scatter( isiR(zind) + randn([sum(zind),1])/5, lfet(bind,1), 10, pz(zind), 'Filled');
caxis([0,2*pi]);colormap(gca(),'hsv');
xlim([-3,2]);ylim([-4,4]);
grid('on')
zind = nind & fd<100;
bind = fd(nind)<100;
subplot(224);scatter( isiR(zind) + randn([sum(zind),1])/5, lfet(bind,1), 10, pz(zind), 'Filled');
caxis([0,2*pi]);colormap(gca(),'hsv');
xlim([-3,2]);ylim([-4,4]);
grid('on')
% <<< FIGURE far/close isiR, lfet, fd <<< -------------------------------------

% >>> FIGURE far/close isiR, ffet, fd >>> -------------------------------------
figure();
nind = nniz(fd);

subplot(221);
    zind = nind & fd>100;
    scatter( isiR(zind) + randn([sum(zind),1])/5, ffet(zind,1), 10, fd(zind), 'Filled');
    caxis([0,200]);colormap(gca(),'jet');
    xlim([-3,2]);ylim([-4,4]);
    grid('on')
subplot(222);
    zind = nind & fd<100;
    scatter( isiR(zind) + randn([sum(zind),1])/5, ffet(zind,1), 10, fd(zind), 'Filled');
    colormap(gca(),'jet'); caxis([0,200]);
    xlim([-3,2]);ylim([-4,4]);
    grid('on')
subplot(223);
    zind = nind & fd>100;
    scatter( isiR(zind) + randn([sum(zind),1])/5, ffet(zind,1), 10, pz(zind), 'Filled');
    caxis([0,2*pi]);colormap(gca(),'hsv');
    xlim([-3,2]);ylim([-4,4]);
    grid('on');
subplot(224);
    zind = nind & fd<100;
    scatter( isiR(zind) + randn([sum(zind),1])/5, ffet(zind,1), 10, pz(zind), 'Filled');
    caxis([0,2*pi]);colormap(gca(),'hsv');
    xlim([-3,2]);ylim([-4,4]);
    grid('on')
% <<< FIGURE far/close isiR, ffet, fd <<< -------------------------------------




figure();
hold('on');
nind = nniz(fd);
scatter( ffet(nind,1), pz(nind)-2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( ffet(nind,1), pz(nind)+0*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( ffet(nind,1), pz(nind)+2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
caxis([-3,1]);colormap('jet');
xlim([-4,4]); ylim([-2*pi,4*pi]);

figure();
subplot(121);
hold('on');
nind = nniz(fd);
nind = nind & fd<100;
scatter( sfet(nind,1), pz(nind)-2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( sfet(nind,1), pz(nind)+0*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( sfet(nind,1), pz(nind)+2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
xlim([-3,3]);ylim([-2*pi,2*pi]);caxis([-3,1.5]);colormap('jet');grid('on');
subplot(122);
hold('on');
nind = nniz(fd);
nind = nind & fd>100;
scatter( sfet(nind,1), pz(nind)-2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( sfet(nind,1), pz(nind)+0*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
scatter( sfet(nind,1), pz(nind)+2*pi, 10, isiR(nind) + randn([sum(nind),1])/5, 'Filled');
xlim([-3,3]);ylim([-2*pi,2*pi]);caxis([-3,1.5]);colormap('jet');grid('on');



figure();
subplot(121);
nind = nniz(fd) & fd<100;
scatter(lfet(nind,1),lfet(nind,2),10,fd(nind),'filled');
xlim([-3,3]);ylim([-3,3]);
grid('on');
caxis([0,300]);colormap('jet');
subplot(122);
nind = nniz(fd) & fd>100;
scatter(lfet(nind,1),lfet(nind,2),10,fd(nind),'filled');
xlim([-3,3]);ylim([-3,3]);
grid('on');
caxis([0,300]);colormap('jet');



figure();
subplot(121);
nind = nniz(fd) & fd<100;
hist2([lfet(nind,1),lfet(nind,2)],linspace(-3,3,20),linspace(-3,3,20));
Lines([],0,'w');Lines(0,[],'w');
colormap('jet');
subplot(122);
nind = nniz(fd) & fd>100;
hist2([lfet(nind,1),lfet(nind,2)],linspace(-3,3,20),linspace(-3,3,20));
Lines([],0,'w');Lines(0,[],'w');
colormap('jet');


figure();
hist2( lfet(:,1:2) * [cos(pi/6), -sin(pi/6); sin(pi/6), cos(pi/6)], linspace(-3,3,20), linspace(-3,3,20));




figure();
subplot(121);
nind = nniz(fd) & fd<100;
hist2([tfet(nind,2), sfet(nind,1)],linspace(-40,40,20),linspace(-3,3,20));
colormap('jet');
Lines([],0,'w');Lines(0,[],'w');
subplot(122);
nind = nniz(fd) & fd>100;
hist2([tfet(nind,2), sfet(nind,1)],linspace(-40,40,20),linspace(-3,3,20));
colormap('jet');
Lines([],0,'w');Lines(0,[],'w');






