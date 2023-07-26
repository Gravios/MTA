

configure_default_args();
MjgER2016_load_data();
sampleRate = 250;

unitsEgo = cell([1,numel(Trials)]);
unitsEgo{1}  = [];
unitsEgo{2}  = [];
unitsEgo{3}  = [ 80, 158, 173];                          % ER06-20130612
unitsEgo{4}  = [ 28,  35,  54,  61,  76, 107, 119, 175]; % ER06-20130613
unitsEgo{5}  = [ 31,  34,  38,  99, 112, 121, 130, 136]; % ER06-20130614
unitsEgo{6}  = [  4,   7,  10,  18,  25,  35,  37,  38,  49,  68, 104]; %Ed10-20140816
unitsEgo{7}  = [  1,  10,  24,  33,  38,  57,  63,  64,  73, 105, 108]; %Ed10-20140817
unitsEgo{8}  = [];   % jg04-20120128
unitsEgo{9}  = [];   % jg04-20120129
unitsEgo{10} = [];   % jg04-20120130
unitsEgo{11} = [];   % jg04-20120131
unitsEgo{12} = [];   % jg04-20120201
unitsEgo{13} = [];   % jg04-20120210
unitsEgo{14} = [20]; % jg04-20120211
unitsEgo{15} = [];   % jg04-20120212
unitsEgo{16} = [5];  % jg04-20120213
unitsEgo{17} = [];   % jg05-20120309
unitsEgo{18} = [ 11,  29,  33,  42,  49,  52,  54,  60,  75,  78,  80];          % jg05-20120310
unitsEgo{19} = [ 10,  13,  27,  28,  33,  50,  53,  63,  66,  67,  97, 108,  ... % jg05-20120311
                 115,142, 143, 146  153, 160, 172];
unitsEgo{20} = [ 20,  21,  25,  31,  35,  41,  61,  72,  79,  81,  85, 103, ...  % jg05-20120312
                104, 110, 111, 119, 138, 139]; 
unitsEgo{21} = [  6,  22,  24,  25,  37,  43,  44,  61,  63,  68,  77];          % jg05-20120315
unitsEgo{22} = [ 13,  19,  30,  41,  42,  48,  56,  58,  61,  65];               % jg05-20120316
unitsEgo{23} = [ 29,  50,  63,  72];           % jg05-20120317
unitsEgo{24} = [ 24,  26,  48];                % jg05-20120323
unitsEgo{25} = [ 10,  29,  55];                % jg05-20120324
unitsEgo{26} = [  4,  70,  82, 140, 145];           % ER06-20130624
unitsEgo{27} = [ 39,  51,  83,  98, 100, 102]; % Ed10-20140815
unitsEgo{28} = [];                             % er01-20110722
unitsEgo{29} = [ 11,  24,  27,  33,  49,  63,  64,  68,  70,  71,  99, 106, ...
                114, 115, 116, 132];
unitsEgo{30} = [  4,  20,  22,  23,  29,  32,  56,  63,  64,  71,  76,  83, ...
                 84,  91, 102];
unitsEgo = cf(@(T,U)  T.spk.set_unit_set(T,'egocentric',U),  Trials, unitsEgo); 



% theta phase spike histogram check phase offset
Trials{1}.meta.correction.thetaPhase  = pi;
Trials{6}.meta.correction.thetaPhase  = 2*pi/3;% 
Trials{26}.meta.correction.thetaPhase = 3*pi/4;% ER06-20130624
Trials{27}.meta.correction.thetaPhase = 2*pi/3;% Ed10-20140815
Trials{29}.meta.correction.thetaPhase = pi/2; % FS03-20201222
Trials{30}.meta.correction.thetaPhase = pi/4; % jg05-20120329


figure();
t = 4;
phz = load_theta_phase(Trials{t},sampleRate);
subplot(211);
    spkm = Trials{t}.load('spk',sampleRate,'gper&theta-groom-sit-rear',unitsInts{t},'');
    histcirc(phz(spkm.res))
subplot(212);
    spkm = Trials{t}.load('spk',sampleRate,'gper&theta-groom-sit-rear',units{t},'');
    histcirc(phz(spkm.res))
