% EgoProCode2D_load_data - loads the folowing variables and functions
%
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      figBasePath
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfs
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector


project.name = 'EgoProCode2d';



%% GLOBAL VARIABLES ------------------------------------------------------------------------------
global MTA_FIGURES_PATH
global PROJECT_FIGURE_PATH

PROJECT_FIGURE_PATH = create_directory(getenv(['MTA_PROJECT_' upper(project.name)]));


%% LOCAL PROJECT VARIABLES -----------------------------------------------------------------------
% LOAD Trials and metadata
sessionListName = 'MjgER2016';
sessionList = get_session_list_v2(sessionListName);
Trials = af(@(s) MTATrial.validate(s), sessionList);

% LOAD place cell ids -> one cell per Trial
units = cf(@(T)  T.spk.get_unit_set(T,'placecells'),  Trials);
% REMOVE place cells with poor sampling or low rates
units = cf(@(T,U) remove_bad_units(T,U), Trials,units);

% LOAD interneuron ids -> one cell per Trial
unitsInts = cf(@(T)  T.spk.get_unit_set(T,'interneurons'),  Trials); 

% LOAD place cell ids for egocentric project -> one cell per Trial
unitsEgo = cf(@(T)  T.spk.get_unit_set(T,'egocentric'),  Trials); 


% er01-20110719
unitsEgoHvf{1} = [15,42,99];
% er01-20110721
unitsEgoHvf{2} = [23,75,86,99];
% ER06-20130612
unitsEgoHvf{3} = [38,61,80,151,158];
% ER06-20130613
unitsEgoHvf{4} = [23,33,38,49,51,69,94,107,117,126,140,145,154,166,175];
% ER06-20130614
unitsEgoHvf{5} = [15,34,35,38,40,44,86,90,99,112,121,159];
% Ed10-20140816
unitsEgoHvf{6} = [4,7,10,13,18,25,35,37,49,89];
% Ed10-20140817
unitsEgoHvf{7} = [1,10,17,23,33,38,57,63,64,66,73,82,99,104,105,108];% 2
% jg04-20120128
unitsEgoHvf{8} = [10];
% jg04-20120129
unitsEgoHvf{9} = [39];
% jg04-20120130
unitsEgoHvf{10} = [];
% jg04-20120131
unitsEgoHvf{11} = [];
% jg04-20120201
unitsEgoHvf{12} = [];
% jg04-20120210
unitsEgoHvf{13} = [];
% jg04-20120211
unitsEgoHvf{14} = [10,17,20];
% jg04-20120212
unitsEgoHvf{15} = [];
% jg04-20120213
unitsEgoHvf{16} = [5];
% jg05-20120309
unitsEgoHvf{17} = [70];
% jg05-20120310
unitsEgoHvf{18} = [9,11,17,18,19,20,21,22,24,25,29,33,34,35,42,...
                  44,49,52,56,60,61,72,74,75,78,81,82,84];
% jg05-20120311
unitsEgoHvf{19} = [10,13,15,27,63,66,67,97,143,150,152,153,160];
% jg05-20120312
unitsEgoHvf{20} = [20,21,25,31,35,41,44,52,59,61,72,79,80,81,85,86,89,103,...
                   105,110,111,119,129,139,141,144,151];
% jg05-20120315
unitsEgoHvf{21} = [6,22,24,25,37,43,61,63,68,77,97];
% jg05-20120316
unitsEgoHvf{22} = [8,11,13,18,19,22,30,41,42,43,48,50,56,58,61,65];
% jg05-20120317
unitsEgoHvf{23} = [9,12,18,29,31,36,40,50,51,54,56,63,69,70,72];
% jg05-20120323
unitsEgoHvf{24} = [22,26,32];
% jg05-20120324
unitsEgoHvf{25} = [10,29];
% ER06-20130624
unitsEgoHvf{26} = [82,140,175,218];
% Ed10-20140815
unitsEgoHvf{27} = [6,10,39,51,92,98,100,102];
% er01-20110722
unitsEgoHvf{28} = [15,45,74,90];
% FS03-2020122
unitsEgoHvf{29} = [11,12,24,27,33,43,63,64,65,69,70,71,76,116,120,126];
% jg05-20120329
unitsEgoHvf{30} = [20,23,29,30,55,56,63,82,83,84,97,102,107];
                   
                   
%location = 'CA1'
%isInCA1 = any(arrayfun(@(probe) strcmp(probe.location,location),subject.probe));
sessionsCA1 = [3:5,8:12,17:25,29];
uehCount = cellfun(@numel,unitsEgoHvf);
uehCsum = [0,cumsum(uehCount)];


unitsEgoHvfCA1 = [];
for up = [uehCsum(sessionsCA1)'+1,uehCsum(sessionsCA1+1)']'
    unitsEgoHvfCA1 = cat(1,unitsEgoHvfCA1,[up(1):up(2)]');
end


sessionsCA3 = [1,2,6,7,13:16,26:28,30];
unitsEgoHvfCA3 = [];
for up = [uehCsum(sessionsCA3)'+1,uehCsum(sessionsCA3+1)']'
    unitsEgoHvfCA3 = cat(1,unitsEgoHvfCA3,[up(1):up(2)]');
end



% SET index arry for each anatomical location in oder of cumulative
% index within cells
% $$$ unitsEgoCA1 = [1:19,44:149];
% $$$ unitsEgoCA3 = [20:41,123:127,150:164];
unitsEgoCA1 = [1:19,44:122,134:149];
unitsEgoCA3 = [20:43,123:133,150:164];





% MAP 
cluSessionMap = [];
for u = 1:numel(units)
    cluSessionMap = cat(1,cluSessionMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

egoCluSessionMap = [];
for u = 1:numel(unitsEgo)
    egoCluSessionMap = cat(1,egoCluSessionMap,[u*ones([numel(unitsEgo{u}),1]),unitsEgo{u}(:)]);
end

egoHvfCluSessionMap = [];
for u = 1:numel(unitsEgo)
    egoHvfCluSessionMap = cat(1,egoHvfCluSessionMap,[u*ones([numel(unitsEgoHvf{u}),1]),unitsEgoHvf{u}(:)]);
end


% SET the reference trial for computing the inter-subject head pitch correction 
pitchReferenceTrial = 'Ed05-20140529.ont.all';


% SET states
clear('states')
states.name = 'msnn_ppsvd_raux';
states.labels = {{'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta'}};

interpParPfs = struct('bins',{{linspace(-500,500,100)',linspace(-500,500,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

interpParDfs = struct('bins',{{linspace(-2,2,100)',linspace(-2,2,100)'}},...
                   'nanMaskThreshold', 0,...
                   'methodNanMap',     'linear',...
                   'methodRateMap',    'linear');

sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,units);    


EgoProCode2D_load_ratemaps();



hbaBin.edges = [-1.2,-0.2,0.2,1.2];
hbaBin.centers = mean([hbaBin.edges(1:end-1);hbaBin.edges(2:end)]);
hbaBin.count = numel(hbaBin.centers);
hbaBin.color = [0,1,0;...
                0,0,1;...
                1,0,0];
hbaBin.key = 'LCR';
hbaBin.label = {'Left','Center','Right'};



phzBin.edges = linspace(0.5,2*pi-0.5,4);
phzBin.centers = (phzBin.edges(1:end-1)+phzBin.edges(2:end))./2;
phzBin.count = numel(phzBin.centers);
phzBin.color = cool(3);
phzBin.key = 'DTA';
phzBin.label = {'Descending','Trough','Ascending'};


hvfBin.edges = [-25,-5,5,25,80];
hvfBin.centers = (hvfBin.edges(1:end-1)+hvfBin.edges(2:end))./2;
hvfBin.count = numel(hvfBin.centers);
hvfBin.color = bone(3);
hvfBin.key = 'RISF';
hvfBin.label = {'Reverse','Immobile','Slow','Fast'};


% Head Body Angle
bins.hba.name = 'hba';
bins.hba.description = 'Head Body Angle';
bins.hba.edges = [-1.2,-0.2,0.2,1.2];
bins.hba.centers = mean([bins.hba.edges(1:end-1);bins.hba.edges(2:end)]);
bins.hba.count = numel(bins.hba.centers);
bins.hba.color = [0,1,0;...
                0,0,1;...
                1,0,0];
bins.hba.key = 'LCR';
bins.hba.label = {'Left','Center','Right'};

% Theta Phase
bins.phz.name = 'phz';
bins.phz.description = 'CA1 Theta Phase';
bins.phz.edges = linspace(0.5,2*pi-0.5,4);
bins.phz.centers = (bins.phz.edges(1:end-1)+bins.phz.edges(2:end))./2;
bins.phz.count = numel(bins.phz.centers);
bins.phz.color = cool(3);
bins.phz.key = 'DTA';
bins.phz.label = {'Descending','Trough','Ascending'};

% Speed along the Anteroposterior Axis
bins.phz.name = 'hvf';
bins.phz.description = 'Speed along Anteroposterior Axis of the Head';
bins.hvf.edges = [-25,-5,5,25,80];
bins.hvf.centers = (bins.hvf.edges(1:end-1)+bins.hvf.edges(2:end))./2;
bins.hvf.count = numel(bins.hvf.centers);
bins.hvf.color = bone(3);
bins.hvf.key = 'RISF';
bins.hvf.label = {'Reverse','Immobile','Slow','Fast'};








fwd = 1;
lat = 2;

% Egofield offset
offset = [-2,0.5];


%% HELPER FUNCTIONS ----------------------------------------------------------------------------
% SET helper function to reshape eigenvectors
reshape_eigen_vector = @(V,p) reshape(V(:,1),p{1}.adata.binSizes')';
