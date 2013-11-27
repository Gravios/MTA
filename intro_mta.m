

%% Section 1 - Load Session

% load session                    
Session = MTASession('jg05-20120310');


% Create Session from scratch
Session = MTASession('jg05-201203XX',... name
                     'cof',       ... mazeName
                     1,           ... overwrite
                     '0x0040',    ... TTLValue found in .all.evt
                     119.881035,  ... xyzSampleRate
                     'head_front',... trackingMarker 
                     [],          ... spath 
                     32552);      ... sampleRate

% $$$ Session = 
% $$$ 
% $$$   MTASession
% $$$ 
% $$$   Properties:
% $$$               path: [1x1 struct]
% $$$              spath: [1x1 struct]
% $$$           filebase: 'jg05-20120310.cof.all'
% $$$               name: 'jg05-20120310'
% $$$          trialName: 'all'
% $$$               Maze: [1x1 MTAMaze]
% $$$              Model: [1x1 MTAModel]
% $$$                xyz: [1028817x9x3 double]
% $$$       xyzSegLength: [1x12 double]
% $$$         xyzPeriods: [12x2 double]
% $$$      xyzSampleRate: 119.8810
% $$$        syncPeriods: [12x2 double]
% $$$         sampleRate: 32552
% $$$      lfpSampleRate: 1250
% $$$     trackingMarker: 'head_front'
% $$$                ang: []
% $$$                Bhv: [1x1 MTABhv]
% $$$                Pfs: []
% $$$ 
% $$$ Session.path = 
% $$$ 
% $$$          root: '/lustre/home/zdv/kn/knajg01/data/'
% $$$       MTAPath: '/lustre/home/zdv/kn/knajg01/data/config/MTA/'
% $$$      analysis: '/lustre/home/zdv/kn/knajg01/data/analysis/'
% $$$           nlx: '/lustre/home/zdv/kn/knajg01/data/nlx/'
% $$$           xyz: '/lustre/home/zdv/kn/knajg01/data/xyz/'
% $$$       cheetah: '/lustre/home/zdv/kn/knajg01/data/config/cheetah/'
% $$$     blackrock: '/lustre/home/zdv/kn/knajg01/data/config/blackrock/'
% $$$ 
% $$$ Session.spath = 
% $$$ 
% $$$     analysis: '/lustre/home/zdv/kn/knajg01/data/analysis/jg05-20120310/'
% $$$          nlx: '/lustre/home/zdv/kn/knajg01/data/nlx/jg05-20120310/'
% $$$          xyz: '/lustre/home/zdv/kn/knajg01/data/xyz/jg05-20120310/'
                     
                     




%% Section 2 - Loading Stuff

% Load LFP
channels = 65:96;
lfp = Session.loadlfp(channels);

% Load Units
[Res,Clu,Map] = Session.loadCluRes();

% Calulate all marker angles
Session.ang = Session.markerAngles(); % last dim is direction,pitch,distance

% Create center of mass estimates from marker subsets
rbcom = zeros(size(Session.xyz,1),2,3);
rb = Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
rbcom(:,1,:) = sq(Session.com(rb));
rb = Session.Model.rb({'head_back','head_left','head_front','head_right'});
rbcom(:,2,:) = sq(Session.com(rb));



%% Section 3 - BHV

% Auto Labeling needs a xml file in analysis e.g. jg05-201203XX.cof.all.bhv.auto.xml
Session = Session.autoLableBhv('auto');

% Calulate individual features - See MTA/heuristics/rear.m for example
[~,rfet] = rear(Session,'com');

% Get BHV state periods with xyzSampleRate
rper = Session.Bhv.getState('rear').state;

% Get State position at times rper
rpos = sq(Session.xyz(rper(:,1),Session.Model.gmi(Session.trackingMarker),[1,2]));



%% Section 4 - Place Fields

% Calculate Place Fields, automatically saves a pf file
Pfs = pf(Session,[],'rear',1);

% Load Place Fields
Session = Session.loadPfs;

% Search struct for selecting place field maps
pf_search.mazeName = 'cof';
pf_search.trialName = Session.trialName;
pf_search.trackingMarker = Session.trackingMarker;
pf_search.stateLabel = 'head';
pf_search.spk_shuffle = 'n';
pf_search.pos_shuffle = 0;
pf_search.numBSiterations = 1;
pf_search.numZslices = 1;
pf_search.nbins = 50;
pf_search.smooth = 0.03;

% Get search item
Pfs = Session.getPfs(pf_search);

