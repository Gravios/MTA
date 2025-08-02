
% PROTOCOL - Spot session synchronization
% Generate an MTASession object for the position "spot" files 
%     ASSOCIATE - Connect Recording Periods to Trials
%        - SubjectId-YYYYMMDD.all.evt
%        - SubjectId-YYYYMMDD.pos
%
% The neuralynx-ndm evt file contains records of events.
% The primary events are the "start recording" and "stop recording"
% strings. Events during the recording periods are either a custom
% or a TTL Input string.
%
% examples:
%
% 2581488.000000 Starting Recording
% 4382605.000000 Stopping Recording
%
% 1847929.000000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0002).
% 1848039.688000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).
% 1879375.469000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0001).
% 1879485.750000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).
%
% 4458336.469000 Right reward
% 4458444.875000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).
% 4459439.562000 Left reward
% 4459536.094000 TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).
%

% ER06-20130613
% HC    -        0.000000 Starting Recording             1822013.000000 Stopping Recording
% CIR   -  1822016.000000 Starting Recording             2581476.000000 Stopping Recording
% HC    -  2581488.000000 Starting Recording             4382605.000000 Stopping Recording
% TMZ   -  4382608.000000 Starting Recording             5349154.000000 Stopping Recording
% HC    -  5349168.000000 Starting Recording             7167343.000000 Stopping Recording
% VCN   -  7167344.000000 Starting Recording      
%   COF -  7173335.500000 Vicon start recording COF      7777002.407000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
%   COF -  7783597.469000 Vicon start recording COF      8385594.375000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000). 
% VCN   -                                                8390537.000000 Stopping Recording 
% VCN   -  8390544.000000 Starting Recording
%   HC  -  8393248.063000 Vicon start recording HC       8997879.969000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
%   HC  -  9000360.032000 Vicon start recording HC       9606422.000000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
% VCN   -                                                9613103.000000 Stopping Recording
% VCN   -  9613104.000000 Starting Recording
%   COF -  9615203.032000 Vicon start recording COF     10220774.969000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).   
%   COF - 10223200.000000 Vicon start recording COF     10825206.875000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
% VCN   -                                               10830359.000000 Stopping Recording
% VCN   - 10830368.000000 Starting Recording
%   HC  - 10832934.750000 Vicon start recording HC      11437931.719000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
%   HC  - 11440441.750000 Vicon start recording HC      12046778.813000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
% VCN   -                                               12050168.000000 Stopping Recording
% VCN   - 12050176.000000 Starting Recording
%   COF - 12052055.719000 Vicon start recording COF     12661062.875000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
%   COF - 12663667.907000 Vicon start recording COF     13264945.032000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
% VCN   -                                               13269415.000000 Stopping Recording
% VCN   - 13269424.000000 Starting Recording
%   HC  - 13271718.032000 Vicon start recording HC      13877470.094000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
%   HC  - 13879650.125000 Vicon start recording HC      14481382.157000 TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).
% VCN   -                                               14485206.000000 Stopping Recording


% SETUP - Spots MTASession
% SubjectId-YYYYMMDD.pos
%    contains the 2D positions
name = 'ER06-20130613';
overwrite = true;
sync_peridos = [0.000000, 14485206.000000];
maze_name = 'cpm';
overwrite = true;
TTLValue = nan;
data_loggers = {{'nlx','spots'}};
bhv_samplerate = 50.0;

dpaths = struct('nlx',   '/storage/evgeny/data/processed/ER06/ER06-20130613/', ...
                'spots', '/storage/evgeny/data/processed/ER06/ER06-20130613/'  ...
);

link_session_Dpath( name, dpaths);

session = MTASession( name, maze_name, overwrite, TTLValue, data_loggers, bhv_samplerate);

session = MTASession.validate('ER06-20130613.cpm.all');

events = load_nlx_evt(fullfile(session.spath, [session.name '.all.evt']));

% CONCEPT MTAEvt -> MTAEvtNlx
% METHODS
%    load <- load_nlx_evt.m
%    get
%
% sync_per = events.first_and_last();
%
% sync_per = events.get('recording start','recording stop',[1,3,5]);
%
% cw_per = events.get('TTL Input on AcqSystem1_0 board 0 port 0 value (0x0001).',...
%                     'TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).',...
%                     [],'adjacent');
%
% ccw_per = events.get('TTL Input on AcqSystem1_0 board 0 port 0 value (0x0001).',...
%                      'TTL Input on AcqSystem1_0 board 0 port 0 value (0x0000).',...
%                      [],'adjacent');
%

Trial = MTATrial( name, maze_name, 'cir', overwrite, [1,1000;1100,2100])
Trial = MTATrial( name, maze_name, 'hc',  overwrite, [1,1000;1100,2100])
Trial = MTATrial( name, maze_name, 'tmz', overwrite, [1,1000;1100,2100])
Trial = MTATrial( name, maze_name, 'cof', overwrite, [1,1000;1100,2100])

session_list = get_session_list_v3('ER06-20130613-spots');


Trials = af(                                                                ...
    @( S )                                                                  ...
    MTATrial( session, session.maze.name , S.trialName),                    ...
    session_list                                                            ...
);

cf( @( T ) set(T.spk, 'parent', T), Trials);

units = session.spk.map(:,1)';


Pft = cf( @( T )  pfs_2d_theta( T, units, 'theta' ), Trials);


global MTA_FIGURES_PATH

hfig = figure();
index = 1;
groups =                                                   ...
    {                                                      ...
        struct('key', 'g', 'label',    'good', 'data',[]), ...
        struct('key', 'b', 'label',     'bad', 'data',[]), ...
        struct('key', 'c', 'label',  'center', 'data',[]), ...
        struct('key', 'e', 'label',    'edge', 'data',[]), ...
        struct('key', 'i', 'label',     'int', 'data',[]), ...
        struct('key', 'p', 'label',     'pyr', 'data',[]), ...
        struct('key', 'w', 'label',   'start', 'data',[]), ...
        struct('key', 'r', 'label',  'reward', 'data',[]), ...
        struct('key', 's', 'label', 'strange', 'data',[]), ...
    };
groups = containers.Map( cf(@(g) g.key, groups), groups);
autoIncr = false;
figname = create_directory(fullfile( MTA_FIGURES_PATH, 'unit_selection', session.name, session.name ));
while index ~= -1
    for tid = 1:numel(Pft)
        subplot(1,2,tid)
        plot( Pft{tid}, index, '', 'text', [], false);
    end
    index = figure_controls(hfig, index,  units, autoIncr, figname, 'v', groups);
end

Trial = Trials{2};
Trial.meta = session_list(2).subject;
pft = Pft{2};
samplerate = 500;
theta = 0.225;
load_features = true;
load_spk_waveforms = true;

xy = Trial.load('xyz');
xy.resample(samplerate);
phz = load_theta_phase( Trial, xy, 1);
spk = Trial.load('spk', samplerate, '', groups('g').data, '', load_features, load_spk_waveforms);

zfet = spk.fet;
for c = unique(spk.clu)'
    zfet(spk.clu==c,:) = zscore(zfet(spk.clu==c,:));
end



uid = groups('g').data(6);
figure();
plot(pft,uid,[],[],[],false);
% uid = 1; field_center = [30, 430];
field_center = [-470,-930];
%[~,field_center] = pft.maxRate(uid, false);
ego = fet_ego(Trial, xy, {'head_left','head_right'}, field_center, theta);

mres = spk(uid);

mres = mres(abs(ego(mres,2))<200);

figure,
plot(ego(mres,1),phz(mres),'.');
hold('on');
plot(ego(mres,1),2*pi+phz(mres),'.b');
plot(ego(mres,1),-2*pi+phz(mres),'.b');
xlim([-300,300]);
ylim([-2*pi,4*pi]);




cind = spk.clu==uid & ismember(spk.res,mres);
fd = sqrt(sum(ego(mres,:).^2,2));
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
        sfet = MedianFilter(sfet * V(:,pscrI),50);
    end
end


% EOF










