function turn_state = turn(Session,varargin)

Session = MTASession('jg05-20120315');
Session = MTATrial(Session,'crt1');

[Res,Clu,Map] = LoadCluRes([Session.spath.nlx Session.name]);
Res = round(Res.*Session.lfpSampleRate./Session.sampleRate);
[Res,ind] = SelectPeriods(Res,[Session.syncPeriods(1),Session.syncPeriods(end)],'d',1,1);
Clu = Clu(ind);


fet = unwrap(sq(Session.ang(:,Session.Model.gmi('spine_upper'),Session.Model.gmi('head_back'),1)-Session.ang(:,Session.Model.gmi('head_back'),Session.Model.gmi('head_front'),1)));
%fet = unwrap(sq(Session.ang(:,Session.Model.gmi('spine_upper'),Session.Model.gmi('head_back'),2)- Session.ang(:,Session.Model.gmi('head_back'),Session.Model.gmi('head_front'),2)));
%fet = sq(Session.ang(:,Session.Model.gmi('spine_upper'),Session.Model.gmi('head_back'),3));

fet(isnan(fet)) = 0;


fetmh = ButFilter(fet,2,[5 12]/59,'bandpass');
fetmm = ButFilter(fet,2,[1 5]/59,'bandpass');
fetml = ButFilter(fet,2,[0.1 1]/59,'bandpass');
tfet  = d2t(fetml,Session.xyzSampleRate,0);
tfetv = d2t(diff(fetml),Session.xyzSampleRate,0);



%% Spike Triggered Average
sta = cell(size(Map,1),1);
for i=1:size(Map,1),
sta{

%% Use LocalMinima to find turning points 

figure,
sp4 = subplot(211);
plot(tfetv,unity(diff(fetml)),tfetv,unity(diff(fetmm)),tfetv,unity(diff(fetmh)));
Lines([],0,'k')
sp5 = subplot(212);
plot(tfet,fet);
linkaxes([sp4,sp5],'x');


hist(fetmt,20000)


fetlt = unity(diff(fetml));
fetmt = unity(diff(fetmm));
fetht = unity(diff(fetmh));

trpointsl = LocalMinima(fetlt);
trpointsm = LocalMinima(fetmt);
trpointsh = LocalMinima(fetht);

tlpointsl = LocalMinima(-fetlt);
tlpointsm = LocalMinima(-fetmt);
tlpointsh = LocalMinima(-fetht);

trintl = sum(GetSegs(abs(fetlt),trpointsl-5,10));
trintm = sum(GetSegs(abs(fetmt),trpointsm-5,10));
trinth = sum(GetSegs(abs(fetht),trpointsh-5,10));


tlintl = sum(GetSegs(abs(fetlt),tlpointsl-5,10));
tlintm = sum(GetSegs(abs(fetmt),tlpointsm-5,10));
tlinth = sum(GetSegs(abs(fetht),tlpointsh-5,10));

hist(tlintm,20000)

figure,
sp1 = subplot(211);
plot(tfet,fet);
sp2 = subplot(212);
plot(tfetv,fetlt)
hold on 
linkaxes([sp1,sp2],'x');


%% Set Thresholds
ut = 3;
lt = 2;
gprl = find(abs(trintl)>lt&abs(trintl)<ut);
gprm = find(abs(trintm)>lt&abs(trintm)<ut);
gprh = find(abs(trinth)>lt&abs(trinth)<ut);

gpll = find(abs(tlintl)>lt&abs(tlintl)<ut);
gplm = find(abs(tlintm)>lt&abs(tlintm)<ut);
gplh = find(abs(tlinth)>lt&abs(tlinth)<ut);


%% CCG of each turn
stsprl = round((trpointsl(gprl)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
stspll = round((tlpointsl(gpll)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);

stsprm = round((trpointsm(gprm)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
stsplm = round((tlpointsm(gplm)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);

stsprh = round((trpointsh(gprh)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);
stsplh = round((tlpointsh(gplh)-1)./Session.xyzSampleRate.*Session.lfpSampleRate);


[ccgl tbinl ] = Trains2CCG({stsprl, stspll, Res},{1,2,Clu},10,100,Session.lfpSampleRate,'hz');
[ccgm tbinm ] = Trains2CCG({stsprm, stsplm, Res},{1,2,Clu},4,150,Session.lfpSampleRate,'hz');
[ccgh tbinh ] = Trains2CCG({stsprh, stsplh, Res},{1,2,Clu},1,150,Session.lfpSampleRate,'hz');


uClu = unique(Clu);

agccgl = zeros(size(ccgl,1),2,size(Map,1)); 
agccgl(:,:,uClu) = sq(ccgl(:,1:2,3:end));

agccgm = zeros(size(ccgm,1),2,size(Map,1)); 
agccgm(:,:,uClu) = sq(ccgm(:,1:2,3:end));

agccgh = zeros(size(ccgh,1),2,size(Map,1)); 
agccgh(:,:,uClu) = sq(ccgh(:,1:2,3:end));




%% Display Everything                                                                                                                                                                                                                       
figure(102)
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
clf
set(gcf,'Name',num2str(unit));
subplot2(6,3,[1,3],1);
    pf_search.stateLabel = 'head';
    Pfs = Session.getPfs(pf_search);
    ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})

subplot2(6,3,[4,6],1);
bar(atbin,accg(:,unit));axis tight;grid on;

subplot2(6,3,1,[2,3]);
bar(tbinl,agccgl(:,1,unit)),axis tight;grid on;
subplot2(6,3,2,[2,3]);
bar(tbinl,agccgl(:,2,unit)),axis tight;grid on;

subplot2(6,3,3,[2,3]);
bar(tbinm,agccgm(:,1,unit)),axis tight;grid on;
subplot2(6,3,4,[2,3]);
bar(tbinm,agccgm(:,2,unit)),axis tight;grid on;

subplot2(6,3,5,[2,3]);
bar(tbinh,agccgh(:,1,unit)),axis tight;grid on;
subplot2(6,3,6,[2,3]);
bar(tbinh,agccgh(:,2,unit)),axis tight;grid on;

waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #: ');
      case double('t')
        thresh = input('Enter threshold: ');
      case double('m')
        ppf_type = input('Enter placefield plot mode: ');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
    end
end




fs = MTASession(Session.name,Session.Maze.name);
[accg atbin] = autoccg(fs);


%% Load Place Fields                                                                                                                                                                                                                        
Session.Pfs = {};
Session = Session.loadPfs;

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


%stateLabel = 'turn_slow';
stateLabel = 'turn_med';
%stateLabel = 'turn_fast';

save([Session.spath.analysis Session.filebase '.ccg.' stateLabel '.mat'],'tbin','sgccg');


