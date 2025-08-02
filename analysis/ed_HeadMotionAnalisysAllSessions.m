% This script is for the analisys of the Head motion (Old script)
% Using the detected motion of the head and the relationship with the
% respiration.
% Inside we can generate the plots for figure 5 of the Behavior paper

%% Initial parameters and folder locations
%% Loading data
%% Select a specific trial
%% Low pass filtering the xyz data before angle calculations
%% Select what behaviors to include in the analysis
%% Calculating  Variables for all the sessions
%% Movement Detection Points
%%% Initialize variables
%% Bivariate Histogram Counts (Separated Detection)
%% Pitch Azimuth Head Displacement Summary Visualization only (Separated Detection)
%% Angular and Head Displacement Histograms Summary (Separated Detection)
%% Angular and Head Displacement Histograms Summary (Low speed Detected Motion)
%% Resume JPDF Azimuth and Pitch movements 
%% Spherical distance of Azimuth and Pitch Pair points to get the 2D Displacement (Low speed Detected Motion)
%% Resume plotting for spherical distances
%% Spherical distance Statistics for all events
%% Pressure sensor Phase relationship with movement (Low speed Detected)
%% Plotting 2D histograms Pressure sensor Phase
%% 1D Histograms for X and Y independent
%% Direction of motion of Azimuth and Pitch Pair points
%% Statistics of polar histograms
%% Statistics for direction


clear all
close all
%% Initial parameters and folder locations
Trim = 10;
% V_sr = 119.881035;
cd /storage/eduardo/data/Analysis/Behavior_paper/Anton-Eduardo-rhm/NewData
% corrected data
% cd /storage/eduardo/data/Analysis/Behavior_paper/Anton-Eduardo-rhm/Head-center-corrected/Head-corrected-Spline-markers/
%% Loading data

files = dir('*.mat*');
XYZ=Summary('xyz');
RHM=Summary('rhm');
PS =Summary('ps');
BEHA = Summary('Beha');
V_SR=Summary('V_sr');
Labels=Summary('Label');
Data = struct('Session',[],'gSamples',[],'TPeriods',[],'xyz',[],'rhm',[],'ps',[],'Beha', [],'V_sr',[],'Labels',[]);

for s=1: size(XYZ,2);
    [t1,t2,t3] = RemovNonVicon(XYZ{s},Trim,V_SR(s));
    Data.Session {s}=files(s).name;
    Data.gSamples {s} =t1;
    Data.TPeriods {s}=t2;
    Data.xyz {s}=t3;
    Data.V_sr{s}=V_SR(s);
    Data.Labels{s}=Labels(s,:);
    Data.rhm{s}=RHM{s}(t1);
    Data.ps{s}=PS{s}(t1);
    Data.Beha{s}=BEHA{s}(t1,:);
end

clear 't1' 't2' 't3' 's' 'XYZ' 'RHM' 'PS'

%% Select a specific trial

% Trial=GetTrial(Data,1);
Trial=Data;

%% Low pass filtering the xyz data before angle calculations


temp = {};
for ii=1:size(Trial.xyz,2);
    
    tmp=ButFilter(Trial.xyz{ii},[4],[20/(Trial.V_sr{ii}/2)],'low');
    Trial.xyz{ii}=tmp;
end

%% Select what behaviors to include in the analysis
% Acording to the order of labels by using the first 4 will include:
% pause, walk, rear, turn (the more likelly active behaviors)
% GoodP are all behavior of interest minus errors in the data

GoodP={};
BehaP={};
ErrorP={};

% Add all the data during the first four behaviours for all sessions
BehaP = cellfun(@(x) sum(x(:,[1]),2),Trial.Beha,'UniformOutput',false);
ErrorP = cellfun(@(x) sum(x(:,[7]),2),Trial.Beha,'UniformOutput',false);

GoodP = cellfun(@minus,BehaP,ErrorP,'UniformOutput',false);
% All data points to binary format
GoodP = cellfun(@(x) x>0,GoodP,'UniformOutput',false);

%% Calculating  Variables for all the sessions


% data summary
SAngles = {};
SNxyz = {};
SHeadMov = {};
SCentroid = {};
SHead_Vectors = {};
T={};
SFastAngles={};
STPitch={};


% SmAnglesFh={};
SFastHeadMotion = {};
SHeadMotion_Speed = {};
% STrans={};
% SmTransMF={};
SFastHead_Vectors={};


for ii=1:size(Trial.xyz,2);
    
    % variables
    [Angles,Nxyz,TPitch]=Head_body_motion(Trial.xyz{ii},Trial.V_sr{ii});
    [HeadMov,Centroid,Head_Vectors] = Head_mov(Trial.xyz{ii},Trial.V_sr{ii});
    t = linspace(0,length(Angles)/Trial.V_sr{ii},length(Angles));
    
    % Filtering the xyz data of the body-head relationship to get only the fast
    % components
    FastAngles = ButFilter(unwrap(Angles),[4],[5/(Trial.V_sr{ii}/2) 25/(Trial.V_sr{ii}/2)],'bandpass');
    FastHead_Vectors=ButFilter(Head_Vectors,[4],[5/(Trial.V_sr{ii}/2),25/(Trial.V_sr{ii}/2) ],'bandpass');
    
    
    %!!! The motion is detected over "FastHeadMotion", is posible to do it in the
    % xyz data or the Nxyz.
    
    FastHeadMotion=ButFilter(Nxyz,[4],[5/(Trial.V_sr{ii}/2) 20/(Trial.V_sr{ii}/2)],'bandpass');
    SlowHeadMotion=ButFilter(Nxyz,[4],[.5/(Trial.V_sr{ii}/2) 4/(Trial.V_sr{ii}/2)],'bandpass');
    

    dt = 1/ Trial.V_sr{ii};
    FastHeadMotion_Speed = sqrt(sum(diff(FastHeadMotion).^2,2))/dt;
    HeadMotion_Speed = sqrt(sum(diff(Nxyz).^2,2))/dt;
    

    %Summary Variables
    STPitch{ii}=TPitch;

    SAngles {ii}= Angles; % Non filtered angles
    SNxyz {ii}= Nxyz; % Non filtered xyz data (Head Body difference)
    SHeadMov {ii}= HeadMov; % Head markers after Head centroid substraction (2.5 low pass filtered)
    SCentroid {ii}= Centroid; % Head Centroid low pass filtered (2.5 low)
    SHead_Vectors {ii} = Head_Vectors; %  up, back and front Vector magnitude (band pass 3-20)
    SHeadMotion_Speed {ii}= HeadMotion_Speed;% Speed of Nxyz data (Non filtered!!)
    T{ii}=t; % Time
%     STrans{ii}= Trans;
%     STransM{ii}=TransM;
%     STransMF{ii}= TransMF;
    
    
    % data filtered
    SFastAngles {ii} =FastAngles; % Angles filtered 5-20
    SFastHeadMotion {ii}= FastHeadMotion; %  Nxyz filtered data (5-20)
    SFastHeadMotion_Speed {ii}= FastHeadMotion_Speed;% Derivative of the filtered Nxyz (5-20 bandpass filtered)
    SFastHead_Vectors{ii}=FastHead_Vectors;
    % Uniform distributed data
%     SmAnglesFh{ii}=mAnglesFh; % Uniform distributed angles
%     SmTransMF{ii}=mud(TransMF);
    
end

clear 'Angles' 'Nxyz' 'HeadMov' 'Centroid' 'Head_Vectors' 'AnglesFh' 'AnglesFl' 'Head_VectorsF' 't' 'xyzf' 'i' 'AnglesF'

%% Movement Detection Points
% In here we extract both:
% The Peak and throughts of the angles and Head vector Independently
% The lower speed points in 3d
% NOTE that the Angles variable are only unwraped and not filtered to avoid any distortion 

thresholdSpeed=35;
thresholdPitch=0.03;
thresholdAzimuth=0.03;
thresholdHeadVector=0.45;

%%% Initialize variables

% Detected motion points with low speed
SmovEnd={};
SmovEndF={};

%
SPitchStepLS={};
SAzimuthStepLS={};
SHvStepLS={};

% Variables from periods of Behavior of interest
SgT={};
SgAng={};
SgAngF={};
SgHv={};
SgPs={};
SgTPitch={};

% variables for separated detection
SPitchStep={};
STimeStepPitch={};

STimeStepAzimuth={};
SAzimuthStep={};

STimeStepHv={};
SHvStep={};


SPairHvIndex={};
SPairAzimuthIndex={};
SPairPitchIndex={};
SpairPsIndex={};

%Random points
SRanStepA={};
SRanStepP={};
SRanStepH={};


for ii=1:size(Trial.xyz,2);
    
    
    % to have the correct length vector variable we added 0 at the end
    VectSpeed = [ SHeadMotion_Speed{ii};0];
    VectSpeedF = [ SFastHeadMotion_Speed{ii};0];
    
    
    % localize the ends of motion (lower speed) for the NON-filtered signal
    [pk, movEnd] = findpeaks(-VectSpeed(GoodP{ii},1),'MinPeakProminence',thresholdSpeed);
    
    
    % localize the ends of motion (lower speed) for the filtered signal
    
    [pk, movEndF] = findpeaks(-VectSpeedF(GoodP{ii},1),'MinPeakProminence',thresholdSpeed);
    
    
    
    
    
    %%%%%%
    % Pitch
    
    [val,PPitch,prom ]=findpeaks(SFastAngles{ii}(GoodP{ii},2),'MinPeakProminence',thresholdPitch);
    [val,TPitch,prom ]=findpeaks(-SFastAngles{ii}(GoodP{ii},2),'MinPeakProminence',thresholdPitch);
    
    
    % Azimuth
    [val,PAzimuth,prom ]=findpeaks(SFastAngles{ii}(GoodP{ii},1),'MinPeakProminence',thresholdAzimuth);
    [val,TAzimuth,prom ]=findpeaks(-SFastAngles{ii}(GoodP{ii},1),'MinPeakProminence',thresholdAzimuth);
    
    
    % HeadVector front
    [val,PHv,prom ]=findpeaks(SFastHead_Vectors{ii}(GoodP{ii},3),'MinPeakProminence',thresholdHeadVector);
    [val,THv,prom ]=findpeaks(-SFastHead_Vectors{ii}(GoodP{ii},3),'MinPeakProminence',thresholdHeadVector);
    
    % Pressure Sensor
    pstmp= unity(ButFilter(Trial.ps{ii},[4],[20/(Trial.V_sr{ii}/2)],'low'));
    
    [val,Pps,prom ]=findpeaks(pstmp(GoodP{ii}),'MinPeakProminence',.5);
    [val,Tps,prom ]=findpeaks(-pstmp(GoodP{ii}),'MinPeakProminence',.5);
    
    
    % GoodPeriods of Behavior
    gT=T{ii}(GoodP{ii});
    gAng=(unwrap(SAngles{ii}(GoodP{ii},:)));
    gAngF=unwrap(SFastAngles{ii}(GoodP{ii},:));
    gHv=(SHead_Vectors{ii}(GoodP{ii},:));
    gPs=(pstmp(GoodP{ii}));
    gTPitch=(STPitch{ii}(GoodP{ii}));
    %%%
%     
% figure
% plot(gT,gPs)
% hold on
% plot(gT(movEndF),gPs(movEndF),'or')
% plot(gT,unity(VectSpeedF(GoodP{ii})))
    
    
    
    %%%%
    % find matching pairs
    
    pxi=[];pyi=[];axi=[];ayi=[];
    
    [p, pxi ,pyi] = NearestNeighbour(gT(PPitch),gT(TPitch),'right', Trial.V_sr{ii}*0.1);
    [a, axi ,ayi] = NearestNeighbour(gT(PAzimuth),gT(TAzimuth),'right', Trial.V_sr{ii}*0.1);
    [h, hxi ,hyi] = NearestNeighbour(gT(PHv),gT(THv),'right', Trial.V_sr{ii}*0.1);
    [P, Pxi ,Pyi] = NearestNeighbour(gT(Pps),gT(Tps),'right', Trial.V_sr{ii}*0.001);
    
    
    %%% Pitch
    PairPitchIndex = [PPitch(pxi), TPitch(pyi)];
    PairPitchTime=[gT(PairPitchIndex(:,1))', gT(PairPitchIndex(:,2))'];
    
    TimeStepPitch= abs(diff(PairPitchTime')*1000);
    PitchStep= diff([gAng(PairPitchIndex(:,1),2), gAng(PairPitchIndex(:,2),2)]');
    
    %%% Azimuth
    PairAzimuthIndex =[PAzimuth(axi) TAzimuth(ayi)];
    PairAzimuthTime=[gT(PairAzimuthIndex(:,1))' gT(PairAzimuthIndex(:,2))'];
    
    
    TimeStepAzimuth= abs(diff(PairAzimuthTime')*1000);
    AzimuthStep= diff([gAng(PairAzimuthIndex(:,1),1), gAng(PairAzimuthIndex(:,2),1)]');
    
    %%% Head vector
    PairHvIndex=[PHv(hxi), THv(hyi)];
    PairHvTime=[gT(PairHvIndex(:,1))',gT(PairHvIndex(:,2))' ];
    
    TimeStepHv = abs(diff(PairHvTime')*1000);
    HvStep = diff([gHv(PairHvIndex(:,1),3), gHv(PairHvIndex(:,2),3)]');
    
    %%% Pressure Sensor
    PairPsIndex =[Pps(Pxi) Tps(Pyi)];% Peak on Inhalation first
    PairPsTime=[gT(PairPsIndex(:,1))' gT(PairPsIndex(:,2))'];
    
    TimeStepPs= abs(diff(PairPsTime')*1000);

      %%% Movement detected from low speed
    
    PitchStepLS =diff(gAng(movEndF,2),1);
    AzimuthStepLS= diff(gAng(movEndF,1),1);
    HvStepLS=diff(gHv(movEndF,3),1);
    
    
    
    %%%%% jitter the movement points needs checking
    
    %     jitter=5;
    
    % Head vector
    % Shift1 = randi([-jitter jitter],[1,size(PairHvIndex(:,1))]); % -+5 was used because it will shift the detected points between -+41.71ms
    % movEndShift1=[PairHvIndex(:,1)+ Shift1' PairHvIndex(:,2)+Shift1'];
    %
    % movEndShift1((movEndShift1(:,2)<0 | movEndShift1(:,1)>length(gHv)),:)=[]; % remove outlayers
    % RandPoint=[gHv(movEndShift1(:,1),3), gHv(movEndShift1(:,2),3)];
    % RanStepH=diff(RandPoint'); % Random angular step
    % % tRanStep= abs(diff(T{ii}(movEndShift')))*1000; % Random time step
    %
    %
    %
    % % Pitch
    % Shift2 = randi([-jitter jitter],[1,size(PairPitchIndex(:,1))]); % -+5 was used because it will shift the detected points between -+41.71ms
    % movEndShift2=[PairPitchIndex(:,1)+Shift2' PairPitchIndex(:,2)+Shift2'];
    %
    % movEndShift2((movEndShift2(:,2)<0 | movEndShift2(:,1)>length(gAng)),:)=[]; % remove outlayers
    % RandPoint=[gAng(movEndShift2(:,1),2), gAng(movEndShift2(:,2),2)];
    % RanStepP=diff(RandPoint'); % Random angular step
    % %
    % %
    % %
    % % Azimuth
    % Shift3 = randi([-jitter jitter],[1,size(PairAzimuthIndex(:,1))]); % -+5 was used because it will shift the detected points between -+41.71ms
    % movEndShift3=[PairAzimuthIndex(:,1)+ Shift3' PairAzimuthIndex(:,2)+Shift3'];
    %
    % movEndShift3((movEndShift3(:,2)<0 | movEndShift3(:,1)>length(gAng)),:)=[]; % remove outlayers
    % RandPoint=[gAng(movEndShift3(:,1),1), gAng(movEndShift3(:,2),1)];
    % RanStepA=diff(RandPoint'); % Random angular step
    
    
    % collecting variables
    SgT{ii}=gT;
    SgAng{ii}=gAng;
    SgAngF{ii}=gAngF;
    SgHv{ii}=gHv;
    SgPs{ii}=gPs;
    SgTPitch{ii}=gTPitch;
    
    SPitchStepLS{ii}=PitchStepLS;
    SAzimuthStepLS{ii}=AzimuthStepLS;
    SHvStepLS{ii}=HvStepLS;
    
    
    SPitchStep{ii}=PitchStep;
    STimeStepPitch{ii}=TimeStepPitch;
    
    STimeStepAzimuth{ii}=TimeStepAzimuth;
    SAzimuthStep{ii}=AzimuthStep;
    
    STimeStepHv{ii}=TimeStepHv;
    SHvStep{ii}=HvStep;
    
    STimeStepPs{ii}=TimeStepPs;
    
    SPairHvIndex{ii}=PairHvIndex;
    SPairAzimuthIndex{ii}=PairAzimuthIndex;
    SPairPitchIndex{ii}=PairPitchIndex;
    SPairPsIndex{ii}=PairPsIndex;
    
    
    
    %Random points
    % SRanStepA{ii}=RanStepA;
    % SRanStepP{ii}=RanStepP;
    % SRanStepH{ii}=RanStepH;
    
    %
    SmovEnd{ii}=movEnd;
    SmovEndF{ii}=movEndF;
    
    
end

%% Bivariate Histogram Counts (Separated Detection)
edx = linspace(20,120,40);
edy = linspace(0,.5,40);
edyH = linspace(0,5,40);
OutA=[];
OutP=[];
OutH=[];

for ii=1:size(Trial.xyz,2);
    
    sig = [STimeStepAzimuth{ii}'+(randn(length(STimeStepAzimuth{ii}),1,1)*4), (abs(SAzimuthStep{ii})')];
    [out]=hist2(sig,edx,edy);
    OutA(:,:,ii)=out;
    
    sig = [STimeStepPitch{ii}'+(randn(length(STimeStepPitch{ii}),1,1)*4), (abs(SPitchStep{ii})')];
    [outp]=hist2(sig,edx,edy);
    OutP(:,:,ii)=outp;
    
    sig = [STimeStepHv{ii}'+(randn(length(STimeStepHv{ii}),1,1)*4), (abs(SHvStep{ii})')];
    [outh]=hist2(sig,edx,edyH);
    OutH(:,:,ii)= outh;
    
    
end

%% Pitch Azimuth Head Displacement Summary Visualization only (Separated Detection)


edxP = linspace(20,120,39);
edyP = linspace(0,.5,39);
edyHP = linspace(0,5,39);


figure
subplot(3,3,1)
imagesc(edx,edy,((sum(OutA,3))./sum(sum(sum(OutA,3))))')
colorbar
axis xy
title ('Azimuth')
subplot(3,3,2)
plot(edxP,bsxfun(@rdivide, sq(sum(OutA,2)),sum(sq(sum(OutA,2)))),'-o')
grid on
subplot(3,3,3)
plot(edyP,bsxfun(@rdivide, sq(sum(OutA,1)),sum(sq(sum(OutA,1)))),'-o')
grid on

%
subplot(3,3,4)
imagesc(edx,edy,((sum(OutP,3))./sum(sum(sum(OutP,3))))')
colorbar
axis xy
title ('Pitch')
subplot(3,3,5)
plot(edxP,bsxfun(@rdivide, sq(sum(OutP,2)),sum(sq(sum(OutP,2)))),'-o')
grid on

subplot(3,3,6)
plot(edyP,bsxfun(@rdivide, sq(sum(OutP,1)),sum(sq(sum(OutP,1)))),'-o')
grid on

%
subplot(3,3,7)
imagesc(edx,edyH,((sum(OutH,3))./sum(sum(sum(OutH,3))))')
colorbar
axis xy
title ('Head Vector')
subplot(3,3,8)
plot(edxP,bsxfun(@rdivide, sq(sum(OutH,2)),sum(sq(sum(OutH,2)))),'-o')
grid on

subplot(3,3,9)
plot(edyHP,bsxfun(@rdivide, sq(sum(OutH,1)),sum(sq(sum(OutH,1)))),'-o')
grid on

%% Angular and Head Displacement Histograms Summary (Separated Detection)

Az=bsxfun(@rdivide, sq(sum(OutA,1)),sum(sq(sum(OutA,1))));
Pt=bsxfun(@rdivide, sq(sum(OutP,1)),sum(sq(sum(OutP,1))));
Hv=bsxfun(@rdivide, sq(sum(OutH,1)),sum(sq(sum(OutH,1))));

figure
subplot(311)
bar(edyP,mean(Pt,2),'hist','k')
hold on
errorbar(edyP,mean(Pt,2),std(Pt,[],2),'+k')
xlim ([-.02 .4 ])
ylim ([0 .3])
xlabel('Angular Displacement (rad)')
ylabel ('Probability')
title ('Pitch')

subplot(312)
bar(edyP,mean(Az,2),'hist','k')
hold on
errorbar(edyP,mean(Az,2),std(Az,[],2),'+k')
xlim ([-.02 .4 ])
ylim ([0 .153])
xlabel('Angular Displacement (rad)')
ylabel ('Probability')
title('Azimuth')

subplot(313)
bar(edyHP,mean(Hv,2),'hist','k')
hold on
errorbar(edyHP,mean(Hv,2),std(Hv,[],2),'+k')
xlim ([-.2 4 ])
ylim ([0 .13])
xlabel('Displacement (mm)')
ylabel ('Probability')
title ('Head extension')

%% Angular and Head Displacement Histograms Summary (Low speed Detected Motion)

edy = linspace(0,.5,40);
edyH = linspace(0,5,40);

SCountsP=[];
SCountsA=[];
SCountsH=[];


for ii=1:length(Trial.xyz);
    
    
    [countsP,binsP]=histcounts(abs(SPitchStepLS{ii}),edy,'Normalization','Probability');
    [countsA,binsA]=histcounts(abs(SAzimuthStepLS{ii}),edy,'Normalization','Probability');
    [countsH,binsH]=histcounts(abs(SHvStepLS{ii}),edyH,'Normalization','Probability');
    
    SCountsP(:,ii)=countsP;
    SCountsA(:,ii)=countsA;
    SCountsH(:,ii)=countsH;
    
    
end

ypl = linspace(0,.5,39);
yHpl= linspace(0,5,39);

figure
subplot(311)
bar(ypl,mean(SCountsP,2),'hist','k')
hold on
errorbar(ypl,mean(SCountsP,2),std(SCountsP,[],2),'+k')
xlim ([-.02 .4 ])
ylim ([0 .13])
xlabel('Angular Displacement (rad)')
ylabel ('Probability')
title ('Pitch')

subplot(312)
bar(ypl,mean(SCountsA,2),'hist','k')
hold on
errorbar(ypl,mean(SCountsA,2),std(SCountsA,[],2),'+k')
xlim ([-.02 .4 ])
ylim ([0 .13])
xlabel('Angular Displacement (rad)')
ylabel ('Probability')
title ('Azimuth')

subplot(313)
bar(yHpl,mean(SCountsH,2),'hist','k')
hold on
errorbar(yHpl,mean(SCountsH,2),std(SCountsH,[],2),'+k')
xlim ([-.2 4 ])
ylim ([0 .13])
xlabel('Displacement (mm)')
ylabel ('Probability')
title ('Head extension')

%% Resume JPDF Azimuth and Pitch movements 
azinBin=linspace(min(edx),max(edx),bin-1);
pitchBin=linspace(min(edy),max(edy),bin-1);

figure
plot(pitchBin,mean(sq(sum(Sjpd,1)),2),'k')
hold on
plot(pitchBin,mean(sq(sum(Sjpd,1)),2)+(std(sq(sum(Sjpd,1)),[],2)),'--k')
plot(pitchBin,mean(sq(sum(Sjpd,1)),2)-(std(sq(sum(Sjpd,1)),[],2)),'--k')
title('Pitch')
xlabel('Pitch Normalized (rad)')
ylabel('Probability')
xlim ([-.3 1.2])

figure
plot(azinBin,mean(sq(sum(Sjpd,2)),2),'k')
hold on
plot(azinBin,mean(sq(sum(Sjpd,2)),2)+(std(sq(sum(Sjpd,2)),[],2)),'--k')
plot(azinBin,mean(sq(sum(Sjpd,2)),2)-(std(sq(sum(Sjpd,2)),[],2)),'--k')
title('Azimuth')
xlabel('Azimuth Normalized (rad)')
ylabel('Probability')

figure
imagesc(edx,edy,(Sjpd(:,:,3)/sum(sum(Sjpd(:,:,3))))')
axis xy
xlabel ('Azimuth Normalized (rad)')
ylabel('Pitch Normalized (rad)')
title(Trial.Session{3})


%% Spherical distance of Azimuth and Pitch Pair points to get the 2D Displacement (Low speed Detected Motion)
yAx = linspace(0,.5,40); % bins for histogram

SvalsG=[];
SvalsR=[];
ScountsG=[];
ScountsR=[];

SMovPairsIdx={};
SSpheDistG={};
for ii= 1: size(Trial.xyz,2);
    
    AngG=(SgAng{ii});
    AngR=flip(SgAng{ii}); % inverted signal to get random distribution
    
    % MovPairs= [SgAng{ii}(SmovEndF{ii},1), SgAng{ii}(SmovEndF{ii},2)];
    
    % Selecting arbitrary extrema points (Given that the motion is ocurring in the 3d space)
    Seq = repmat([1,2],1,length(SmovEndF{ii}));
    Seq =Seq(1:length(SmovEndF{ii}));
    
    Ini = SmovEndF{ii}(Seq==1);
    End = SmovEndF{ii}(Seq==2);
    
    % Finding pairs
    [p, pxi ,pyi] = NearestNeighbour(Ini,End,'right', Trial.V_sr{ii}*1);
    
    
    MovPairsIdx=[Ini(pxi),End(pyi)];
    
    PeakPairR = [AngR(Ini(pxi),1),AngR(Ini(pxi),2)];
    TroughPairR = [AngR(End(pyi),1),AngR(End(pyi),2)];
    
    PeakPairG = [AngG(Ini(pxi),1),AngG(Ini(pxi),2)];
    TroughPairG = [AngG(End(pyi),1),AngG(End(pyi),2)];
    
    
    %     % thresholding Pitch
    %     SpheDistG=AngleSphereDiff(PeakPairG(PeakPairG(:,2)<-.7,:),TroughPairG(PeakPairG(:,2)<-.7,:));
    %     SpheDistR=AngleSphereDiff(PeakPairR(PeakPairG(:,2)<-.7,:),TroughPairR(PeakPairG(:,2)<-.7,:));
    %
    % full Range
    SpheDistG=AngleSphereDiff((PeakPairG),(TroughPairG));
    SpheDistR=AngleSphereDiff((PeakPairR),(TroughPairR));
    
    
    %Histograms
    yAx = linspace(0,.5,40);
    
    [valsG,xx] =histcounts(SpheDistG,yAx,'Normalization','Probability');
    
    %     [valsR,xx] =histcounts(SpheDistR,yAx,'Normalization','Probability');
    
    
    SvalsG(ii,:)=valsG;
    %     SvalsR(ii,:)=valsR;
    
    SMovPairsIdx{ii}=MovPairsIdx;
    SSpheDistG{ii}=SpheDistG;
end

% yAxp = linspace(0,.5,39);
%
% figure
% % bar(yAxp,mean(SvalsG),'hist','k')
% % hold on
% errorbar(yAxp,mean(SvalsG),std(SvalsG,[],1),'-k')
% hold on
% errorbar(yAxp,mean(SvalsR),std(SvalsR,[],1),'-r')
% xlim ([-.02 .5])
% xlabel ('Spherical Distance')
% ylabel ('Probability')
% title ('Spherical distance')

%%%% Pressure sensor

for ii=1:length(SPairPsIndex);
    
    
    AngG=(SgAng{ii});
    
    PeakPairPs = [AngG(SPairPsIndex{ii}(:,1),1),AngG(SPairPsIndex{ii}(:,1),2)];
    TroughPairPs = [AngG(SPairPsIndex{ii}(:,2),1),AngG(SPairPsIndex{ii}(:,2),2)];
    SpheDistPs=AngleSphereDiff((PeakPairPs),(TroughPairPs));
    
    %%%%Permutation
    rng ('default');
    PermNum=1000;
    SpheDistPsR = zeros(length(SPairPsIndex{ii}),PermNum);
    for rr=1:PermNum;
        
        AngR=circshift(SgAng{ii},randperm(length(SgAng{ii}),1));
        
        PeakPairPsRp = [AngR(SPairPsIndex{ii}(:,1),1),AngR(SPairPsIndex{ii}(:,1),2)];
        TroughPairPsRp = [AngR(SPairPsIndex{ii}(:,2),1),AngR(SPairPsIndex{ii}(:,2),2)];
        SpheDistPsRp=AngleSphereDiff((PeakPairPsRp),(TroughPairPsRp));
        SpheDistPsR(:,rr)=SpheDistPsRp;
    end
    
    SSpheDistPsR{ii}=SpheDistPsR;
    SSpheDistPs{ii}=SpheDistPs;
    
end
%% Resume plotting for spherical distances

yAx = linspace(0,.5,40); % binning

MeanRanSphe=cellfun(@(x) mean(x,2),SSpheDistPsR,'UniformOutput',false);

% permutation values
PermSphD = reshape(cell2mat(SSpheDistPsR'),length(cell2mat(SSpheDistPsR')),PermNum);
SvalsR=zeros(length(yAx)-1,length(SSpheDistPsR));

for ss=1:length(SSpheDistPsR);
    
    tem=[];
    tem2=[];
    tem=SSpheDistPsR{ss};
    for rr=1:PermNum;
        
        [valR,xx] = histcounts(tem(:,rr),yAx,'Normalization','Probability');
        tem2(:,rr)=valR;
    end
    
SvalsR(:,ss)=mean(tem2,2);    
end


% spherical distances pressure sensor
Svals=zeros(length(yAx)-1,length(SSpheDistPs));

for ss=1:length(SSpheDistPs);
    
    [val,xx] = histcounts(SSpheDistPs{ss},yAx,'Normalization','Probability');
    Svals(:,ss)=val;
end

yAxp = linspace(0,.5,39);

%
AllEventsPs = cat(1,SSpheDistPs{:});
AllEventsMov= cat(1,SSpheDistG{:});
AllEventsRand=cat(1,PermSphD(:));

[valPs,xx] = histcounts(AllEventsPs,yAx,'Normalization','Probability');
[valMov,xx] = histcounts(AllEventsMov,yAx,'Normalization','Probability');
[valRand,xx] = histcounts(AllEventsRand,yAx,'Normalization','Probability');

figure
% Ps detected
subplot (1,2,1)
plot(yAxp,[(mean(Svals,2)+std(Svals,[],2)),(mean(Svals,2)-std(Svals,[],2)) ],'--b');
hold on
plot(yAxp,mean(Svals,2),'b')

plot(yAxp,[(mean(SvalsR,2)+std(SvalsR,[],2)),(mean(SvalsR,2)-std(SvalsR,[],2)) ],'k');
plot(yAxp,mean(SvalsR,2),'k')

% Low Speed detected points
plot(yAxp,[(mean(SvalsG',2)+std(SvalsG',[],2)),(mean(SvalsG',2)-std(SvalsG',[],2)) ],'--m');
plot(yAxp,mean(SvalsG',2),'m')


xlabel ('Spherical Distance')
ylabel ('Probability')
title ('Spherical distance from Ps Detected Points')

%
subplot (1,2,2)
plot(yAxp,cumsum(valPs));
hold on
plot(yAxp,cumsum(valMov))
plot(yAxp,cumsum(valRand))

xlabel ('Cummulative Distibution')
ylabel ('Probability')
title ('Spherical distance from Ps Detected Points')

%% Spherical distance Statistics for all events

% permiutation vs Inhalation exhalation points
[h,p,ks2stat] = kstest2(AllEventsRand(randperm(length(AllEventsRand),length(AllEventsPs))),AllEventsPs);
%h=1
%p=0
%ks2tat=0.1833

% speed detected motion vs permutation
[h,p,ks2stat] = kstest2(AllEventsRand(randperm(length(AllEventsRand),length(AllEventsPs))),AllEventsMov);
% h=1
% p=0
% ks2stat=.2238
%

%% Pressure sensor Phase relationship with movement (Low speed Detected)

BinContA=[];
BinContP=[];

BinProbP=[];
BinProbA=[];
bins=25;

for ii=1:size(SMovPairsIdx,2);
    
    % Linarized detected points
    
    % low speed pairs (see previous section) it might be not needed
    %     InxPs=reshape(f{ii},1,length(SPairPsIndex{ii})*2);
    
    
    % Pressure sensor phase
    gPs = Trial.ps{ii}(GoodP{ii},:);
    hbPs =hilbert(ButFilter(gPs,[1],[3/(Trial.V_sr{ii}/2),13/(Trial.V_sr{ii}/2)],'bandpass'));
    PhPs=(angle(hbPs));
    %     PhPs=mud(PhPs); % case of using mud transformation
    
    % Pitch and Azimuth Phase
    %     HilAng= hilbert(SgAngF{ii});
    %     PhAng=angle(HilAng);
    
    edx = linspace(-3,3,bins);
    edy = linspace(-.1,.1,bins);
    % Azimuth
    [Na]= hist2([ PhPs(SmovEndF{ii}), SgAngF{ii}(SmovEndF{ii},1)], edx,edy);
    %     Pressure sensor
    %     [Np]= hist2([ PhAng(InxPs,1), PhAng(InxPs,2)], 12,12);
    
    
    
    % Pitch
    [Np]= hist2([ PhPs(SmovEndF{ii}), SgAngF{ii}(SmovEndF{ii},2)], edx,edy);
    
    %     [Np]= hist2([ PhPs(InxPs), SgAngF{ii}(InxPs,2)], edx,edy);
    
    %%%%%Random permutation
    rng ('default');
    permNum=1000;
    RPa=zeros(length(edx)-1,length(edy)-1,permNum);
    RPp=zeros(length(edx)-1,length(edy)-1,permNum);
    
    for rr=1:permNum;
        temp= circshift(SgAngF{ii},randperm(length(SgAngF{ii}),1));
        [Ra]= hist2([ PhPs(SmovEndF{ii}), temp(SmovEndF{ii},1)], edx,edy);
        [Rp]= hist2([ PhPs(SmovEndF{ii}), temp(SmovEndF{ii},2)], edx,edy);
        
        RPa(:,:,rr)=Ra/length(SmovEndF{ii});
        RPp(:,:,rr)=Rp/length(SmovEndF{ii});
        
    end
    
    
    pVaA=zeros(length(edx)-1,length(edy)-1,size(SmovEndF{ii},2));
    pVaP=zeros(length(edx)-1,length(edy)-1,size(SmovEndF{ii},2));
    
    for rr=1:1000;
        
        pVaA = pVaA +lt(RPa(:,:,rr), Na/length(SmovEndF{ii}));
        pVaP = pVaP +lt(RPp(:,:,rr), Np/length(SmovEndF{ii}));
        
    end
    
    
    BinContA(:,:,ii)=Na/length(SmovEndF{ii});
    BinContP(:,:,ii)=Np/length(SmovEndF{ii});
    BinProbA(:,:,ii)=pVaA/permNum;
    BinProbP(:,:,ii)=pVaP/permNum;
    
end


%% Plotting 2D histograms Pressure sensor Phase
%Plotting
Xp = linspace(-3,3,bins-1);
Yp = linspace(-.1,.1,bins-1);

% Azimuth
figure
subplot(2,2,1)
imagesc(Xp,Yp,median(BinContA,3)')
axis xy
xlabel ('Pressure sensor phase ')
ylabel ('Azimuth (rad)')
title ('Ps Azimuth mean JPDF ')

subplot(2,2,2)
imagesc(Xp,Yp,median(BinProbA,3)')
axis xy
xlabel ('Pressure sensor phase ')
ylabel ('Azimuth (rad)')
title (' JPDF  vs Permutation')

% Pitch

subplot(2,2,3)
imagesc(Xp,Yp,median(BinContP,3)')
axis xy
xlabel ('Pressure sensor phase ')
ylabel ('Pitch (rad)')


subplot(2,2,4)
imagesc(Xp,Yp,median(BinProbP,3)')
axis xy
xlabel ('Pressure sensor phase ')
ylabel ('Azimuth (rad)')

% Keep in mind that different rats shows some degree of significance in the
% azimuth angles not apereant in the median plot.


%% 1D Histograms for X and Y independent
% histograms
% Pitch
figure
plot(Yp, bsxfun(@plus,std(sq(sum(BinContP,1)),[],2)* [-1 1],mean(sq(sum(BinContP,1)),2)),'--k')
hold on
plot(Yp,mean(sq(sum(BinContP,1)),2))
xlim ([-.1 .1 ])

figure
plot(Xp, bsxfun(@plus,std(sq(sum(BinContP,2)),[],2)* [-1 1],mean(sq(sum(BinContP,2)),2)),'--k')
hold on
plot(Xp,mean(sq(sum(BinContP,2)),2))
xlim ([-3 3 ])


% Azimuth
figure
plot(Yp, bsxfun(@plus,std(sq(sum(BinContA,1)),[],2)* [-1 1],mean(sq(sum(BinContA,1)),2)),'--k')
hold on
plot(Yp,mean(sq(sum(BinContA,1)),2))
xlim ([-.1 .1 ])

figure
plot(Xp, bsxfun(@plus,std(sq(sum(BinContA,2)),[],2)* [-1 1],mean(sq(sum(BinContA,2)),2)),'--k')
hold on
plot(Xp,mean(sq(sum(BinContA,2)),2))
xlim ([-3 3 ])



%% Direction of motion of Azimuth and Pitch Pair points
% High low pitch direction differences Also


%%%%
% Mov Pairs Inhanlation-> Exhalation
SMovPairsExInp={};
SMovPairsExInpF={};

% Mov Pairs Exhalation -> Inhalation
SMovPairsExIn={};
SMovPairsInExpF={};


for ii=1:size(SgAng,2);
    
    % low speed points
    %     MovPairs= [SgAng{ii}(SmovEndF{ii},1), SgAng{ii}(SmovEndF{ii},2)];
    %       MovPairsF= [SgAngF{ii}(SmovEndF{ii},1), SgAngF{ii}(SmovEndF{ii},2)];
    
    % Variable to adjust Pitch (normalized or different values)

    

    PichVar=[];
    PitchNorm=SgAng{ii}(:,2)-SlowPeakPitch(ii);
    
    %%%%% Selecting PS peaks above Pitch threshold
    PitchTreshold=.22;
    RangePitch=@(x,y) gt(x,y); %function to change depending to load values lower or grater than threshold
    PitchVar=RangePitch(PitchNorm(SPairPsIndex{ii}(:,1)),PitchTreshold); 
    %     PitchVar=':'; % Full 
    
    
    % Pressure sensor
    %%% note: SPairPsIndex has store the (peaks,trough) of pressure sensor
    %%% with pairs oriented to the right
    
    tmpExp  = [SgAng{ii}(SPairPsIndex{ii}(PitchVar,2),1), SgAng{ii}(SPairPsIndex{ii}(PitchVar,2),2)];% Peak of Exhalation
    tmpInp  = [SgAng{ii}(SPairPsIndex{ii}(PitchVar,1),1), SgAng{ii}(SPairPsIndex{ii}(PitchVar,1),2)];% Peak of inhalation
    
    tmpExpf  = [SgAngF{ii}(SPairPsIndex{ii}(PitchVar,2),1), SgAngF{ii}(SPairPsIndex{ii}(PitchVar,2),2)];% Peak of Exhalation
    tmpInpf  = [SgAngF{ii}(SPairPsIndex{ii}(PitchVar,1),1), SgAngF{ii}(SPairPsIndex{ii}(PitchVar,1),2)];% Peak of Inhalation
    
    % Inhalation->Exhalation % Direction peak of Inhalation
    MovPairsEIp= minus(tmpInp,tmpExp);
    MovPairsEIpF= minus(tmpInpf,tmpExpf);
    
    % Exhalation->Inhalation % Direction peak of Exhalation
    MovPairsIEp= minus(tmpExp,tmpInp);
    MovPairsIEpF= minus(tmpExpf,tmpInpf);
    
    %%%% Permutation
    
    rng ('default');
    PermNum=1000;
    tmpExp =[];
    tmpInp =[];
    tmpExpf=[];
    tmpInpf =[];
    MovR=zeros(length(SPairPsIndex{ii}(PitchVar,:)),2,PermNum);
    MovRf=zeros(length(SPairPsIndex{ii}(PitchVar,:)),2,PermNum);
    for rr=1: PermNum;
        
        if ~ischar(PitchVar)
            
            AngPer= circshift(repmat(SgAng{ii}(RangePitch(PitchNorm,PitchTreshold),:),4,1),randperm(length(SgAng{ii}),1));
            AngPerf= circshift(repmat(SgAngF{ii}(RangePitch(PitchNorm,PitchTreshold),:),4,1),randperm(length(SgAngF{ii}),1));
            
            tmpExp  = [AngPer(SPairPsIndex{ii}(PitchVar,2),1), AngPer(SPairPsIndex{ii}(PitchVar,2),2)];% Inhalation
            tmpInp  = [AngPer(SPairPsIndex{ii}(PitchVar,1),1), AngPer(SPairPsIndex{ii}(PitchVar,1),2)];% Exhalation
            
            tmpExpf  = [AngPerf(SPairPsIndex{ii}(PitchVar,2),1), AngPerf(SPairPsIndex{ii}(PitchVar,2),2)];% inhalation
            tmpInpf  = [AngPerf(SPairPsIndex{ii}(PitchVar,1),1), AngPerf(SPairPsIndex{ii}(PitchVar,1),2)];% Exhalation
            
            MovR(:,:,rr)= minus(tmpInp,tmpExp);
            MovRf(:,:,rr)= minus(tmpInpf,tmpExpf);
            
            
        else
            
            
            AngPer= circshift(SgAng{ii},randperm(length(SgAng{ii}),1));
            AngPerf= circshift(SgAngF{ii},randperm(length(SgAngF{ii}),1));
            
            tmpExp  = [AngPer(SPairPsIndex{ii}(:,2),1), AngPer(SPairPsIndex{ii}(:,2),2)];% Inhalation
            tmpInp  = [AngPer(SPairPsIndex{ii}(:,1),1), AngPer(SPairPsIndex{ii}(:,1),2)];% Exhalation
            
            tmpExpf  = [AngPerf(SPairPsIndex{ii}(:,2),1), AngPerf(SPairPsIndex{ii}(:,2),2)];% inhalation
            tmpInpf  = [AngPerf(SPairPsIndex{ii}(:,1),1), AngPerf(SPairPsIndex{ii}(:,1),2)];% Exhalation
            
            MovR(:,:,rr)= minus(tmpInp,tmpExp);
            MovRf(:,:,rr)= minus(tmpInpf,tmpExpf);
        end
        
    end
    
    MovPairsPr{ii}= MovR;
    MovPairsPrF{ii}= MovRf;
    
    SMovPairsExInp{ii}=MovPairsEIp;
    SMovPairsExInpF{ii}=MovPairsEIpF;
    
    SMovPairsInExp{ii}=MovPairsIEp;
    SMovPairsInExpF{ii}=MovPairsIEpF;
end

%% Direction resume
% just because the data is an sphere (see below )

% Direction resume

Bins=25;

SPolarCount=[];
SPolarCountF=[];
SMovDirection={};
SAngDisp={};
SAngDispF={};
SMovDirectionF={};

SPolarCountIEp=[];
SPolarCountIEpF=[];
SMovDirectionIE={};
SMovDirectionIEpF={};
SAngDispIEp={};
SAngDispIEpF={};

SPolarCountEIp=[];
SPolarCountEIpF=[];
SMovDirectionEIp={};
SMovDirectionEIpF={};
SAngDispEIp={};
SAngDispFEIp={};

for ii=1:length(Trial.xyz);
    
    % Taking the angular displacement for every pair of azimuth and pitch
    % (After this we can work with them as vectors)
    
    % peak Exhalation - peak Inhalation 
    AngDispEIp=(SMovPairsExInp{ii}(:,:));
    AngDispEIpF=(SMovPairsExInpF{ii}(:,:));
    
    % peak Inhalation - peak Exhalation
    AngDispIEp=(SMovPairsInExp{ii}(:,:));
    AngDispIEpF=(SMovPairsInExpF{ii}(:,:));
    
    
    
    % In case to select values of Pitch
%         AngDisp=diff(SMovPairs{ii}((SMovPairs{ii}(:,2)<-.7),:),1); % for values of pitch below .7 rad
    %     AngDispF=diff(SMovPairsF{ii}((SMovPairs{ii}(:,2)<-.7),:),1); % for values of pitch below .7 rad
    %
    %
    
    % peak Exhalation - peak Inhalation
    
    % Extracting the angle of angle pairs
    % mod is use to correct the output of atan2
    MovDirectionEIp = mod(atan2(AngDispEIp(:,2), AngDispEIp(:,1)),2*pi);
    MovDirectionEIpF = mod(atan2(AngDispEIpF(:,2), AngDispEIpF(:,1)),2*pi);
    
    % Polar counts
    [Polarbin,PolarcountEIp]=rose(MovDirectionEIp,Bins);
    [PolarbinF,PolarcountEIpF]=rose(MovDirectionEIpF,Bins);
    
    % peak Inhalation - peak exhalation
    MovDirectionIEp = mod(atan2(AngDispIEp(:,2), AngDispIEp(:,1)),2*pi);
    MovDirectionIEpF = mod(atan2(AngDispIEpF(:,2), AngDispIEpF(:,1)),2*pi);
    
    % Polar counts
    [Polarbin,PolarcountIEp]=rose(MovDirectionIEp,Bins);
    [PolarbinF,PolarcountIEpF]=rose(MovDirectionIEpF,Bins);
    
    
    
    %%%% Permutation test
    MovDirectionPr=zeros(length(PolarcountEIp),PermNum);
    MovDirectionPrf=zeros(length(PolarcountEIp),PermNum);
    
    for rr=1:PermNum;
        tempMovD=[];
        tempMovDf=[];
        tempMovD = mod(atan2(MovPairsPr{ii}(:,2,rr), MovPairsPr{ii}(:,1,rr)),2*pi);
        tempMovDf = mod(atan2(MovPairsPrF{ii}(:,2,rr), MovPairsPrF{ii}(:,1,rr)),2*pi);
        
        [Polarbin,PolarcountPr]=rose(tempMovD,Bins);
        [PolarbinF,PolarcountPrF]=rose(tempMovDf,Bins);
        
        MovDirectionPr(:,rr)= PolarcountPr/length(MovPairsPrF{ii});% Normalize all to probability
        MovDirectionPrf(:,rr)=PolarcountPrF/length(MovPairsPrF{ii});% Normalize all to probability
    end
    
    % Summary Permutation
    SMovDirectionPr{ii}=MovDirectionPr;
    SMovDirectionPrf{ii}=MovDirectionPrf;
    
    % Summary Inhalation-Exhalation
    SPolarCountEIp(ii,:)=PolarcountEIp/length(MovDirectionEIp);% Direction peak of inhalation
    SPolarCountEIpF(ii,:)=PolarcountEIpF/length(MovDirectionEIpF);
    SMovDirectionEIp{ii}=MovDirectionEIp; % Direction peak of exhalation
    SMovDirectionEIpF{ii}= MovDirectionEIpF;
    SAngDispEIp{ii}=AngDispEIp;
    SAngDispFEIp{ii}=AngDispEIpF;
    
    % Summary Exhalation-Inhalation
    SPolarCountIEp(ii,:)=PolarcountIEp/length(MovDirectionIEp);
    SPolarCountIEpF(ii,:)=PolarcountIEpF/length(MovDirectionIEpF);
    SMovDirectionIE{ii}=MovDirectionIEp;
    SMovDirectionIEpF{ii}= MovDirectionIEpF;
    SAngDispIEp{ii}=AngDispIEp;
    SAngDispIEpF{ii}=AngDispIEpF;
    
end
%% Statistics of polar histograms

% Random permutation test of normality (all are normal distributed)
% H = zeros(size(SMovDirectionPr{1},1),length(SMovDirectionPr));
% for ii=1:length(SMovDirectionPr);
%
%     for rr=size(SMovDirectionPr{ii},1);
%
%         h=kstest(SMovDirectionPr{ii}(rr,:));
%         H(rr,ii)=h;
%     end
%
% end

%% Total Resume polar plot

% chossing between fast or slow permutated data
PrVar=SMovDirectionPrf;

% chossing between fast slow or IE or EI

varEIp=SPolarCountEIpF;% Direction Peak of Inhalation
varIEp=SPolarCountIEpF;% Direction Peak of Exhalation

ResumePer = (cell2mat(cellfun(@(x) median(x,2),  PrVar,'UniformOutput',0)))';

% Statistics

ZscoreEIp=[];
ZscoreIEp=[];
for ii=1:length(Trial.xyz);
    
    
    ZscoreEIp(:,ii) = (varEIp(ii,:)'-mean(PrVar{ii},2,'omitnan'))./std(PrVar{ii},[],2,'omitnan');% Direction Peak of Inhalation
    ZscoreIEp(:,ii) = (varIEp(ii,:)'-mean(PrVar{ii},2,'omitnan'))./std(PrVar{ii},[],2,'omitnan');% Direction Peak of Exhalation
    
    
end


[x indbin]=unique(Polarbin);

%
% figure
% bar(rad2deg(Polarbin(indbin)),Zscore(indbin,:))
% hold on
% line(linspace(0,360,length(indbin)),repmat(3,1,length(Polarbin(indbin)))*1)
% line(linspace(0,360,length(indbin)),repmat(3,1,length(Polarbin(indbin)))*-1)
%
Nsize=cellfun(@(x) length(x),SMovDirectionIE);

[x,indexN]=sort(Nsize);



figure
subplot(2,2,1) % Inhalation -Exhalation

% Random distribution
% polar(Polarbin,mean(varIE,1)+std(varIE,1),'-.r')
% polar(Polarbin,mean(varIE,1)-std(varIE,1),'-.r')
polar(Polarbin,median(varEIp(1:11,:),1),'r')

hold on

% Data Distribution
% polar(Polarbin,mean(ResumePer,1)+std(ResumePer,1),':k')
% polar(Polarbin,mean(ResumePer,1)-std(ResumePer,1),':k')
polar(Polarbin,median(ResumePer,1),'k')
title ('Direction to peak of Inhalation')

subplot(2,2,2)

bar(rad2deg(Polarbin(indbin)),ZscoreEIp(indbin,1:11))
title ('Zscore from permutated distibution')

subplot(2,2,3) % Exhalation-Inhalation

% Random distribution
% polar(Polarbin,mean(varEI,1)+std(varEI,1),'-.r')
% polar(Polarbin,mean(varEI,1)-std(varEI,1),'-.r')
polar(Polarbin,median(varIEp(1:11,:),1),'r')
hold on
% Data Distribution
% % polar(Polarbin,mean(ResumePer,1)+std(ResumePer,1),':k')
% % polar(Polarbin,mean(ResumePer,1)-std(ResumePer,1),':k')
polar(Polarbin,median(ResumePer,1),'k')
title ('Direction to peak of Exhalation')

subplot(2,2,4)
bar(rad2deg(Polarbin(indbin)),ZscoreIEp(indbin,1:11))
title ('Zscore from permutated distibution')

figure
plot(rad2deg(Polarbin(indbin)),ZscoreEIp(indbin,1:11),'-b')
hold on
plot(rad2deg(Polarbin(indbin)),ZscoreIEp(indbin,1:11),'-g')
title('Resume Nose Direction')
xlabel ('Direction (degrees)')
ylabel ('Zscore (from permutation)')
ylim ([-15 27] )

%% Statistics for direction
% The zscore of the permutated data are normaly distributed acording to the Kstest
% then calculating the corrected p-val with normcdf would be fine

%  p-val
PvalEI=zeros(length(indbin),length(Trial.Session));
PvalIE=zeros(length(indbin),length(Trial.Session));

for ss=1: length(Trial.Session);
    for ii=1:length(indbin);
        PvalEI(ii,ss)= 2 * (normcdf(-abs(ZscoreEIp((indbin(ii)),ss)),0,1));% two side
        PvalIE(ii,ss)= 2 * (normcdf(-abs(ZscoreIEp((indbin(ii)),ss)),0,1));% two side
    end
end

PvalEIFDR=zeros(length(indbin),length(Trial.Session));
PvalIEFDR=zeros(length(indbin),length(Trial.Session));

crit_pVal=zeros(2,length(Trial.Session));

for  ss=1: length(Trial.Session);
    
[hfdrEI, crit_pEI, adj_pEI]=fdr_bh(PvalEI(:,ss),0.05,'pdep','yes');
[hfdrIE, crit_pIE, adj_pIE]=fdr_bh(PvalIE(:,ss),0.05,'pdep','yes');

PvalEIFDR(:,ss)=hfdrEI;
crit_pVal(1,ss)=crit_pEI;

PvalIEFDR(:,ss)=hfdrIE;
crit_pVal(2,ss)=crit_pIE;
end
PvalEIFDR(isnan(PvalEIFDR))=0;
PvalIEFDR(isnan(PvalIEFDR))=0;


%%
figure
for ii=1:length(Trial.Session);

subplotfit(ii,length(Trial.Session));
plot(rad2deg(Polarbin(indbin)),ZscoreEIp(indbin,ii),'-ob')
hold on
plot(rad2deg(Polarbin(indbin(logical(PvalEIFDR(:,ii))))),ZscoreEIp(indbin(logical(PvalEIFDR(:,ii))),ii),'*r')

plot(rad2deg(Polarbin(indbin)),ZscoreIEp(indbin,ii),'-og')
plot(rad2deg(Polarbin(indbin(logical(PvalIEFDR(:,ii))))),ZscoreIEp(indbin(logical(PvalIEFDR(:,ii))),ii),'*r')

title(Trial.Session{ii})
xlabel ('Direction (degrees)')
ylabel ('Zscore (from permutation)')
ylim ([-15 28] )
end

%% Spherical distance for each bin
% Selecting radiand bins
% RadBin=[];
% tmp=([0,Polarbin(Polarbin>0)]);
% RadBin(:,1)= tmp(logical(mod(1:length(tmp),2)));
% RadBin(:,2)= tmp(~logical(mod(1:length(tmp),2)));
% 
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% S=3;
% figure
% subplot(2,3,1)
% polar(Polarbin,SPolarCount(4,:))
% title ('Movement Direction')
% 
% subplot(2,3,2)
% compass(SAngDisp{S}(:,1),SAngDisp{S}(:,2),'.k')
% title ('Movement Direction')
% 
% 
% subplot(2,3,4)
% polar(PolarbinF,SPolarCountF(S,:))
% title ('Fast Movement Direction')
% 
% subplot(2,3,5)
% compass(SAngDispF{S}(:,1),SAngDispF{S}(:,2),'.k')
% title ('Fast Movement Direction')
% 
% % Join distribution of directions
% 
% subplot(2,3,3)
% 
% edx = linspace(0,6.2,10);
% edy = linspace(0,6.2,10);
% 
% [outTT]=hist2([(SMovDirection{S}) (SMovDirectionF{S}) ],edx,edy);
% 
% imagesc(edx,edy, outTT)
% axis xy
% xlabel 'Movement Direction'
% ylabel 'Fast Movement Direction'
% 
% 
% % History dependent motion
% subplot(2,3,6)
% hist2([SMovDirection{S},circshift(SMovDirection{S},1)],10,10)
% xlabel ('Previous motion direction')
% ylabel ('Current motion direction')
% title ('History dependence motion direction')
% 
% %% Example of movement trayectory
% 
% S=4;
% InxPs=reshape(SPairPsIndex{S},1,length(SPairPsIndex{S})*2);
% 
% % Example for few trayectories
% span=1:100;
% figure
% plot(SMovPairsInEx{S}(span,1),SMovPairsInEx{S}(span,2),'-or') % peak points
% hold on
% quiver(SMovPairsInEx{S}(span,1),SMovPairsInEx{S}(span,2),SAngDisp{S}(span,1),SAngDisp{S}(span,2),.5,'-k','LineWidth',1)
% xlabel 'Azimuth'
% ylabel 'Pitch'
% 
% 
% %% Inhalation tiggered motion
% 
% myFreqRange = [7 12];
% % finding the phase values for a given peak and trough
% ps=Trial.ps{4}(GoodP{4});
% 
% [v,ind] = findpeaks(-ps,'MinPeakProminence',1000);
% 
% 
% 
% 
% 
% tchunk=round(1/mean(diff(SgT{4})))/2; % defining the number of bins (1/2 second in this case)
% 
% time= linspace(-250,250,tchunk);
% 
% TgIn=[];
% % TgPeak = zeros(length(PairPitchIndex),tchunk);
% 
% 
% for ii= 1:length(ind);
%     
%     % fragment to extract
%     piece = ind(ii)-floor(tchunk/2): ind(ii)+floor(tchunk/2);
%     
%     if min(piece)<1 || max(piece)> length(ps)
%         
%         continue
%     else
%         TgIn(ii,:) = SgAngF{4}(piece,2);
%     end
%     
% end
% 
% % removing ceros if any
% 
% TgIn(sum(TgIn,2)==0,:)=[];
% 
% % Normalized
% TgInN= bsxfun(@rdivide,TgIn,max(TgIn,[],2));
% 
% 
% 
% 
% %% Cross correlation of Angles (Check)
% Scorr=[];
% lag=[];
% for ii=1:length(Trial.xyz)
%     %Cross correlation of the signals
%     [corr,lag] = xcorr(SgAngF{ii}(:,1:3),ceil(Trial.V_sr{ii}*1),'coeff');
%     Scorr(:,:,ii)=corr;
% end
% 
% ScorrR=[];
% lagR=[];
% for ii=1:length(Trial.xyz)
%     
%     [corrR,lagR] = xcorr(SgAngF{ii}(:,2),flip(SgAngF{ii}(:,1)),ceil(Trial.V_sr{ii}*1),'coeff');
%     
%     ScorrR(:,ii)=corrR;
% end
% 
% 
% figure
% for  ii=1:length(Trial.xyz)
%     
%     plot(lag,Scorr(:,2,ii))
%     hold on
%     grid on
% end
% 
% hold on
% errorbar(lagR,mean(ScorrR,2),std(ScorrR,[],2),'+k')
% xlabel ('Time')
% ylabel ('Correlation ')
% 
% 
% % Cross correlation of the single events
% % [acorE,lagE] = xcorr(SgAngF{S}(SmovEndF{S},1:3),ceil(Trial.V_sr{ii}*2),'coeff');
% %
% % figure
% % plot(lagE,acorE(:,2))
% % hold on
% % plot(lagE,acorE(:,8))
% 
% 
% 
% %% polar histograms (Separated Detection) not used
% 
% TpolarA=[];
% TpolarP=[];
% for ii=1:4;
%     
%     [binPa,countPa]=rose(MovDirecA{ii},50);
%     [binPp,countPp]=rose(MovDirecP{ii},50);
%     
%     TpolarA(:,:,ii)=[binPa',countPa' ];
%     TpolarP(:,:,ii)=[binPp',countPp' ];
%     
% end
% 
% mTpolarA=mean(TpolarA,3);
% mTpolarP=mean(TpolarP,3);
% 
% 
% 
% 
% figure
% polar(TpolarA(:,1),mTpolarA(:,2))
% 
% hold on
% polar(TpolarP(:,1),mTpolarP(:,2))
% 
% %%
% 
% 
% 
% % [val,Ind]=sort(gAng(:,3));
% %
% % figure
% % plot(val,gAng(Ind,3),'-')
% % hold on
% % plot(val,gAng(Ind,2),'.k')
% % plot(val,gAng(Ind,1),'.m')
% %
% % figure
% % plot(gAng(PairAzimuthIndex(:,2),1),gAng(PairAzimuthIndex(:,2),3),'.')
% 
% %%
% 
% % [val,idx]=sort(SAnglesFh{1}(:,3));
% %
% % [v,Nidx]=sort(idx);
% %
% %
% %
% % figure
% % plot(val,SAnglesFh{1}(idx,1),'.')
% % hold on
% % plot(val,SAnglesFh{1}(idx,3),'.k')
% %
% %
% %
% % test = SAnglesFh{1}(idx,1)-SAnglesFh{1}(idx,3);
% %
% %
% % Nsignal=test(Nidx);
% %
% %
% % figure
% % plot(val,Nsignal(idx),'.')
% % hold on
% % plot(val,SAnglesFh{1}(idx,3),'.k')
% 
% 
% %%
% figure
% imagesc(edxP,[],mean(bsxfun(@rdivide, sq(sum(OutP,2)),sum(sq(sum(OutP,2)))),2)')
% axis xy
% 


