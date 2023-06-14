%% Spectrum and coherence computation for all sessions
% This Script computes the Spectra of motion variables and plots summary of
% Density spectra and lag analysis of RHM and Pressure sensor
%%
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
    Data.Labels{s}=Labels{s};
    Data.rhm{s}=RHM{s}(t1);
    Data.ps{s}=PS{s}(t1);
    Data.Beha{s}=logical(BEHA{s}(t1,:));
end

clear 't1' 't2' 't3' 's' 'XYZ' 'RHM' 'PS'

%% Select a specific trial

% Trial=GetTrial(Data,1);
Trial=Data;
%% Low pass filtering the xyz data before angle calculations


temp = {};
for ii=1:size(Trial.xyz,2);
    
    tmp=ButFilter(Trial.xyz{ii},[4],[15/(Trial.V_sr{ii}/2)],'low');
    Trial.xyz{ii}=tmp;
end
%% Select behaviors to include in the analysis
% Acording to the order of labels by using the first 4 will include:
% pause, walk, rear, turn (the more likelly active behaviors)
% GoodP are all behavior of interest minus errors in the data

GoooP={};
BehaP={};
ErrorP={};
% Add all the data during the first four behaviours for all sessions
BehaP = cellfun(@(x) double(sum(x(:,[1:6]),2)>0),Trial.Beha,'UniformOutput',false);
ErrorP = cellfun(@(x) sum(x(:,[7]),2),Trial.Beha,'UniformOutput',false);

% Keep Behaviors minus the Errors in xyz data
GoodP = cellfun(@(x,y) minus(x,y)>0,BehaP,ErrorP,'UniformOutput',false);
%% NOT USED Extracting the Behavior periods for all data

% 
% for ii = 1:size(GoodP,2);
%     
%     Trial.xyz{ii}=Trial.xyz{ii}(GoodP{ii},:,:);
%     Trial.rhm{ii}=Trial.rhm{ii}(GoodP{ii});
%     Trial.ps{ii}=Trial.ps{ii}(GoodP{ii});
%     Trial.Beha{ii}=Trial.Beha{ii}(GoodP{ii},:);
% end    


%% Variables to Analize  
%    

% data summary
% SNxyz = {};
% SHeadMov = {};
% SCentroid = {};
% SHead_Vectors = {};

% SAnglesF={};
% SAngleFh={};
% SAngleFl={};
% SHead_VectorsF={};
T={};
Svariables={};
SAngles = {};


for ii=1:size(Trial.xyz,2);
% variables
% xyzf =ButFilter(Data.xyz{ii},[2],[.1/(V_sr/2) 50/(V_sr/2)],'bandpass'); % note that the xyz raw data is filtered before the calculations

[Angles,Nxyz,x,PitchT]=Head_body_motion(Trial.xyz{ii},Trial.V_sr{ii});

[HeadMov,Centroid,Head_Vectors] = Head_mov(Trial.xyz{ii}(:,:,:),Trial.V_sr{ii});

t = linspace(0,length(Angles)/Trial.V_sr{ii},length(Angles));
    
% Filtering

% Broad filtering
AnglesF = ButFilter(unwrap(Angles),[4],[.1/(Trial.V_sr{ii}/2) 50/(Trial.V_sr{ii}/2) ],'bandpass');
% high frequency
% AnglesFh = ButFilter(unwrap(Angles),[4],[6/(Trial.V_sr{ii}/2) 20/(Trial.V_sr{ii}/2)],'bandpass');
% 
% % low frequency
% AnglesFl = ButFilter(unwrap(Angles),[4],[.1/(Trial.V_sr{ii}/2) 6/(Trial.V_sr{ii}/2)],'bandpass');

Head_VectorsF = ButFilter(Head_Vectors,[4],[.1/(Trial.V_sr{ii}/2) 50/(Trial.V_sr{ii}/2)],'bandpass');    

variables = [AnglesF  Head_VectorsF Trial.rhm{ii} Trial.ps{ii} PitchT];




SAngles {ii}= Angles;
% SNxyz {ii}= Nxyz;
% SHeadMov {ii}= HeadMov;
% SCentroid {ii}= Centroid;
% SHead_Vectors {ii} = Head_Vectors;
T{ii}=t;

% data filtered
% SAnglesF {ii}=AnglesF;
% SAngleFh {ii} =AnglesFh;
% SAngleFl {ii} =AnglesFl;
% SHead_VectorsF{ii}=Head_VectorsF;

Svariables {ii}= variables;

end

clear 'Angles' 'Nxyz' 'HeadMov' 'Centroid' 'Head_Vectors' 'AnglesFh' 'AnglesFl' 'Head_VectorsF' 't' 'xyzf' 'ii' 'AnglesF'
%% Calculate the Spectrum and Coherence
Sys = {};
Syc ={};
Sts={};
Sfs={};
Stc={};
Sfc={};
Sycl={};
Sysl={};

SyslnW={};
SysnW={};
for c= 1: size(Svariables,2);
    
    sig = Svariables{c};


[sigW, ARmodel] = WhitenSignal(sig, ceil(500*Trial.V_sr{c}),0,[], [], []);
padding = 2; 
nw = 3; %  Can be 1.5, 2, 2.5 or 3

% Important parameters
windowsize = 2.5; % seconds
overlap = 90; % percent
freqrange = [.5 20];
% NOTE frequency is top 60 becauase is the maximun that can be resolved
% with vicon sampling rate. the chosed window size is to observe the
% sniffing pattterns (10Hz, 100msec,minwindow is .2) and shaking that
% apears to be around 25Hz (40msec,minwindow is .08) 

Ntapes =3;
% Boring calculations
WinLength = 2^nextpow2(round(windowsize*Trial.V_sr{c}));
nFFT=2^(nextpow2(windowsize*Trial.V_sr{c}) + padding);
nOverlap   = ceil(WinLength*(overlap/100));



[ysl, fs, ts] = mtcsglong(sigW, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear', 2*nw-1, freqrange);
[yslnW, fs, ts] = mtcsglong(sig, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear', 2*nw-1, freqrange);

[ycl, fc, tc] = mtchglong(sig, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear', 2*nw-1, freqrange);

[ys, fs] = mtcsd(sigW, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear', 2*nw-1, freqrange);
[ysnW, fs] = mtcsd(sig, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear',2*nw-1, freqrange);

[yc, fc] = mtchd(sig, nFFT, Trial.V_sr{c}, WinLength, nOverlap, nw, 'linear', 2*nw-1, freqrange);

%Correct the timing
ts = ts + windowsize/2;
tc= tc + windowsize/2;

Sys{c}=ys;
Syc{c}=yc;
Sts{c}=ts;
Sfs{c}=fs;
Stc{c}=tc;
Sfc{c}=fc;
Sysl{c}=ysl;
SyslnW{c}=yslnW;
SysnW{c}=ysnW;
Sycl{c}=ycl;
end 

%% Periods of interest Vivon + Spectrum Errors (Spectrum Time)

GoodPR={};
SSpecErrors={};
for ii=1:length(SyslnW);
    
 % Define the time vector within the ranges of time in the spectrum
tmp = T{ii}-Sts{ii}(1);
[x,i]=min(abs(tmp));
in= i;    
tmp=[];
tmp = T{ii}-Sts{ii}(end);
[x,i]=min(abs(tmp));
en= i;  

tn=[];
tn= T{ii}>T{ii}(in) & T{ii}<T{ii}(en);

    
% Resampling the periods of interest
GoodPresampled = interp1(1:sum(tn),double(GoodP{ii}(tn)),linspace(1,sum(tn),length(SyslnW{ii})));
freqInt1 = Sfs{ii}>=.5 & Sfs{ii}<=16;

% Signal errors in the spectrum in the Pressure sensor (The threshold was set by visual inspection of the data, check session 5)
SpecErrors = (zscore(medfilt1(nanstd(tiedrank(SyslnW{ii}(:,~freqInt1,9)),[],2),15)))<-2.2;
% SpecErrors= zscore(medfilt1(nanstd(temporal(:,freqInt1),[],2),20))

SSpecErrors{ii}=~SpecErrors;
%Total Errors 
GoodPR {ii} = (GoodPresampled-double(SpecErrors)')>0;
% figure
% subplot(311)
% imagesc(Sts{ii},Sfs{ii},10*log10(SyslnW{ii}(:,:,9))');axis xy
% subplot(312)
% imagesc(Sts{ii},[],SpecErrors')
% subplot(313)
% imagesc(Sts{ii},[],((GoodPresampled-double(SpecErrors)')>0))
% axis xy
% link_axes


end
%% Excluding only  Spectrum errors (For mean PS RHM Coh All Sessions)

% StsG={};
% StcG={};
% 
% SyslG={};
% SyslnWG={};
% SyclG={};
% for ii=1:length(Sysl);
%     
% StsG{ii}=Sts{ii}(SSpecErrors{ii});
% StcG{ii}=Stc{ii}(SSpecErrors{ii});
% SyslG{ii}=Sysl{ii}(SSpecErrors{ii},:,:);
% SyslnWG{ii}=SyslnW{ii}(SSpecErrors{ii},:,:);
% SyclG{ii}=Sycl{ii}(SSpecErrors{ii},:,:,:);
% end
%% Spectrogram variables with good periods of behavior and without errors

StsG={};
StcG={};

SyslG={};
SyslnWG={};
SyclG={};
for ii=1:length(Sysl);
    
StsG{ii}=Sts{ii}(GoodPR{ii});
StcG{ii}=Stc{ii}(GoodPR{ii});
SyslG{ii}=Sysl{ii}(GoodPR{ii},:,:);
SyslnWG{ii}=SyslnW{ii}(GoodPR{ii},:,:);
SyclG{ii}=Sycl{ii}(GoodPR{ii},:,:,:);
end


% Inspection after cleaning

% figure
% for ii=1:16;
% subplotfit(ii,16);
% imagesc(StsG{ii},Sfs{ii},log10(SyslnWG{ii}(:,:,9))');
% axis xy
% end
%% Extract times of Interest Of Spectrum/Coherogram (Coherence)
Treshold=.6; % Treshold(For the figure of phase coupling the threshold is .6) for all session =

SgRHM= {};

for ii = 1: size(SyslG,2);
 
% temporal = (SyslG{ii}(:,:,8));% RHM
% temporal = (SyclG{ii}(:,:,8,9));% PS
freqInt1 = Sfs{ii}>=7 & Sfs{ii}<=12;
Specstd = (medfilt1(nanmedian(smoothn(SyclG{ii}(:,freqInt1,8,9),5),2),30));

% figure
% subplot(211)
% imagesc(StsG{ii},Sfs{ii}(freqInt1),(SyclG{ii}(:,freqInt1,8,9))')
% axis xy
% subplot(212)
% imagesc(StsG{ii},[],(Specstd>Treshold)')
% link_axes

% Threshold RHM high pow
SgRHM{ii} = Specstd>Treshold;
end
%% Extract times of Interest Of Spectrum/Coherogram (Ps power)
Treshold=-.7; % Treshold(For the figure of phase coupling the threshold is .6)

SgPspow= {};

for ii = 1: size(SyslG,2);
 
% temporal = (SyslG{ii}(:,:,8));% RHM
% temporal = (SyclG{ii}(:,:,8,9));% PS
freqInt1 = Sfs{ii}>=7 & Sfs{ii}<=12;
Specstd = zscore(medfilt1(nanmedian(log10(SyslnWG{ii}(:,freqInt1,9)),2),5));

% figure
% subplot(211)
% imagesc(StsG{ii},Sfs{ii}(freqInt1),log10(SyslnWG{ii}(:,freqInt1,9))')
% axis xy
% subplot(212)
% imagesc(StsG{ii},[],(Specstd>Treshold)')
% link_axes

% Threshold RHM high pow
SgPspow{ii} = Specstd>Treshold;
end


%% Spectrum densities for RHM periods 

FreqBinNum=length(Sfs{1}); % number of frequency bins of most of the sessions
Schsd={};
Spsd={};
SpsdnW={};


for ii=1:size(Sysl,2);
    tempPSD=[];
    tempPSDnW=[];
    tempfs=[];
    tempCHSD=[];
    
    tempPSD= sq(sum(SyslG{ii}(SgRHM{ii},:,:),1));
    tempPSDnW=sq(sum(SyslnWG{ii}(SgRHM{ii},:,:),1));
    tempCHSD=sq(sum(SyclG{ii}(SgRHM{ii},:,:,:),1)/sum(SgRHM{ii}));
    tempfs=Sfs{ii};
    
    % to correct the number of points due to difference in Vicon sampling rate interp1
    % is use to resample
    if length(tempPSD)>FreqBinNum
        tempresam=[];
        tempresamnW=[];
        tempresamCoh=[];
        for var=1:size(Sysl{1},3);
       tempresam(:,var)=interp1(1:length(tempPSD(:,var)),tempPSD(:,var)',linspace(1,length(tempPSD(:,var)),FreqBinNum));
       tempresamnW(:,var)=interp1(1:length(tempPSDnW(:,var)),tempPSDnW(:,var)',linspace(1,length(tempPSDnW(:,var)),FreqBinNum));
               for var2=1:size(Sysl{1},3);
       tempresamCoh(:,var,var2)=interp1(1:length(tempCHSD(:,var,var2)),tempCHSD(:,var,var2)',linspace(1,length(tempCHSD(:,var,var2)),FreqBinNum));
               end       
        end
        tempPSD=tempresam;
        tempPSDnW=tempresamnW;
        tempCHSD=tempresamCoh;
    end
    
    Spsd{ii}=tempPSD;
    SpsdnW{ii}=tempPSDnW;
    Schsd{ii}=tempCHSD;
end
%% Inspection of the signals Spectrum time of Interest (Coherence)


SyslRhmPow={};
SyslnWRhmPow={};
SyclRhmPow={};
StsRhm={};

for ii=1:size(Sysl,2);

    SyslRhmPow{ii}= SyslG{ii}(SgRHM{ii},:,:);
    SyslnWRhmPow{ii}=SyslnWG{ii}(SgRHM{ii},:,:);
    SyclRhmPow{ii}=SyclG{ii}(SgRHM{ii},:,:,:);
    StsRhm{ii}=StsG{ii}(SgRHM{ii});
    
end


figure
for ii=1:16;
subplotfit(ii,16);
imagesc(StsRhm{ii},Sfs{ii}(freqInt1),(SyclRhmPow{ii}(:,freqInt1,9,8))');
axis xy
title( Trial.Session{ii})
end
%% %% Inspection of the signals Spectrum time of Interest (Ps Power)


SyslSgPspow={};
SyslnWSgPspow={};
SyclSgPspow={};
StsSgPspow={};

for ii=1:size(Sysl,2);

    SyslSgPspow{ii}= SyslG{ii}(SgPspow{ii},:,:);
    SyslnWSgPspow{ii}=SyslnWG{ii}(SgPspow{ii},:,:);
    SyclSgPspow{ii}=SyclG{ii}(SgPspow{ii},:,:,:);
    StsSgPspow{ii}=StsG{ii}(SgPspow{ii});
    
end
%% AllSession PSD all Data Explore

% Azimuth
AzimuthPSD= cell2mat(cellfun(@(x) x(:,1),Spsd,'UniformOutput',false));
AzimuthChSD= cell2mat(cellfun(@(x) x(:,1,9),Schsd,'UniformOutput',false));

% Pitch
PitchPSD= cell2mat(cellfun(@(x) x(:,2),Spsd,'UniformOutput',false));
PitchChSD= cell2mat(cellfun(@(x) x(:,2,9),Schsd,'UniformOutput',false));

% Roll
RollPSD= cell2mat(cellfun(@(x) x(:,3),Spsd,'UniformOutput',false));
RollChSD= cell2mat(cellfun(@(x) x(:,3,9),Schsd,'UniformOutput',false));

% RHM
RHMPSD= cell2mat(cellfun(@(x) x(:,8),Spsd,'UniformOutput',false));
RHMChSD= cell2mat(cellfun(@(x) x(:,8,9),Schsd,'UniformOutput',false));

% Ps
PsPSD= cell2mat(cellfun(@(x) x(:,9),SpsdnW,'UniformOutput',false));
% Pitch Azimuth
PitchAzimuthChSD= cell2mat(cellfun(@(x) x(:,1,2),Schsd,'UniformOutput',false));



figure
subplot(3,4,1)
imagesc(Sfs{1},[],zscore(10*log10(AzimuthPSD))')
axis xy
title 'Azimuth'
set(gca,'YTick',[1:length(Trial.xyz)])
set(gca,'YTickLabel',Trial.Session)
caxis ([-2 2])
colorbar

subplot(3,4,2)
imagesc(Sfs{1},[],((AzimuthChSD))')
axis xy
title 'Azimuth PsCoh'
caxis ([0 1])
colorbar

subplot(3,4,3)
imagesc(Sfs{1},[],zscore(10*log10(PitchPSD))')
axis xy
title 'Pitch'
caxis ([-2 2])
colorbar

subplot(3,4,4)
imagesc(Sfs{1},[],((PitchChSD))')
axis xy
title 'Pitch PsCoh'
caxis ([0 1])
colorbar

subplot(3,4,5)
imagesc(Sfs{1},[],zscore(10*log10(RollPSD))')
axis xy
title 'Roll'
set(gca,'YTick',[1:length(Trial.xyz)])
set(gca,'YTickLabel',Trial.Session)
caxis ([-2 2])
colorbar

subplot(3,4,6)
imagesc(Sfs{1},[],((RollChSD))')
axis xy
title 'Roll PsCoh'
caxis ([0 1])
colorbar

subplot(3,4,7)
imagesc(Sfs{1},[],zscore(10*log10(RHMPSD))')
axis xy
title 'RHM'
caxis ([-2 2])
colorbar

subplot(3,4,8)
imagesc(Sfs{1},[],((RHMChSD))')
axis xy
title 'RHM PsCoh'
caxis ([0 1])
colorbar

subplot(3,4,9)
imagesc(Sfs{1},[],zscore(10*log10(PsPSD))')
axis xy
title 'Ps'
caxis ([-2 2])
colorbar
set(gca,'YTick',[1:length(Trial.xyz)])
set(gca,'YTickLabel',Trial.Session)

subplot(3,4,10)
imagesc(Sfs{1},[],(PitchAzimuthChSD)')
axis xy
title 'Pitch Azimuth ChSD '
caxis ([0 1])
colorbar

%% Power correlations of PSD (Supplementary figure)


freqInterest = Sfs{1}>0 &Sfs{1}<20;
AzPitPowRatio = bsxfun(@rdivide, bsxfun(@minus,AzimuthPSD,PitchPSD),bsxfun(@plus, AzimuthPSD,PitchPSD));

figure
subplot(3,1,1)
plot(Sfs{1},(AzPitPowRatio))
hold on
errorbar(Sfs{1},mean(AzPitPowRatio,2),std(AzPitPowRatio,[],2),'y')
xlim ([0 20])
xlabel ('Frequency')
ylabel ('Power Ratio Azimuth-Pitch ')


subplot(3,1,2)
freqInterest = Sfs{1}>0 &Sfs{1}<13;

for ii=1: size(PitchPSD,2);
    scatter(tiedrank(10*log10(AzimuthPSD(freqInterest,ii))),tiedrank(10*log10(PitchPSD(freqInterest,ii))),20,Sfs{1}(freqInterest),'Filled')
    hold on
    colormap 'jet'
    colorbar
    title ('Power correlation of Azimuth vs Pitch')
    xlabel ('Azimuth Rank Power')
    ylabel ('Pitch Rank Power')
    
end    

subplot(3,1,3)

Srho=[];
Spval=[];
for ii=1: size(PitchPSD,2);
    
    [rho,pval]=corr((AzimuthPSD(freqInterest,ii)),(PitchPSD(freqInterest,ii)),'type','Spearman');
   Srho(ii)=rho;
   Spval(ii)=pval;
end  

boxplot(Srho)
hold on
plot(ones(1,16),Srho,'ok')
ylim ([0 1])
ylabel ('Correlation')
title('Spearman Correlation')

%% Spectra and coherence Angles and Ps (Example plot for figure)
sm =5; %smooth factor
figure
subplot(7,1,1)
imagesc(Sts{1},Sfs{1},smoothn(10*log10(Sysl{1}(:,:,1)'),sm));
axis xy
colorbar
title ('Azimuth')

subplot(7,1,2)
imagesc(Sts{1},Sfs{1},smoothn(10*log10(Sysl{1}(:,:,2)'),sm));
axis xy
colorbar
title ('Pitch')

subplot(7,1,3)
imagesc(Sts{1},Sfs{1},smoothn(10*log10(Sysl{1}(:,:,8)'),sm));
axis xy
colorbar
title ('RHM')


subplot(7,1,4)
imagesc(Sts{1},Sfs{1},smoothn(10*log10(SyslnW{1}(:,:,9)'),sm));
axis xy
colorbar
title ('Pressure sensor')


subplot(7,1,5)
imagesc(Sts{1},Sfs{1},smoothn((Sycl{1}(:,:,1,9)'),sm));
axis xy
colorbar
title ('Coh Ps Azimuth')


subplot(7,1,6)
imagesc(Sts{1},Sfs{1},smoothn((Sycl{1}(:,:,2,9)'),sm));
axis xy
colorbar
xlabel ('Time (sec)')
ylabel ('Frequency')

title ('Coh Ps Pitch')

subplot(7,1,7)
imagesc(Sts{1},Sfs{1},smoothn((Sycl{1}(:,:,7,9)'),sm));
axis xy
colorbar
xlabel ('Time (sec)')
ylabel ('Frequency')

title ('Coh Ps RHM')

link_axes
%% Summary PSD and CHSD for figure 

%PSD
figure
plot(Sfs{1},[mean(zscore(PitchPSD),2)+std(zscore(PitchPSD),[],2),mean(zscore(PitchPSD),2)-std(zscore(PitchPSD),[],2)],'b')
hold on
plot(Sfs{1},mean(zscore(PitchPSD),2),'b')

plot(Sfs{1},[mean(zscore(AzimuthPSD),2)+std(zscore(AzimuthPSD),[],2),mean(zscore(AzimuthPSD),2)-std(zscore(AzimuthPSD),[],2)],'r')
plot(Sfs{1},mean(zscore(AzimuthPSD),2),'r')

plot(Sfs{1},[mean(zscore(RollPSD),2)+std(zscore(RollPSD),[],2),mean(zscore(RollPSD),2)-std(zscore(RollPSD),[],2)],'g')
plot(Sfs{1},mean(zscore(RollPSD),2),'g')

plot(Sfs{1},[mean(zscore(PsPSD),2)+std(zscore(PsPSD),[],2),mean(zscore(PsPSD),2)-std(zscore(PsPSD),[],2)],'k')
plot(Sfs{1},mean(zscore(PsPSD),2),'k')
title ('PSD Angles')
xlabel ('Frequency')
ylabel ('Power (z-score)')

% CHSD
figure
plot(Sfs{1},[mean((PitchChSD),2)+std((PitchChSD),[],2),mean((PitchChSD),2)-std((PitchChSD),[],2)],'b')
hold on
plot(Sfs{1},mean((PitchChSD),2),'b')

plot(Sfs{1},[mean((AzimuthChSD),2)+std((AzimuthChSD),[],2),mean((AzimuthChSD),2)-std((AzimuthChSD),[],2)],'r')
plot(Sfs{1},mean((AzimuthChSD),2),'r')

plot(Sfs{1},[mean((RollChSD),2)+std((RollChSD),[],2),mean((RollChSD),2)-std((RollChSD),[],2)],'g')
plot(Sfs{1},mean((RollChSD),2),'g')

plot(Sfs{1},[mean((RHMChSD),2)+std((RHMChSD),[],2),mean((RHMChSD),2)-std((RHMChSD),[],2)],'k')
plot(Sfs{1},mean((RHMChSD),2),'k')


title ('ChSD Angles')
xlabel ('Frequency')
ylabel ('Coherence')

%% Extracting Period of interest Spectrum time to Normal time (Coherence) 
%!!! Note that from now on the time vector is different !!!
SclipTime={};
SNGoodPR={};
SNSgRHM={};
SNgPspow={};
RhmT={};
SvariablesRhm={};
for ii=1:length(Trial.Session);

% finding the begin and end of the signal corresponding to the spectrogram

tmp = T{ii}-Sts{ii}(1);% note that we are using the time after error subtraction(GoodPR Vicon+Spectrum errors)
[x,ind1]=min(abs(tmp));
in= ind1;    
tmp=[];
tmp = T{ii}-Sts{ii}(end);
[x,ind2]=min(abs(tmp));
en= ind2;  

% Define the time vector within the above ranges to match spectrum time
clipTime= T{ii}>T{ii}(ind1) & T{ii}<T{ii}(ind2);
tmpVarSpec_time=[];
% Variables within spectrum time
tmpVarSpec_time = Svariables{ii}(clipTime,:);
tempTime1=T{ii}(clipTime);
% interpolate the vector of periods of Vicon + Spectrum errors
tmp=[];
tmp = interp1(1:length(GoodPR{ii}),double(GoodPR{ii}),linspace(1,length(GoodPR{ii}), length(T{ii}(clipTime))));
tmp(isnan(tmp))=0;
NGoodPR= logical(tmp);

tmpVarClean=tmpVarSpec_time(NGoodPR,:);
tempTime2=tempTime1(NGoodPR);

% Interpolate the RHM periods
tmp=[];
tmp = interp1(1:length(SgRHM{ii}),double(SgRHM{ii}),linspace(1,length((SgRHM{ii})), length(tmpVarClean)));
NSgRHM=logical(tmp);
 %Interpolate Ps periods
tmp=[];
tmp = interp1(1:length(SgPspow{ii}),double(SgPspow{ii}),linspace(1,length((SgPspow{ii})), length(tmpVarClean)));
NSgPspow=logical(tmp);



SvariablesRhm{ii}=tmpVarClean(NSgRHM,:);
RhmT{ii}=tempTime2(NSgRHM);
SclipTime{ii}=clipTime;
SNGoodPR{ii}=NGoodPR;
SNSgRHM{ii}=NSgRHM;
SNgPspow{ii}=NSgPspow;
end
%% Variables for Phase analysis


SBackV = cellfun(@(x) x(:,7), SvariablesRhm,'UniformOutput',false);
SRhm =cellfun(@(x) x(:,8), SvariablesRhm,'UniformOutput',false);
SPs =cellfun(@(x) x(:,9), SvariablesRhm,'UniformOutput',false);
%% Find Peaks pairs RHM PS 

PairsVec={};
PairsRhm={};

ValuesVec={};
ValuesRhm={};

StepsVec={};
StepsRhm={};

IndPairsVec={};
IndPairsRhm={};

PeaksVec={};
PeaksRhm={};

TroughsVec={};
TroughsRhm={};





for ii = 1:size(SBackV,2); 
    
VectF= ButFilter(SBackV{ii},[2],[6/(Trial.V_sr{ii}/2) 20/(Trial.V_sr{ii}/2)],'bandpass');
RhmF = ButFilter(SRhm{ii},[2],[6/(Trial.V_sr{ii}/2) 20/(Trial.V_sr{ii}/2)],'bandpass');

% Back vector    
TroughVec = LocalMinima(VectF, Trial.V_sr{ii}*0.05,0); 
PeakVec = LocalMinima(-VectF, Trial.V_sr{ii}*0.05,0);
% RHM
TroughRhm = LocalMinima(RhmF, Trial.V_sr{ii}*0.05,0); 
PeakRhm = LocalMinima(-RhmF, Trial.V_sr{ii}*0.05,0);


% find matching pairs

[xv, xvi ,yvi] = NearestNeighbour(PeakVec,TroughVec,'right', Trial.V_sr{ii}*0.1);
[xr, xri ,yri] = NearestNeighbour(PeakRhm,TroughRhm,'right', Trial.V_sr{ii}*0.1);

% Store the corresponding indexes

IndPairVec =[xvi yvi];
PairVec = [PeakVec(xvi) TroughVec(yvi)];

IndPairRhm =[xri yri];
PairRhm = [PeakRhm(xri) TroughRhm(yri)];

% Get the values

ValVec = zeros(size(PairVec));
ValRhm = zeros(size(PairRhm));


ValVec(:,1) = VectF(PairVec(:,1));
ValVec(:,2) = VectF(PairVec(:,2));

ValRhm(:,1) = RhmF(PairRhm(:,1));
ValRhm(:,2) = RhmF(PairRhm(:,2));

% Get the step

StepVec = diff(ValVec,1,2);
StepRhm = diff(ValRhm,1,2);

% GoodMv = Step>1;
% 
% Pair = Pair(GoodMv,:);
% Val =Val(GoodMv,:);
% Step = Step(GoodMv,:);

%Store all the data
PairsVec{ii}=PairVec;
PairsRhm{ii}=PairRhm;

ValuesVec{ii}=ValVec;
ValuesRhm{ii}=ValRhm;

StepsVec{ii}=StepVec;
StepsRhm{ii}=StepRhm;

IndPairsVec{ii}=IndPairVec;
IndPairsRhm{ii}=IndPairRhm;

PeaksVec{ii}=PeakVec;
PeaksRhm{ii}=PeakRhm;

TroughsVec{ii}=TroughVec;
TroughsRhm{ii}=TroughRhm;


end
%% Time Shift Parameters for Raylagh statistics

bs=2;
n=21;

d1 = log(1)/log(bs);
d2 = log(25)/log(bs);
y = (bs).^ [d1+(0:n-2)*(d2-d1)/(floor(n)-1), d2];
FreqRange=[];
FreqRange(:,1) = ((y(1:end-1)+y(2:end))/2)';
FreqRange(:,2) = diff(y(:))*2;
V_nq = Trial.V_sr{ii}/2;

Shift=zeros(length(Trial.xyz),n);
for ii=1:length(Trial.xyz);
timebin=.05;%in seconds (.50 it correspond to the time evaluated in rats with V_sr=199.9770)
Nsamples=round(timebin*Trial.V_sr{ii});
Shift(ii,:) =-Nsamples*((n-1)/2):Nsamples:Nsamples*((n-1)/2);
end
%% Time ShiftBack Vector Statistics
out=[];
RayStatsVec={};
for c = 1:size(SBackV,2);
    
DnTime = reshape(PairsVec{c}',[],1);% create a single vector with the Pairs of detected motion
DnClu = repmat([1 2]', size(PairsVec{c},1),1); % crate a vector with the identity of data (peak or trough)



for ii=1:size(Shift,2)
    for jj=1:size(FreqRange,1)
        
        myFreqRange = FreqRange(jj,1)+[-0.5 0.5]*FreqRange(jj,2);
        fps = ButFilter(SPs{c},2,myFreqRange/(Trial.V_sr{c}/2),'bandpass');
        hlb = hilbert(fps);
        ph = angle(hlb);
        DnTimeShifted= DnTime+Shift(c,ii); 
        DnTimeGood = DnTimeShifted>0 & DnTimeShifted<length(SPs{c});
        rt = RayleighTest(ph(DnTime(DnTimeGood)+Shift(c,ii)),DnClu(DnTimeGood),2);
        
        out.r(ii,jj,:) = rt.r;
        out.th0(ii,jj,:) = rt.th0;
        
    end
end
RayStatsVec{c}=out;

end
%% Time Shift Rhm Statistics
out=[];
RayStatsRhm={};
for c = 1:size(SRhm,2);
    
DnTime = reshape(PairsRhm{c}',[],1);
DnClu = repmat([1 2]', size(PairsRhm{c},1),1);



for ii=1:size(Shift,2)
    for jj=1:size(FreqRange,1)
        
        myFreqRange = FreqRange(jj,1)+[-0.5 0.5]*FreqRange(jj,2);
        fps = ButFilter(SPs{c},2,myFreqRange/(Trial.V_sr{c}/2),'bandpass');
        hlb = hilbert(fps);
        ph = angle(hlb);
        DnTimeShifted= DnTime+Shift(c,ii); 
        DnTimeGood = DnTimeShifted>0 & DnTimeShifted<length(SPs{c});
        rt = RayleighTest(ph(DnTime(DnTimeGood)+Shift(c,ii)),DnClu(DnTimeGood),2);
        
        out.r(ii,jj,:) = rt.r;
        out.th0(ii,jj,:) = rt.th0;
        
    end
end
RayStatsRhm{c}=out;

end
%% Time Shift Resume of statistics 

evt = 1; % 1= detected peak of head motion 2= trough of head motion
signal=[];
signal = RayStatsRhm; % chosse either frontal vector or RHM

temp=cell2mat(signal);
TRayRstats =zeros(n,n-1,length(temp));
for ii = 1:length(temp);
    
    TRayRstats(:,:,ii)= temp(ii).r(:,:,evt);
end    
%% Time Shift All Sessions Phase coupling  

figure

for ii=1:16;
subplotfit(ii,16);
pcolor(Shift(ii,:)/Trial.V_sr{ii}*1000, FreqRange(:,1), TRayRstats(:,:,ii)');
shading flat
colorbar
% caxis ([0 .6])
xlabel ('Time shift (msec)')
ylabel ('Frequency')
title(Trial.Session{ii})
end
%% Time Shift maximun of coupling strength

timemax=zeros(n-1,size(TRayRstats,3));
freqmax=zeros(n,size(TRayRstats,3));
Smrow=zeros(1,size(TRayRstats,3));
Smcol=zeros(1,size(TRayRstats,3));
for ii=1:size(TRayRstats,3);
tmp=[];    
tmp = sq(TRayRstats(:,:,ii));
[v,in]=max(tmp(:));% Here the matrix is linearized to find the max

[mrow,mcol]=ind2sub(size(tmp),in);% Extracting indexes in matrix form

timemax(:,ii)= tmp(mrow,:);
freqmax(:,ii)= tmp(:,mcol);
Smrow(ii)=mrow;
Smcol(ii)=mcol;
end
% test = round(linspace(1,size(FreqRange(:,1),1),11));
%% Time Shift Summary Phase couplig strength Plot
Ses=1:16;
ExcluSes=(Ses == 13); 

figure
h = subplot(2,2,1);
pcolor(Shift(ii,:)/Trial.V_sr{ii}*1000, FreqRange(:,1), TRayRstats(:,:,ii)');
shading flat
colorbar
xlabel ('Time shift (msec)')
ylabel ('Frequency')
title('Phase coupling strength ')
% hold on
% plot(Shift(Smrow)/Trial.V_sr{ii}*1000,FreqRange(Smcol,1),'.k','MarkerSize',4)


h2= subplot(2,2,2);
plot(FreqRange(:,1),median(timemax(:,~ExcluSes),2),'k', 'LineWidth', 2)
hold on
plot(FreqRange(Smcol(~ExcluSes),1),linspace(.001,.1,length(Smcol(~ExcluSes))),'*k', 'MarkerSize', 6)
plot(FreqRange(:,1),median(timemax(:,~ExcluSes),2)+std(timemax(:,~ExcluSes),[],2),'k')
plot(FreqRange(:,1),median(timemax(:,~ExcluSes),2)-std(timemax(:,~ExcluSes),[],2),'k')
hold off
xlim ([1 23])
ylim ([0 .9 ])
view ([90 -90])


h3 = subplot(2,2,3);
plot(Shift(ii,:)/Trial.V_sr{ii}*1000,mean(freqmax(:,~ExcluSes),2),'k', 'LineWidth', 2)
hold on
% for pp=1:16;
% plot(Shift(pp,Smrow(pp))/Trial.V_sr{pp}*1000,.01*pp,'*k', 'MarkerSize', 6)
% end

plot(Shift(ii,Smrow(~ExcluSes))/Trial.V_sr{ii}*1000,linspace(.001,.1,length(Smrow(~ExcluSes))),'*k', 'MarkerSize', 6)

plot(Shift(ii,:)/Trial.V_sr{ii}*1000,mean(freqmax(:,~ExcluSes),2)+std(freqmax(:,1:16),[],2),'k')
plot(Shift(ii,:)/Trial.V_sr{ii}*1000,mean(freqmax(:,~ExcluSes),2)-std(freqmax(:,1:16),[],2),'k')

xlim ([-500 500 ])

set(h2,'position',[.78 .05 .15 .7]);
set(h3,'position',[.05 .8 .6 .15]);
set(h,'position',[.05 .05 .6 .7]);


%% Respiration Phase Distribution of frontal marker motion 

myFreqRange = [7 12];
% finding the phase values for a given peak and trough

SPhase={};
for ii=1:4;

fps = ButFilter(SPs{ii},2,myFreqRange/V_nq,'bandpass'); 
hlb = hilbert(fps);
ph = angle(hlb);
SPhase{ii}=[ph(PairsVec{ii})];
end

% concatenate all data
AllPhase = cat(1,SPhase{:});


% % overlay scatter
% figure
% for ii=1:size(SPhase,2);
% scatter(SPhase{ii}(:,1),SPhase{ii}(:,2),'.')
% hold on
% end

% find pressure sensor peaks to plot the mean shape of pressure sensor

fps = ButFilter(SPs{1},2,myFreqRange/V_nq,'bandpass'); 
[v,i] = findpeaks(fps,'MinPeakProminence',1000);

% to avoid the borders
i=i(i>200 & i<i(end)-200);

% boundaries of Ps cycle
limt = 7;
id=[i-limt i+limt];

% vector to plot over the phase
tps = linspace(-3,3, limt*2+1);

% collect all data 

% Ps
psP=[];
for ii=1:size(id,1);
psP(ii,:) = fps(id(ii,1):id(ii,2));
end

% RHM
rhmP = [];
for ii=1:size(id,1);
rhmP(ii,:) = SBackV{1}(id(ii,1):id(ii,2));
end


% overlay note that the Ps signal is scaled to the ratio beetween the two
%signals

% in case that we need to plot the probability
OutT=[];
for ii= 1:size(SPhase,2);
    
[Out, XBins, YBins, Pos] = hist2(SPhase{ii}, 50, 50);

%
OutT(:,:,ii)= (Out/(size(SPhase{ii},1)));
end
  


% ploting
figure
h1=subplot(2,1,1);
imagesc(XBins,YBins,mean(OutT,3)');
axis xy
xlabel ('Peak Head motion (Respiration Phase)')
ylabel ('Trough Head motion (Respiration Phase)')
title ('Mean JPD Head Motion and Respiration Phase')
colorbar

hold on
%ps
plot(tps,mean(psP)/max(mean(psP)/max(ph)),'w','linewidth',2)
% rhm
plot(tps,mean(rhmP)/max(mean(rhmP)/max(ph)),'r','linewidth',2)
hold off

% countours for to overlay JPD
h2=subplot(2,1,2);
tmp=mean(OutT,3);

contour(smoothn(tmp,1)',[prctile(tmp(:),[ 70 80 95]) prctile(tmp(:),[ 70 80 95])]);

Lev90=smoothn(tmp,1)>prctile(tmp(:),95);
Lev80=smoothn(tmp,1)>prctile(tmp(:),80);
Lev70=smoothn(tmp,1)>prctile(tmp(:),70);

PrbLev90 = sum(tmp(Lev90));
PrbLev80 = sum(tmp(logical(Lev80-Lev90)));
PrbLev70 = sum(tmp(logical(Lev70-Lev80)));
ylabel ([PrbLev90; PrbLev80; PrbLev70])
set(h1,'position',[.1 .5 .8 .4]);
set(h2,'position',[.1 .05 .8 .4]);
%% RHM Ps Power relationship
SgCoh= {};

for ii = 1: size(Sysl,2);
    
    
tmp = (Sycl{ii}(:,:,8,9));
maxC = max(max(tmp));
tmp=tmp/maxC;

freqInt1 = Sfs{ii}>=6 & Sfs{ii}<=11;


mtestC = zscore(smoothn(mean(tmp(:,freqInt1),2),5));

SgCoh{ii} = mtestC>-.6;
end



%
SgPs= {};

for ii = 1: size(Sysl,2);
    
    
tmp = (SyslnW{ii}(:,:,9));
maxS = max(max(10*log10(tmp)));
tmp = tmp/maxS;

freqInt1 = Sfs{ii}>=6 & Sfs{ii}<=11;

mtestS = zscore(smoothn(std(tmp(:,freqInt1),[],2),5));
% mtestS1 =zscore(smoothn(std(tmp(:,freqInt1),[],2),5));
SgPs{ii} = mtestS>-1;
end
%% Ploting to check Coherence and PsPow extraction
S =3;
maxS = max(max(10*log10(SyslnW{S}(:,:,9))));
maxC = max(max((Sycl{S}(:,:,8,9))));

figure
subplot(4,1,1)
imagesc(Stc{S},Sfs{S},(Sycl{S}(:,:,8,9)/maxC)')
axis xy
subplot(4,1,2)
imagesc(Stc{S},[],SgCoh{S}')
link_axes

subplot(4,1,3)
imagesc(Sts{S},Sfs{S},(10*log10(SyslnW{S}(:,:,9))/maxS)')
axis xy
subplot(4,1,4)
imagesc(Sts{S},[],SgPs{S}')
axis xy
link_axes
%% Time Proportion of RHM Ps Coherence 

 SpecTSum(:,1)= cellfun(@sum,SgRHM,'UniformOutput',false ); % RHM
 SpecTSum(:,2)= cellfun(@sum,SgPs,'UniformOutput',false ); % Ps
 SpecTSum(:,3)= cellfun(@sum,SgCoh,'UniformOutput',false ); % Coherence

SpecTSum = cell2mat(SpecTSum);
 
PrhmPS = (SpecTSum(:,1)*100)./SpecTSum(:,2); % proportion of time of RHM respect to Ps
PcohRhmPs =(SpecTSum(:,3)*100)./SpecTSum(:,2); % proportion of coherence respect to Ps

% 
Resume=[];
for ii=1:size(SgRHM,2);
A = SgPs{ii}+SgRHM{ii};
B = A+SgCoh{ii};
% C = A+B;

[v,in]=hist(double(B),unique(B)');

Resume(:,:,ii)=[v;in];
end

