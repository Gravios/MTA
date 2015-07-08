
function CSDTrigAverTheta(Session, Channels, ThChannel, varargin)
%CSDTrigAverTheta is a function which calculates theta trough-triggered averaged LFP traces 
%and corresponding CSD map. 
%
%USAGE:     CSDTrigAverTheta(Session, Channels, ThChannel, <SignalType>, <PeriodTitle>, <TrigWinDur_ms>, <CSDstep>)
%
%INPUT:
% Session        is a name of the session (AnimalID-YYYMMDD).
% Channels        is a list of channels for which CSD must be computed.
% ThChannel       is the channel on which theta must be detected
% <SignalType>    is a type of signal used ('lfp' or 'interplfp'). Default='lfp.'
% <PeriodTitle>   is a name of the periods  from which data must be taken. Default='ALL'.
% <TrigWinDur_ms> is a duration (ms) of the window for triggered averaging. Default=375 ms
% <CSDstep>       is a spatial step (channels) of CSD. Use  CSDstep=2 for probes with 50-100 mkm ntersite distance; 
%                 and CSDstep=8 for intersite distance <50mkm. Default=2.
%
%OUTPUT FILES:
% Session.CSDTrigAverTheta.PeriodTitle.jpg  
%
%EXAMPLE:     CSDTrigAverTheta(Session, 65:96, 72, 'lfpinterp', 'RUNTHETA')
%
% Evgeny Resnik
% version 02.08.2012




%Check number of parameters
if nargin < 3
    error('USAGE:   CSDTrigAverTheta(Session, Channels, ThChannel, <SignalType>, <PeriodTitle>, <TrigWinDur_ms>, <CSDstep>)');
end

% Parse input parameters
[SignalType, PeriodTitle, TrigWinDur_ms, CSDstep ] = DefaultArgs(varargin,{ 'lfp', 'ALL', 250, 2 });


%--------------------------------------------------------------------------------------%


fprintf('=======================================================================================================\n')
fprintf(['                   %s - %s \n'], Session.name, mfilename )
fprintf('=======================================================================================================\n')


%Load spike and LFP sampling rates from .xml file
par = LoadXml(fullfile(Session.spath,[Session.name '.xml']));
lfpSamplingRate = par.lfpSampleRate;
nChan = par.nChannels;


%Load start/end timestamps (.lfp samples at 1250Hz) of episodes that must be included/excluded
if ~strcmp(PeriodTitle,'ALL')
    stsFilesIn = fullfile(Session.spath,[Session.name '.sts.' PeriodTitle]);

     if exist(stsFilesIn,'file'),
%    if cellfun(@prod,cellfun(@exist,cellfun(@fullfile,repmat({Session.spath},numel(stsFilesIn),1),stsFilesIn,'UniformOutput',false),repmat({'file'},numel(stsFilesIn),1),'UniformOutput',false)),
        Periods = loadrangefiles({stsFilesIn});
    elseif sum(~cellfun('isempty',{Session.stc{PeriodTitle}}))==1,
        Periods = Session.stc{PeriodTitle,Session.lfp.sampleRate}; 
    end
else
    Periods = [];
end

%Load LFP data 
fprintf(['Loading %s from periods %s  ...'], upper(SignalType), PeriodTitle)
Session.lfp.load(Session,Channels);%,[],Periods); 
lfp = Session.lfp.data;
%LoadBinary([Session.name '.' SignalType], Channels, nChan,[],[],[], Periods);
fprintf('DONE\n')


%Load separately LFP from theta channel (may differ from Channels)
fprintf(['Loading %s on theta channel-%d from periods %s  ...'], upper(SignalType), ThChannel, PeriodTitle)
Session.lfp.clear;
Session.lfp.load(Session,ThChannel);%,[],Periods); 
thlfp = Session.lfp.data;
%thlfp = LoadBinary([Session.name '.' SignalType], ThChannel, nChan,[],[],[], Periods);
fprintf('DONE\n')


%--------------------------------------------------------------------------------------%

%Filter LFP in the theta range (5-12 Hz) to remove low-frequency fluctuations
%and high-frequency events which can biase CSD
FreqRange = [5 12];
fprintf(['Filtering LFP in theta range (%1.0f-%1.0f Hz) ...'], FreqRange(1), FreqRange(2))
Nyquist = 0.5*lfpSamplingRate;
[n, wn] = buttord([5 12]/Nyquist, [3 15]/Nyquist, 3, 13);
flfp = ButFilter(lfp, n, wn, 'bandpass');
fthlfp = ButFilter(thlfp, n, wn, 'bandpass');
clear Nyquist n wn samples
fprintf('DONE\n')


fthlfp = MTADlfp('data',fthlfp,'sampleRate',Session.lfp.sampleRate);
fthlfp = fthlfp(Periods,:);

flfp = MTADlfp('data',flfp,'sampleRate',Session.lfp.sampleRate);
flfp = flfp(Periods,:);


%Detect theta troughs on the specified channel of the filtered LFP
fprintf(['Detecting theta troughs in filtered LFP ...'])
% ind = find(Channels==ThChannel);
% ThTroughs = negativepeak(flfp(ind,:));
% clear ind
ThTroughs = negativepeak(fthlfp);
fprintf('DONE\n')


%Calculate triggered averaged LFP traces
fprintf(['Calculating theta trough-triggered averaged traces ...'])
HalfTrigWinDur = round(0.5*TrigWinDur_ms/10^3*lfpSamplingRate);
TrigAverLFP = TriggeredAv(flfp, HalfTrigWinDur, HalfTrigWinDur, ThTroughs');
fprintf('DONE\n')


%Compute CSD
TrigAverCSD = CurSrcDns(TrigAverLFP, [-HalfTrigWinDur:HalfTrigWinDur]/lfpSamplingRate*10^3, [],[],[],[],CSDstep);
oTrigAverCSD = CurSrcDns(flfp,[], [],[],[],[],CSDstep);
oThTroughs = negativepeak(oTrigAverCSD(:,30));
TrigAverLFP = TriggeredAv(oTrigAverCSD, HalfTrigWinDur*8, HalfTrigWinDur*8, oThTroughs');
% $$$ imagesc(TrigAverLFP')
% $$$ 
% $$$ keyboard



%create a time vector (ms)
t = [-HalfTrigWinDur:HalfTrigWinDur]/lfpSamplingRate*10^3;


%Save data into a file
% $$$ FileOut = sprintf(['%s.%s.%s.%s.%d-%d.mat'], Session.name, mfilename, SignalType , PeriodTitle, Channels([1 end]) );
% $$$ fprintf(['Saving data into a file %s ...'], FileOut )
% $$$ Params = struct('SignalType', SignalType, 'PeriodTitle',PeriodTitle, 'ThChan', ThChannel, 'Channels', Channels, 'lfpSamplingRate', lfpSamplingRate, 'TrigWinDur_ms', TrigWinDur_ms);
% $$$ save(FileOut ,'t', 'TrigAverLFP', 'TrigAverCSD', 'Params','-v7.3');
% $$$ fprintf('DONE\n')



%--------------------------------------------------------------------------------------%
%Prepare some parameters for plotting

%font size
FontSize=12;

%scale traces
scale = 2;
% $$$ shiftedTrigAverLFP = TrigAverLFP*scale;
% $$$ %shift traces by the mean of the first two and last two samples
% $$$ shiftedTrigAverLFP = shiftedTrigAverLFP - repmat(mean(TrigAverLFP([1:2 end-1:end],:),1),size(TrigAverLFP,1),1);
% $$$ shift = max(max(TrigAverLFP)-min(TrigAverLFP));
% $$$ shiftedTrigAverLFP = shiftedTrigAverLFP - shift/2-repmat([0:length(Channels)-1]*shift, size(TrigAverLFP,1),1);
%create Y ticks values
% $$$ yticks = sort(shiftedTrigAverLFP(1,1:1:end));
% $$$ yticklabels =  flipud(Channels(:));
%yticks for CSD
newchrange = linspace(Channels(1)+2-0.5, Channels(end)-2+0.5, size(TrigAverCSD,2))';
%Y axis increment for plotting traces
nCh = size(TrigAverLFP,2);
AxesIncr = 2*2/(nCh-2*2);

%Plot figure
figure; orient landscape
%plot triggered averaged traces
subplot(1,2,1); cla; hold on
% $$$ plot(t, shiftedTrigAverLFP,'color','k','linestyle','-','linewidth',1)
% $$$ if ismember(ThChannel, Channels)
% $$$     plot(t, shiftedTrigAverLFP(:,find(Channels==ThChannel)),'color','k','linestyle','-','linewidth',2)
% $$$ end
% $$$ plot([0 0],[-shift*length(Channels) shift*0.3],'color','k','linestyle',':','linewidth',1)
% $$$ set(gca,'YLim',[-shift*length(Channels) shift*0.3],'YTick', yticks,'YTicklabel', yticklabels,'FontSize',FontSize);
% $$$ set(gca,'XLim',t([1 end]),'FontSize',FontSize)
imagesc(t, newchrange, TrigAverLFP'); axis xy,axis tight
xlabel('Time, (ms)','FontSize',FontSize)
ylabel('Channels','FontSize',FontSize)
title('Triggered averaged trace','FontSize',FontSize)
set(gca,'YTick', Channels(3:end-2),'YTicklabel', Channels(3:end-2) , 'YDir','reverse', 'FontSize',FontSize);
set(gca,'XLim',t([1 end]),'FontSize',FontSize)

%plot triggered averaged CSD
h=subplot(1,2,2); cla; hold on
xlabel('Time, (ms)','FontSize',FontSize)
ylabel('Channels','FontSize',FontSize)
title('Triggered averaged CSD','FontSize',FontSize)
imagesc(t, newchrange, TrigAverCSD'); axis tight
set(gca,'YTick', Channels(3:end-2),'YTicklabel', Channels(3:end-2) , 'YDir','reverse', 'FontSize',FontSize);
set(gca,'XLim',t([1 end]),'FontSize',FontSize)
%scale colors
cx = caxis;
cxmax = max(abs(cx));
caxis([-cxmax cxmax]);
%traces
% $$$ PlotManyCh(TrigAverLFP, t, 1 , 'k', 0, AxesIncr);

%---------------------------- Suptitle ---------------------------------------------%
titlestr{1} = ['Theta trough-triggered averaged CSD: ' Session.name  ', ' SignalType  ', ' PeriodTitle ', ThChan=' num2str(ThChannel)   ];
titlestr{2} = sprintf('(lfpSamplingRate=%0.0f Hz, TrigWinDur_ms=%0.0f)', lfpSamplingRate,  TrigWinDur_ms);
suptitle2(titlestr, .98, 1)
clear titlestr





%Save figure into a file
% $$$ FileOut = sprintf(['%s.%s.%s.%s.%d-%d'], Session.name, mfilename, SignalType , PeriodTitle, Channels([1 end]) );
% $$$ fprintf(['Saving figure  into a file %s ...'], FileOut)
% $$$ AxisPropForIllustrator
% $$$ FigPathCommon = [pwd '/figures/'];
% $$$ mkdir2([ FigPathCommon mfilename ])
% $$$ SaveForIllustrator( [FigPathCommon mfilename '/' FileOut],'jpg')
% $$$ fprintf('DONE\n')






%-------------------------------------------------------------------------------%
%----------------------- Supplementary functions -------------------------------%
%-------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------

function [varargout] = negativepeak(varargin);
%NEGATIVEPEAK is a function which finds negative peaks between adjacent 
%zero crossings in a signal. If a vector of zerocrossings is provided as 
%the third parameter, zerocrossing computing is skipped here.
% 
%INPUTS:
% SAMPLES  is a 1xN vector of signal samples
% ZEROCROSS  is a 1xM vector of zero crossing indicies in the signal 
%
%OUTPUT: 
% PEAK_ID is a vector with indicies of found peak in the signal
% PEAK_VAL is a vector with actual peak values
%
% Evgeny Resnik
% version 3.11.2006


%extract provided input parameters
if nargin==1
    samples = varargin{1};
elseif nargin==2
    samples      = varargin{1};
    zerocross_id = varargin{2};
else
    error('Wrong number of input parameters!')
end

if ~isvector(samples)
    error('Input parameter must be a vector!')
end

if size(samples,1)>size(samples,2)
    samples=samples';
end

%find zero crossing if it was not provided
if nargin<2
    zerocross_id = zerocrossing(samples);
end


%find negative peaks in a signal as a minimum between adjacent zero crossings
peak_val = [];   
peak_id  = [];  
i=1;
for lp=1:length(zerocross_id)-1  
    if (sign(samples(round((zerocross_id(lp) + zerocross_id(lp+1))/2))) == -1)
        %if interval is negative
        [peak_val(i),peak_id(i)] = min(samples(zerocross_id(lp) : zerocross_id(lp+1)));
        peak_id(i) = peak_id(i) + zerocross_id(lp)-1;
        i=i+1;
    end
end


%extract specified output parameters
if nargout==0
    figure('name','Positive peak detection');
    plot(samples); hold on; axis tight
    plot(xlim,[0 0],'k:')
    plot(peak_id,samples(peak_id),'r*')
elseif nargout==1
    varargout{1} = peak_id;
elseif nargout==2
    varargout{1} = peak_id;
    varargout{2} = samples(peak_id);
else
    error('Wrong number of output parameters!')
end


function zerocross_id = zerocrossing(samples)
%ZEROCROSS is a function which finds zero crossing of the provided signal.
%The function searches for sample elements with different sign.
%
% Evgeny Resnik
% version 1. 12.01.2005

if size(samples,1)>size(samples,2)
    samples=samples';
end

zerocross_id = find(sign(pairprod(samples))==-1);


function c=pairprod(a)
%PAIRPROD is a function which calculates pairwise product of the vector element.
%
% Evgeny Resnik
% version 1. 12.01.2005

if isvector(a)
    a2 = a(2:end);
    c = a(1:length(a2)).*a2;
else
    error('Input parameter must be a vector!')
end

%------------------------------------------------------------------------------------------






