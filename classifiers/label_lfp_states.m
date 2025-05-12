function label_lfp_states(FileBase,varargin)
% function CheckEegStates(FileBase,State, AuxData,FreqRange,Channels,Window, Action, Overwrite)
% GUI for browsing the LFP spectrum
%   - segmentation of states
%   - precompute the spectrograms using EegSegmentation function. 
%
% Window - in sec size of the spectral window
% FreqRange - [0, 25], frequency range
% AuxData - additional timeseries to plot
%     cellArray where each row is :
%         - xaxis(Time Seconds),
%         - yaxis,
%         - data,
%         - display_func[ plot, imagesc ]
%               'plot': yaxis may be empty.
%
% Action - 'display' or 'compute'
%
% Keymap :
%     t (default) - select position in time
%     n - ADD new State periods (two left clicks: start and stop) 
%     m - SHIFT border (click and drag)
%         toggles t,n,m are not really neded if you have 3-button mouse, much
%         faster is to use the mouse buttons as described below.
%
%     d - delete border. (left click)
%     z - toggle zoom.
%             Mouse : left (zoom in) /right (zoom out)
%             Keybd : f - resets x axis to max. arrows up/down- zoom in/out the current position',...
%     c - and then  up/down - rescale the color axis in spectrograms',...
%     f - and then  up/down - zoom in/out the y axis in spectrograms',...
%     w - and then up/down - decrease/increas the size of the window for the traces display ',...
%     h - shows help box
%     u - update periods. Removes deleted and reorders all lines
%     s - save the results (currently ASCII two columnt, 1250 s.r.)
%     l - load periods from a file
%     L/R Arrows - move in time
%     space - move to the next start of a period', ....
%     a - automatic theta segementation tool.
%         GUI - A menu to select passband/stopband and channel to use 
%               up. The stopband can be more than one set of bands. The algorythm
%               will compute ratio of power in pass/stop bands and perform HMM two
%               stage fit. New segments that correspond to the state with high values
%               of the power ratio will appear. You may have to wait a bit.
% Mouse: click in trace (default) mode: 
%     left - view trace at the point of click
%     right - add new border, middle - move the border
%
%Toolboxes necessary for CheckEegStates:
% addpath(genpath('/storage/share/matlab/Third-Party_Toolboxes/HMM/hmmbox'));
% addpath(genpath('/storage/share/matlab/Third-Party_Toolboxes/netlab'));
% NOTE: A toolbox netlab has a function 'kmeans', which can shadow the built-in matlab function with the same name.
% NOTE: A toolbox chronux_2_12 shadows some unknown yet function used in CheckEegtates.
%


[State,AuxData,FreqRange,Channels,Window,Action, Overwrite] = ...
    DefaultArgs(varargin,{[],[],[1 100],[],1,'display',0});

% auxil. struct for gui
global gCheckEegStates
gCheckEegStates = struct;


UseRulesBefore = 0; % flag to switch heuristics for periods cleaning after automatic segmentation
MinLen=5; %seconds
Compute=0; 

Par = LoadPar([FileBase '.xml']);

if isfield(Par,'lfpSampleRate')
    eSampleRate = Par.lfpSampleRate;
else
    eSampleRate = 1250;
end

if ~isempty(State)
    if FileExists([FileBase '.sts.' State])
        Per = load([FileBase '.sts.' State]);
    else
        Per = [];
    end       
    if isempty(Per)
        fprintf('No segments detected in %s state\n',State)
        Per = [];
    end
else
    Per = [];
    State = '';
end

if isempty(Channels) & FileExists([FileBase '.eegseg.par']);
    Channels = load([FileBase '.eegseg.par'])+1;
end
if ~isempty(Channels)  & ~FileExists([FileBase '.eegseg.par']);
    fprintf('Your selection of Channels is stored in the file %s. Next time you do not have to pass Channels.\n',...
        [FileBase '.eegseg.par']);
    msave([FileBase '.eegseg.par'],Channels-1);
end

if FileExists([FileBase '.eegseg.mat']) 
    switch Action
        case 'display'
            if Overwrite ==0
                over='No';
            else
                over = questdlg('Do you want to overwrite existing spectrogram?','Overwrite');
            end
            switch over
                case 'No'
                    load([FileBase '.eegseg.mat']); % load whitened spectrograms from EegSegmentation results
                case 'Yes'
                    Compute=1;
            end
        case 'compute'
            if Overwrite
                Compute=1;
            else
                return;
            end
    end
end

if ~FileExists([FileBase '.eegseg.mat']) | Compute==1

    if isempty(Channels)
        ch = inputdlg({'Enter channels to use'},'Channels',1,{'1'});
        Channels = str2num(ch{1});
    end
    % now compute the spectrogram
    if FileExists([FileBase '.lfp'])
        try 
            Eeg = LoadBinary([FileBase '.lfp'], Channels,Par.nChannels,4)';
        catch
            Eeg = LoadBinary([FileBase '.lfp'], Channels,Par.nChannels,2)';
        end
    elseif FileExists([FileBase '.eeg'])
        try 
            Eeg = LoadBinary([FileBase '.eeg'], Channels,Par.nChannels,4)';
        catch
            Eeg = LoadBinary([FileBase '.eeg'], Channels,Par.nChannels,2)';
        end
    elseif FileExists([FileBase '.eeg.0'])
        Eeg = bload([FileBase '.eeg.0'],[1 inf]);
        
    else
        error('no eeg file or eeg.0 file! \n');
    end
    fprintf('computing spectrograms, may take time ... \n');
    SpecWindow = 2^round(log2(Window*eSampleRate));% choose window length as power of two
    nFFT = SpecWindow*4;
    weeg = WhitenSignal(Eeg,eSampleRate*2000,1);
    [y,f,t]=mtcsglong(weeg,nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
    save([FileBase '.eegseg.mat'],'y','f','t','Channels','-v7.3');
end

if strcmp(Action,'compute')
    return;
end
t = (t(2)-t(1))/2 +t;

if UseRulesBefore
    switch State
        case 'REM'

    end
end

% fill the global structure for future use
if ~FileExists([FileBase '.lfp']) & FileExists([FileBase '.eeg.0'])
    gCheckEegStates.eegFile  ='eeg.0';
    gCheckEegStates.Channels =1;
    gCheckEegStates.nChannels = 1;
    gCheckEegStates.nSamples = FileLength([FileBase '.eeg.0'])/2;
else
    gCheckEegStates.eegFile  ='eeg';
    gCheckEegStates.Channels = Channels;
    gCheckEegStates.nChannels = length(Channels);
    gCheckEegStates.nSamples = FileLength([FileBase '.lfp'])/Par.nChannels/2;
end

nAuxData = max(size(AuxData,1));

gCheckEegStates.FileBase = FileBase;
gCheckEegStates.Par = Par;
gCheckEegStates.State = State;
gCheckEegStates.t = 10; %is seconds
gCheckEegStates.eFs = eSampleRate;
gCheckEegStates.trange = [t(1) t(end)];
gCheckEegStates.Periods = Per/eSampleRate; % in seconds
gCheckEegStates.Mode = 't';
gCheckEegStates.nPlots=gCheckEegStates.nChannels+1+nAuxData;
gCheckEegStates.lh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.Window = Window*eSampleRate*2;
gCheckEegStates.SelLine=[];
gCheckEegStates.cposh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.FreqRange = [min(f) max(f)];
gCheckEegStates.newl=[];
gCheckEegStates.tstep = t(2)-t(1);
gCheckEegStates.coolln = [];
gCheckEegStates.LastBut = 'normal';
gCheckEegStates.nAuxData = nAuxData;
gCheckEegStates.RefPoints = [];
if nAuxData>0
    gCheckEegStates.AuxDataType = AuxData(:,4);
end

% create and configure the figure
gCheckEegStates.figh = figure('ToolBar','none');
set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
set(gCheckEegStates.figh, 'NumberTitle', 'off');


% put the uitoolbar and uimenu definitions here .. may require rewriting
% some callbacks as actions rather then cases of actions (e.g. key
% pressing)


% PLOT 
% for color scaling use median+/- std across all freq. bins
ry = reshape(log(y),[],1); ry = ry(~isinf(ry) & ~isnan(ry));

med = median(ry);
dev = std(ry); 
if isnan(dev) dev=med/3; end

for ii=1:gCheckEegStates.nChannels
    subplot(gCheckEegStates.nPlots,1,ii);
    imagesc(t,f,log(sq(y(:,:,ii)))');
    axis xy;
    ylim([max(0,f(1)) min(f(end),20)]);    
    caxis(med+[-3 3]*dev); %ORIGINAL VERSION
    hold on
    if ii==1
        title(sprintf('Spectrogram on ch-%d',gCheckEegStates.Channels)); %Evgeny Resnik: added channel index, 03.2016
        ylabel('Frequency (Hz)');
    end
end

if nAuxData>0
    for ii=[1:nAuxData]
        subplot(gCheckEegStates.nPlots,1,ii+gCheckEegStates.nChannels);
        plot_aux_data(AuxData(ii,:));
        xlim(gCheckEegStates.trange);
        hold on
        
    end
end

subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots)
label_lfp_states_keymap('traces'); % plot the current eeg traces
hold('on');

%plots lines
if ~isempty(Per)
    label_lfp_states_keymap('lines');
end

% assign functions for mouse and keyboard click
set(gCheckEegStates.figh,'WindowButtonDownFcn','label_lfp_states_keymap(''mouseclick'')');
set(gCheckEegStates.figh,'KeyPressFcn', 'label_lfp_states_keymap(''keyboard'')');

return




function [thratio, ThePeriods] = auto_theta_delta_ratio(y,f,t,MinLen)

%automatic theta periods detection just using the thetaratio
thfin = find(f>6 & f<9);
thfout = find(f<5 | (f> 12& f<25));
thratio = log(mean(sq(y(:,thfin,1)),2))-log(mean(sq(y(:,thfout,1)),2));

if nargout>1
    nStates =2;
    % fit gaussian mixture and to HMM - experimental version .. uses only thetaratio
    [TheState thhmm thdec] = gausshmm(thratio,nStates,1,0);

    for i=1:nStates
        thratio_st(i) = mean(thratio(TheState==i));
    end

    [dummy TheInd] = max(thratio_st);
    InTh = (TheState==TheInd);
    DeltaT = t(2)-t(1);
    MinSeg = round(MinLen/DeltaT);

    TransTime = ThreshCross(InTh,0.5,MinSeg);
    ThePeriods = t(TransTime);
end
return

function plot_aux_data(Data)
        nEl =  size(Data,2);
        if nEl<4 
            err=1;
        elseif ~isstr(Data{4})
            err=1;
        else
            err=0;
        end
        if err 
            warning('AuxData has to be cell array where each row is : xaxis, yaxis, data, display_func');
            close 
            return;
        end
        
        switch Data{4} %switch by functions
            case 'plot'
                plot(Data{1}, Data{3});
            case 'imagesc'
                if length(Data{1})~=size(Data{3},1) & length(Data{1})~=size(Data{3},2)
                    Data{3}=Data{3}';
                end
                imagesc(Data{1},Data{2}, Data{3}');
                axis xy
            otherwise 
                error('wrong data display function');
        end

return