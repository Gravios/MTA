function remove_emg_main(Session,fileExt,probeChannels,varargin)
% EMG_rm_main(Session, fileExt ,[denoise_frequency_lowerbound,
%                                savedir,denoise_shank,cleanEMGch,...
%                                silencePeriods,sp_loadingfuns,...
%                                rm_linenoise,line_thrd, ...
%                                numOfIC,hp_freq,nChunks, cmp_method,down_sample, ...
%                                nFFT, Winlength, ...
%                                save_together])
%
% The main function to perform denoising. 
% This function is working under the session directory
%
% Inputs:
%   Session: MTASession
%   fileExt: file extension
%   Optional:
%       denoise_frequency_lowerbound: keep slower signal until this frequency
%       savedir: path to save the cleaned files.
%       denoise_shank: the shanks for denoising. defualt: shank number 1,
%                      use [shk_1,...,shk_n] to denoise multiple shanks
%       cleanEMGch: the channels to detect EMG noise. defualt: 5 sampled in
%                   the first shank.
%       silencePeriods: only compute and remove the EMG from non_silence 
%                       periods, give the periods you don't want to remove 
%                       anything or their filenames. default: false 
%       sp_loadingfuns: functions to load the discarded periods. 
%       rm_linenoise: if remove line noise. default: true
%       line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.
%       Parameters used by the subfunctions: (computational details)
%           numOfIC: number of ICs.
%           hp_freq: high_pass_freq, defualt: 100 Hz
%           nChunks: number of chunks, defult: 6 or every chunk less than 15 mins,
%           cmp_method: methods to compute the EMG components. Whiten('w') or
%                       Highpass&Whiten('hw', a bit more stable). defualt: 'hw'
%           down_sample: down sample the data to compute the EMG
%                       components. default: 3.
%           nFFT: nfft to compute spectrum in EMG_rm_viewspec.m
%           WinLength: length of moving window in EMG_rm_viewspec.m
%           save_together: if save all the chunks together. defult: true
%           use_wb: whether to use the wide band signal when the algorithm
%               can't find a proper flat component at high frequency band
%               defualt: false NB: please check the report fig to see if you
%               really need this! 
%           PeriodLengthLimits:the maximum allowed length of your period.
%               default: 1e6
% 
% Related functions: 
% EMG_Cluster.m, EMG_rm_long.m, EMG_rm_pip.m, EMG_rm_viewspec.m, 
% EMG_rm_report.m, EMG_rm_viewnoise.m
%
% NB: please make sure you remove all the files generated previously when you
% try to redo the denoising.
%
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 11.12.2019.

% TEST VARS --------------------------------------------
Session = MTASession.validate('jg05-20120312.cof.all');
fileExt = 'dat';

denoiseFrequencyLowpass = 0;
denoiseFrequencyHighpass = 100;
probeId = [1,2];
probeChannels = {1:64,65:96};
channels = probeChannels{1};
savedir =  fullfile(Session.spath,'emg');
cleanEMGch = fix(linspace(3,length(channels)-3,5));

dataSampleRate = 32552;
silencePeriods =   false;
specialLoadFunc =   [];
rm_linenoise =  true;
line_thrd =  1.8;
numOfIC =  max(fix(length(channels)*.8),min(16, length(channels)));
nChunks  =  6;
cmp_method = 'hw';
down_sample =  3;
nFFT = 2^round(log2(dataSampleRate)-1);
winLength =  dataSampleRate;
save_together =  true;
use_wb =  false;
PeriodLengthLimits =  1e6;

% ------------------------------------------------------


par = LoadPar(fullfile(Session.spath,[Session.name,'.xml']));
channels = probeChannels{1};
switch fileExt
  case 'lfp'
    dataSampleRate = par.lfpSampleRate;
  case 'eeg'
    dataSampleRate = par.lfpSampleRate;
  case 'lfpinterp'
    dataSampleRate = par.lfpSampleRate;
  case 'dat'
    dataSampleRate = par.SampleRate;
  otherwise
    error('MTA:utilities:remove_emg_main:UnknownFileType');
end
    
% DEFARGS -------------------------------------------------------------------------------------------
defargs = struct(                                                                                 ...
                         'savedir',  fullfile(Session.spath,'emg'),                               ...
                      'cleanEMGch',  fix(linspace(3,length(channels)-3,5)),                       ...
                  'silencePeriods',  false,                                                       ...
                  'specialLoadFunc',  [],                                                         ...
                    'rm_linenoise',  true,                                                        ...
                       'line_thrd',  1.8,                                                         ...
                         'numOfIC',  max(fix(length(channels)*.8),min(16, length(channels))),     ...
      'denoiseFrequencyLowpass',  0,                                                              ...    
        'denoiseFrequencyHighpass',  100,                                                         ...
                         'nChunks',  6,                                                           ...
                      'cmp_method',  'hw',                                                        ...
                     'down_sample',  3,                                                           ...
                            'nFFT',  round(log2(dataSampleRate)-1),                               ...
                       'winLength',  dataSampleRate,                                              ...
                   'save_together',  true,                                                        ...
                          'use_wb',  false,                                                       ...
              'PeriodLengthLimits',  1e6                                                          ...
);
[savedir, cleanEMGch, silencePeriods, specialLoadFunc, rm_linenoise,line_thrd, numOfIC,           ...
 denoiseFrequencyLowpass, denoiseFrequencyHighpass,nChunks, cmp_method,down_sample,nFFT,       ...
 winLength, save_together, use_wb, PeriodLengthLimits] = DefaultArgs(varargin,defargs,'--struct');
%------------------------------------------------------------------------------------------

create_directory(savedir);

% LOAD auto regressive model
filePathData = fullfile(Session.spath,[Session.name,'.',fileExt]);
filePathARmodel = fullfile(Session.spath,[Session.name,'.WhitenSignal.ARmodel.',fileExt,'.mat']);

assert(logical(exist(filePathData,'file')),['MTA:utilities:remove_emg_main:DataNotFound: file: ',filePathData]);

% LOAD auto regressive model
armodel = [];
if exist(filePathARmodel,'file')
    armodel = load(filePathARmodel);
end

% COMPUTE the number of samples. 
nSamples = dir(filePathData).bytes/par.nChannels;

% CREATE denoised dat,lfp,...ect FILE as datd,lfp,...ectd 
filePathOutput = fullfile(Session.spath,[Session.name,'.',fileExt,'d']);
if ~exist(filePathOutput,'file')
    fileID = fopen(filePathOutput,'w');
    fwrite(fileID, repmat(int16(0),1,par.nChannels*nSamples(end)),'int16');
    fclose(fileID);
end

% ???? why the cell
save_range = cell(2,1);
save_range{2} = [par.nChannels, nSamples(end)];

%% DETECTING THE HIGH EMG PRIODS
filePathEMG = fullfile(Session.spath,[Session.name,'.EMG_Cluster.mat']);
if exist(filePathEMG, 'file')
    load(filePathEMG, 'EMG_thrd', 'sug_period');
    fprintf('\n\n--------------------------------------------------\n                   Please!!! \n  Remove all the previously generated files \n     when you try to redo the denoising!\n--------------------------------------------------\n\n')
    if silencePeriods
        load(filePathEMG, 'included_periods');
    else
        included_periods=[];
    end
else
    swin = 500;
    lwinms = 20;

    % LOAD the lfp file
    data = LoadBinary(filePathData, cleanEMGch, par.nChannels)';
    data = bsxfun(@minus,data,mean(data,1));
    
    fprintf('\nDetecting EMG periods...\n')
    if silencePeriods
        % SKIP quiescent periods which have little or no emg
        if ischar(silencePeriods) % load from files
            included_periods = LoadAwake(silencePeriods,nSamples(end),specialLoadFunc);
        elseif length(silencePeriods)<2 && silencePeriods % if true load from files
            included_periods = LoadAwake(sprintf('%s.sts.%s', FileBase, 'SWS'),nSamples(end),specialLoadFunc);
        elseif length(silencePeriods)<nSamples(end) % given periods
            included_periods = true(nSamples(end),1);
            for k = 1:size(silencePeriods,1)
                included_periods(silencePeriods(k,1):silencePeriods(k,2)) = false;
            end
        else % binary sample mask
            included_periods = ~silencePeriods;
        end
        silencePeriods = true;% a flag
        [EMG_thrd, ~, sug_period,included_periods] = ...
            remove_emg_cluster_s(Session, data, included_periods, denoiseFrequencyHighpass,          ...
                                 dataSampleRate, lwinms, nChunks,swin, true, PeriodLengthLimits);
    else
        % COMPUTE the whole darn thing
        [EMG_thrd, ~, sug_period] = ...
            remove_emg_cluster(Session, data, denoiseFrequencyHighpass, dataSampleRate, lwinms,      ...
                               nChunks, swin, true);
        silencePeriods = false;
        included_periods = [];
    end
    clear('data')
end

% REMOVE EMG artifacts
fprintf('\nEMG artifacts removing...\n')

%sug_period = reshape(sort(reshape([sug_period(:,1),sug_period(:,1)+fix(diff(sug_period(1,1:2)))./2,sug_period(:,1)+fix(diff(sug_period(1,1:2)))./2+1,sug_period(:,2)],[],1)),2,[])';
%armodel = [1,-0.997452338033815];
sug_period =[  1,152691143;...
             152691144,305382286;...
                 305382287
                 458073430
                 610764573
                 763455716
                 916146859
                1068838002
                1221529145
                1374220288
sug_period = fix(linspace(1,1374220288,16))';
sug_period = [sug_period(1:end-1),[sug_period(2:end-1)-1;sug_period(end)]];
               

nPeriod  = size(sug_period,1);
nProbe   = length(probeChannels);
HPs      = [];
Ws       = cell(nPeriod,nProbe);
As       = cell(nPeriod,nProbe);
EMG_au   = cell(nPeriod,nProbe);
AW       = cell(nPeriod,nProbe);
flatness = cell(nPeriod,nProbe);
for n = 1:nProbe
    HP = probeChannels{n};
    HPs = [HPs;HP(:)];
    
% WRITE speparate dat file for emg component
    myData = int16(zeros(1,nSamples(end)));
    filePathEMG = fullfile(savedir,[Session.name,'.prb-',num2str(n),'.emg']);
    if ~exist(filePathEMG,'file')
        fileID = fopen(filePathEMG,'w');
        fwrite(fileID, myData,'int16');
        fclose(fileID);
        clear myData
    end
    
    for k = 1:nPeriod
        fprintf('\rshank%d, period%d in %d...', n, k, nPeriod)
        tmp_Period = sug_period(k,:);
        save_range{1} = [tmp_Period;HP(1) HP(end)];
        if ~isempty(included_periods)
            tmp_included_periods=included_periods(tmp_Period(1):tmp_Period(2));
        else
            tmp_included_periods=[];
        end
        [~, Ws{k,n}, As{k,n}, EMG_au{k,n}, AW{k,n}, armodel,scaling_factor(k,n),flatness{k,n}]          ...
            = remove_emg_long(Session, fileExt,                                                         ...
                              LoadBinary(filePathData,HP,par.nChannels,2,[],[],tmp_Period)',            ...
                              denoiseFrequencyLowpass, silencePeriods, tmp_included_periods,            ...
                              dataSampleRate, rm_linenoise, line_thrd, denoiseFrequencyHighpass,        ...
                              EMG_thrd(tmp_Period(1):tmp_Period(2)),                                    ...
                              true, armodel, cmp_method, down_sample, numOfIC,true, save_range,         ...
                              save_together,use_wb,probeId(n),[],[], filePathEMG, filePathOutput);
    end
end
par.savedir = savedir;
par.probeChannels = probeChannels;
par.cleanEMGch = cleanEMGch;
par.silencePeriods = silencePeriods;
par.specialLoadFunc = specialLoadFunc;
par.rm_linenoise = rm_linenoise;
par.line_thrd = line_thrd;
par.numOfIC = numOfIC;
par.denoiseFrequencyLowpass = denoiseFrequencyLowpass;
par.denoiseFrequencyHighpass = denoiseFrequencyHighpass;
par.nChunks = nChunks;
par.cmp_method = cmp_method;
par.down_sample = down_sample;
par.nFFT = nFFT;
par.Winlength = winLength;
par.save_together = save_together;
par.use_wb = use_wb;
par.PeriodLengthLimits = PeriodLengthLimits;
par.PYR_Channel = 37; % ?????
if save_together
    shank_names = sprintf('%d-',probeId);
    save(sprintf('%s/%s.EMG_rm.sh%s.mat',savedir,FileBase,shank_names(1:(end-1))), ...
         'Ws','As','AW','EMG_au',                                                  ...
         'armodel','sug_period','par','probeChannels',                             ...
         'filePathData','scaling_factor','flatness');
end

% ???????????????????????????????
% COMPLETE OTHER CHANNELS
other_channels = setdiff(1:par.nChannels,HPs);
nothrch = length(other_channels);
% MEMMAP output file
m = memmapfile(filePathOutput,'Format',{'int16',save_range{2},'x'},'Writable',true);
% MEMMAP data file
mlfp = memmapfile(filePathData,'Format',{'int16',save_range{2},'x'});
for k = 1:nothrch
    m.Data.x(other_channels(k),:) = mlfp.Data.x(other_channels(k),:);
end
clear('m','mflp');
fprintf('\nDone\n')

% CHECK RESULTS:
PYR_Channel = 37;
EMG_rm_report();% ([],PYR_Channel);
EMG_rm_viewnoise();% (PYR_Channel,[])
if nSamples(end)>1.6e7
    nFFT = 2^(32-ceil(log2(nSamples(end))));
    fprintf('\nNB: we are using nFFT: %d, recompute EMG_rm_viewspec() if you want.',nFFT)
end
EMG_rm_viewspec([],[],[],[],nFFT,winLength);

% EOF