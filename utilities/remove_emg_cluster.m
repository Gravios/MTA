function [EMG_thrd, EMG_heavy, sug_period] = remove_emg_cluster(Session,data,varargin)
% [EMG_thrd, EMG_heavy] = EMG_Cluster(x,[denoiseFrequencyHighpass, dataSampleRate, lwinms, nchunks,
%                                        swin, isave, FileName])
% function to detect large EMG activity period. 
% 
% Inputs: 
%   x: signal in nt x nch
%   Optional inputs:
%       denoiseFrequencyHighpass: high pass filter to detect high museual tune, defualt: 100
%       dataSampleRate: the sampling rate, defualt: 1250 Hz
%       lwinms: window length in miliseconds to calculate the high coherence periods. 
%               default: 20 ms
%       nchunks: chunk the data into n chuncks, defualt:4. this directly
%           related to the window to compute the EMG artifact.
%           Heuristically should be less than 15 mins, aka. 1.125e6
%           samples.  
%       swin: the minimum length of the scilence window. default: 500 ms
%       isave:  if save to a file. defualt: false (when the entry is empty)
%               if the entry is not empty, the file is saved to a file
%               given by "Filename".
%       FileBase: the data is going to be aved in FileBase.EMG_Cluster.mat.
%                 defualt: current_dictionary_name
%                 
% Outputs: 
%   EMG_thrd: vector indicates the high EMG periods, used in te EMG removing.
%   EMG_heavy: vector indicates the highly coherent periods. The condition is
%           looser than the EMG_thrd. This trace is used to detect the
%           recommended working periods. 
%   sug_period: recommened working periods. Determined by the nchunks, nx2.
% 
% This is the perperocessing function of the EMG_rm_long.m to extract EMG
% activity condensed periods. 
% This function uses the StartEnd1d.m function in ./util/. 
% 
% Related functions: 
% EMG_rm_long, StartEnd1d.
% This function is a part of the EMG_removing toolbox.
% 
% errors contact: chen at biologie.uni-muenchen.de
% 
% last modified: 01.12.2019

%% ARGUEMENT COMPELETION

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('denoiseFrequencyHighpass',  100,                                               ...
                           'dataSampleRate',  1250,                                              ...
                                   'lwinms',  20,                                                ...
                                  'nChunks',  4,                                                 ...
                                     'swin',  [],                                                ...
                                    'isave',  false                                              ...
);
[denoiseFrequencyHighpass, dataSampleRate, lwinms, nChunks, swin, isave] =   ...
    DefaultArgs(varargin,defargs,'--struct');

[nt,nch] = size(data);
if nt<nch
    data = data';
    nt = size(data,1);
end

nChunks = fix(nt/min(15*60*dataSampleRate, fix(nt/nChunks)));
if isempty(swin)
    swin = fix(dataSampleRate/10);
else
    swin = fix(swin/1000*dataSampleRate);
end
%---------------------------------------------------------------------------------------------------


% PREPROCESS data 
nch = size(data,2);
tmp_data = ButFilter( data, 4, denoiseFrequencyHighpass/(dataSampleRate/2), 'high');
%clear('data');

tmp_data = zscore(tmp_data); % zscore to compute the coherence instead of covariance. 
% I guess it's not a big deal to use the covariance. 

% DETECT high emg periods
% pairwise coherence with sliding window (here I just using the mean). 
cov_data = cell(nch);
lwin = fix(lwinms*dataSampleRate/1000);
Win = ones(lwin, 1)/lwin;
for k = 1:nch
    for n = k:nch
        cov_data{k,n} = conv(tmp_data(:,k).*tmp_data(:,n), Win, 'same');
    end
end
Xtrain = cmp_Xtrain(cov_data);
% Xtrain=[ttv_v,cov_v];

% CLUSTER with ICA. 
[A, W] = fastica(Xtrain','verbose','off');
tmp_A = abs(A(1,:))./max(abs(A(2,:)),0);
if tmp_A(2)>tmp_A(1)
    A = A(:,[2 1]);
    W = W([2 1],:);
end
EMG_heavy = abs(Xtrain*W(2,:)')> abs(Xtrain*W(1,:)');
EMG_thrd =(abs(Xtrain*W(2,:)') > prctile(abs(Xtrain(~EMG_heavy,:)*W(2,:)'),99));

%% RECOMMENDED WORKING PERIOD SECTION

% suggested segment points: 
sug_period = ceil(linspace(1,nt,nChunks+1));% rough section 

% the midpoints of low EMG long periods (LEP) closest to the rough segment
% points.

se = StartEnd1d(~EMG_heavy);
se_cand = se;
se_durations = diff(se,1,2);

se(diff(se,1,2)<min(swin,max(50,prctile(se_durations,50))),:) = [];
se_cand(diff(se_cand,1,2)<min(max(swin-10,1),max(20,prctile(se_durations,10))),:) = [];% 20 ms
for k  =2:(nChunks-1)
    tmp1 = find(se(:,1)<=sug_period(k),1,'last');
    tmp2 = find(se(:,2)>=sug_period(k),1,'first');
    if tmp1==tmp2 
        % when the rough point is covered by LEP
        sug_period(k) = fix(sum(se(tmp1,:))/2);
    elseif abs(sug_period(k)- se(tmp1,2))<abs(sug_period(k)- se(tmp2,1)) 
        % when the rough point is not covered by LEP
        sug_period(k) = fix(sum(se(tmp1,:))/2);
    else
        if ~isempty(tmp2)
            sug_period(k) = fix(sum(se(tmp2,:))/2);
        else
            tmp2 = find(se_cand(:,2)>=sug_period(k),1,'first');
            sug_period(k) = fix(sum(se_cand(tmp2,:))/2);
        end
    end
end
sug_period = unique(sug_period);

% DISCARD SHORT PERIODS
short_period_threshold = 20;
while sum(diff(sug_period)<(short_period_threshold*dataSampleRate))
    sug_period(find(diff(sug_period)<(short_period_threshold*dataSampleRate))+1)=[];
end

sug_period = [[1; sug_period(2:(end-1))'+1] [sug_period(2:(end-1))';nt]];

%% SAVE RESULTS
if isave
    EMG_par.denoiseFrequencyHighpass = denoiseFrequencyHighpass;
    EMG_par.dataSampleRate           = dataSampleRate;
    EMG_par.lwinms                   = lwinms;
    EMG_par.nChunks                  = nChunks;
    EMG_par.swin                     = swin;
    EMG_par.FileBase                 = fullfile(Session.spath,Session.name);
    EMG_par.datetime                 = datetime();
    EMG_par.A                        = A;
    EMG_par.W                        = W;
    save([EMG_par.FileBase,'.EMG_Cluster.mat'], 'EMG_heavy','EMG_thrd','sug_period','EMG_par');
end

function Xtrain = cmp_Xtrain(cov_data)
c_idx = find(triu(ones(size(cov_data))));
[~,rids] = sort(rand(length(c_idx),1));
rc_idx = c_idx(rids);
cov_v = cov_data{c_idx(1)}(:);
ttv_v = zeros(length(cov_data{1}),1);
for k = 2:length(c_idx)
    cov_v = cov_v+cov_data{c_idx(k)}(:);
    % total variance to garentee line shape, against the dipole structure.
    % ttv_v = sum(abs(diff(tmp_x(:,rids),1,2)),2);
    ttv_v = ttv_v+ abs(cov_data{rc_idx(k)}(:) - cov_data{rc_idx(k-1)}(:));
end
Xtrain=[ttv_v,cov_v];