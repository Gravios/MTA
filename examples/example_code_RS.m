

% GET session metadata as a structure
sessionList = get_session_list('ZY0519');

% SET States - just for example
% states = {'locomotion','rear',...};

% LOAD Trials
Trials = af(@(s)  MTATrial.validate(s),    sessionList);


%% XYZ DATA
% LOAD xyz objects for each Trial
xyzs   = cf(@(t)  t.load('xyz'),           Trials);
% DOWNSAMPLE xyz objects to 20 Hz
cf(@(x)  x.resample(20), xyz);
% DECAPSULATE xyz data objects 
cxyz = cf(@(x)  x.data,  xyzs);
% CONCATENATION of xyz data 
cxyz = cat(1,cxyz{:});


%% STATE MATRIX
stcm = cf(@(t,x)  stc2mat(t.stc,x,states),  Trials, xyzs);
stcm = cat(1,stcm{:});

%% SPECTRAL DATA
% LOAD lfp objects for each Trial
channels = 3;
%lfps   = cf(@(t)  t.load('lfp',channels),  Trials);
% OVERWRITE AR model for 
overwriteARM = [{true}, mat2cell( false([1,numel(Trials)-1]), 1, ones([1,numel(Trials)]))];
% COMPUTE lfp spectrum for each Trial
[ys,fs,ts] = cf(@(t,o) fet_spec(t,                              ... Trial
                                t.load('lfp',channels),         ... lfp 
                                'mtchglong',                    ... mode
                                true,                           ... whiten signal
                                [],                             ... sampleRate (ignore)
                                struct('nFFT'     , 2^9, 'Fs', t.lfp.sampleRate,... spectrum parameters
                                       'WinLength', 2^7, 'nOverlap',2^7*.875,   ...
                                       'FreqRange',[1,40]),                     ...
                                'overwrite', o),                ... 
                Trials);

% RESAMPLE spectral data (probably best to down sample xyz first)
cf(@(y,x)  y.resample(x), ys,xyzs);
% DECAPSULATE xyz data objects
ys = cf(@(y) y.data, ys);
% CONCATENATE xyz data 
ys = cat(1,ys{:});

% GET total time for each spectral object
tshifts = cellfun(@(t) t(end), ts);
% CREATE cumulative time shift for concatenation
tshifts = mat2cell(cumsum([0,tshifts(1:end-1)]),1,ones(size(tshifts)));
% ADD time shifts to time vector
ts = cf(@(t,s) t+s, ts, tshifts);
% CONCATENATE time vectors
ts = cat(1,ts{:});


%% RIPPLE STUFF
%



