function [sfet] = bhv_sniff(Trial,varargin)
% function [sper] = bhv_shake(Trial,varargin)
% [mode,thresh,overwrite] = DefaultArgs(varargin,{'per',5,false},false);
% 
%  varargin:
%
%    mode:    string  - def('per'):{mta,raw,per}, output type
%
%    thresh:  numeric - def(5), Threshold for classification
%
%    overwrite: logical - def(false), overwrite the normalization
%                         coeficients
%

[mode,thresh,overwrite] = DefaultArgs(varargin,{'per',-3,false},false);





fet = Trial.xyz.copy;
fet.ext = 'fet';
fet.label = 'sniff';
fet.key = 'f';
fet.data = fet_rhm(Trial,[],'raw');
% $$$ fet.data = [diff(circ_dist(circ_dist(ang(:,3,11,1),ang(:,3,4,1)),pi).*double(ix)),...
% $$$             diff(circ_dist(circ_dist(ang(:,3,5,1),ang(:,3,10,1)),0).*double(ix)),...
% $$$             diff(circ_dist(circ_dist(ang(:,1,11,1),ang(:,2,4,1)),0).*double(ix))];;

sparm = struct('nFFT',2^9,...
               'Fs',fet.sampleRate,...
               'WinLength',2^7,...
               'nOverlap',2^7*.875,...
               'FreqRange',[1,35]);
[ys,fs,ts] = fet_spec(Trial,fet,...
                      'mode','mtchglong',...
                      'wsig',true,...
                      'defspec',sparm,...
                      'overwrite',true);


ang = create(MTADang,Trial,Trial.load('xyz'));
ys.resample(ang);
ys.data(ys.data<0) = 1e-10;


figure(2),
ind = Trial.stc{'w'};
hist2([log10(ys(ind,40)),ang(ind,5,7,2)],linspace(-6,-2,100),100)
title({Trial.filebase,'RHM 10Hz power Vs Head pitch during walking'})
ylabel('Head pitch (radians)')
xlabel('RHM 10Hz power')



% Load Normalization parameters
try,load([mfilename('fullpath') '-MTAC_bhv_model_' fet.label '_' fet.key '.mat']),end
if ~exist('ysMean','var')||overwrite;
    ysMean = [];
    ysStd = [];
end

%FIX THIS you bald monkey!!
[ys.data,ysMean,ysStd] = nunity(cat(2,median(log10(ys(:,:,1,1)),2),...
                                       median(log10(ys(:,:,2,2)),2),...
                                       median(log10(ys(:,:,3,3)),2)),...
                                 [],ysMean,ysStd);

% Set normalization parameters
if ~exist('sMean','var')||overwrite;
    sMean = [];
    sStd = [];
else
    [~,sMean,sStd] = nunity(prod(ys.data,2),[],sMean,sStd);
end

% upsample to xyz sample rate
ys.resample(xyz);
data = nunity(prod(ys.data,2),[],sMean,sStd);

if overwrite,
    save([mfilename('fullpath'),...
        '-MTAC_bhv_model_' fet.label '_' fet.key '.mat'],...
        'ysMean','ysStd','sMean','sStd'),
end

switch mode
    case 'mta'
        sfet = MTADfet('path',       Trial.spath,...
                       'filename',   Trial.filebase,...
                       'data',       data,...
                       'sampleRate', xyz.sampleRate,...
                       'syncPeriods',Trial.sync.copy,...
                       'syncOrigin', Trial.sync.data(1),...
                       'label',      'shake',...
                       'key',        'k');
                       
    case 'raw'
        sfet = data;
    case 'per'
        sfet = MTADepoch('path',        Trial.spath,...
                         'filename',    Trial.filebase,...
                         'data',        ThreshCross(data,thresh,round(.1*xyz.sampleRate)),...
                         'sampleRate',  xyz.sampleRate,...
                         'syncPeriods', Trial.sync.copy,...
                         'syncOrigin',  Trial.sync.data(1),...
                         'label',       'shake',...
                         'key',         'k');
        assert(~isempty(Trial.stc.gsi('rear')),'MTA:classifiers:bhv_shake:REAR_PERIODS_ABSENT');
        rper = Trial.stc{'rear'};
        rper.resample(sfet);
        sfet = sfet-rper.data;        
end

