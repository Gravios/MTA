function [sfet] = bhv_shakeH(Trial,varargin)
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

[mode,thresh,overwrite] = DefaultArgs(varargin,{'per',5,false},false);


xyz = Trial.load('xyz');
ix = nniz(xyz(:,1,1));
xyz.data = ButFilter(xyz.data,3,[55]/(xyz.sampleRate/2),'low');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
hcom = xyz.com(rb);
bcom = xyz.com(rbb);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('fbcom',[.7,1,.7],{{'spine_lower','spine_middle',[0,0,1]}},ButFilter(bcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
ang = create(MTADang,Trial,xyz);



fet = xyz.copy;
fet.ext = 'fet';
fet.label = 'shake';
fet.key = 'k';
fet.data = [diff(circ_dist(circ_dist(ang(:,3,11,1),ang(:,3,4,1)),pi).*double(ix)),...
            diff(circ_dist(circ_dist(ang(:,3,5,1),ang(:,3,10,1)),0).*double(ix)),...
            diff(circ_dist(circ_dist(ang(:,1,11,1),ang(:,2,4,1)),0).*double(ix))];;

sparm = struct('nFFT',2^7,...
               'Fs',ang.sampleRate,...
               'WinLength',2^5,...
               'nOverlap',2^5*.875,...
               'FreqRange',[10,25]);
[ys,fs,ts] = fet_spec(Trial,fet,...
                      'mode','mtchglong',...
                      'wsig',true,...
                      'defspec',sparm,...
                      'overwrite',false);

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

