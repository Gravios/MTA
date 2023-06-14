function [rhm,varargout] = fet_rhm(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% 
% Computes a rostrocaudal motion feature of the head.
% 
% Varargin:
%           NAME : TYPE    : DEFAULT VAL          : Description
%     sampleRate : Numeric : Trial.xyz.sampleRate : computational sampling rate
%           mode : String  : {'mta','mtchglong'}  : Computational mode ( signal or spectra )
%           wsig : Logical : true                 : flag determining if signal is to be whitened
%        defspec : Struct  : def_spec_parm(xyz)   : default spetral args
%      overwrite : Logical : false                : flag for overwriting saved data
%          newSR : Numeric : []                   : Final sampling rate of output spectra
%           type : String  : mta                  : Output type
%

% DEFARGS ------------------------------------------------------------------------------------------
Trial = MTATrial.validate(Trial);
parspec = empty_spec;
%xyz = Trial.load('xyz');
xyz = preproc_xyz(Trial,'trb');
varargout = cell([1,nargout-1]);
defargs = struct('sampleRate',                            'xyz',                                 ...
                 'mode',                                  'mta',                                 ...
                 'wsig',                                  true,                                  ...
                 'defspec',                               def_spec_parm(xyz),                    ...
                 'overwrite',                             false,                                 ...
                 'newSR',                                 [],                                    ...
                 'type',                                  'mta'                                  ...
);
[sampleRate,mode,wsig,defspec,overwrite,newSR,type] = DefaultArgs(varargin,defargs,'--struct');
fs = []; ts = [];
%---------------------------------------------------------------------------------------------------


% create a ridgid body model
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[3]./(xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},hcom);

% if xyz sampling rate is greater than 120 Hz then resample it to 120 Hz
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    xyz.resample(sampleRate);    
elseif sampleRate > 120, 
    xyz.resample(sampleRate); 
    defspec.Fs = sampleRate;
end

%xyz.filter('RectFilter',3,4);
xyz.filter('ButFilter',4,30,'low');
ang = create(MTADang,Trial,xyz);
ang.data(~nniz(xyz(:,1,1)),:,:,:) = 0;

fet = xyz.copy;

rhm = MTADfet.encapsulate(Trial,...
                          [0;diff(ang(:,'head_back','fhcom',3));0],...
                          ang.sampleRate,...
                          'rhythmic head motion',...
                          'rhm',...
                          'r');

rhm.filter('ButFilter',3,[3,20],'bandpass');
rhm.data = diff(rhm.data);
rhm.update_hash();

if strcmp(mode,'mta'),  return;  end

% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
hash = DataHash(struct('sampleRate',                   sampleRate,                               ...
                       'mode',                         mode,                                     ...
                       'wsig',                         wsig,                                     ...
                       'newSR',                        newSR,                                    ...
                       'type',                         type,                                     ...
                       'xyz',                          xyz.hash)                                 ...
);
%---------------------------------------------------------------------------------------------------


switch mode
  case 'raw'
    rhm = fet.data;
  otherwise
    dsf = fieldnames(defspec);
    for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
    data = zeros(rhm.size);
    
    hash = DataHash({hash,parspec});
    
    if wsig,
        try,load(fullfile(Trial.path.arm,[mfilename,'.arm.mat']));end

        if exist('ARmodel','var')||overwrite,
            data(nniz(rhm.data),:) = WhitenSignal(rhm.data(nniz(rhm.data),:),...
                                                      [],...
                                                      true,...
                                                      ARmodel);
        else
            [data(nniz(rhm.data),:),ARmodel] = WhitenSignal(rhm.data(nniz(rhm.data),:),...
                                                                [],...
                                                                true);
            save(fullfile(Trial.path.arm,[mfilename,'.arm.mat']),'ARmodel');
        end
    else
        data(nniz(rhm.data),:) = rhm.data(nniz(rhm.data),:);
    end

     svout = cell([1,nargout-3]);
    [ys,fs,ts,svout{:}] = spec(str2func(mode),data,parspec);


    % Modify time stamps and spec; add padding (0's)
    ts = ts+(parspec.WinLength/2)/rhm.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(rhm.size(1)-round(parspec.WinLength/2),parspec.WinLength)/rhm.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

    ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

    if numel(svout)>0,
        for i = 1:numel(svout),
        svout{i} = cat(1,zeros([pad(1),size(svout{i},2),size(svout{i},3),size(svout{i},4)]),...
                         svout{i},...
                         zeros([pad(2),size(svout{i},2),size(svout{i},3),size(svout{i},4)]));
        end
    end

    if ~isempty(newSR),
        xyz.resample(newSR); % I know ...
        rhm.resample(xyz);
        temp_ts = Trial.xyz.copy;
        temp_ts.data = ts;
        temp_ts.resample(newSR);
        ts = temp_ts.data;
    end

    if strcmp(type,'raw'),
        rhm = rhm.data;
    else
        rhm.update_hash(hash);
    end
    
    
    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));

end






