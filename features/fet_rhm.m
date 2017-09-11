function [rhm,varargout] = fet_rhm(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate
%
% 


% DEFARGS ------------------------------------------------------------------------------------------
parspec = empty_spec;
xyz = Trial.load('xyz');
varargout = cell([1,nargout-1]);
defargs = struct('sampleRate',                            'xyz',                                 ...
                 'mode',                                  'mta',                                 ...
                 'wsig',                                  true,                                  ...
                 'defspec',                               def_spec_parm(xyz),                    ...
                 'overwrite',                             false,                                 ...
                 'newSR',                                 [],                                    ...
                 'type',                                  'mta'                                  ...
);
{'xyz','mta',1,def_spec_parm(xyz),false,[],'mta'}
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


% if xyz sampling rat e is greater than 120 Hz then resample it to 120 Hz
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    xyz.resample(sampleRate);    
elseif sampleRate > 120, 
    xyz.resample(120); 
    defspec.Fs = 120;
end
xyz.filter('ButFilter',3,20,'low');


ang = create(MTADang,Trial,xyz);
ang.data(~nniz(xyz(:,1,1)),:,:,:) = 0;


fet = xyz.copy;

bang = ButFilter(ang(:,'head_back','fhcom',3),3,[.1,20]./(ang.sampleRate/2),'bandpass');

fet.data = [0;ButFilter(diff(bang),3,[.5,30]/(ang.sampleRate/2),'bandpass')];

switch mode
  case 'mta'
    rhm = MTADfet.encapsulate(Trial,...
                              fet.data,...
                              fet.sampleRate,...
                              'rhythmic head motion feature',...
                              'rhm',...
                              'r');
  case 'raw'
    rhm = fet.data;
  otherwise
    dsf = fieldnames(defspec);
    for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
    data = zeros(fet.size);
    
    if wsig,
        try,load(fullfile(Trial.path.arm,[mfilename,'.arm.mat']));end

        if exist('ARmodel','var')||overwrite,
            data(nniz(fet.data),:) = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                      [],...
                                                      true,...
                                                      ARmodel);
        else
            [data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                                [],...
                                                                true);
            save(fullfile(Trial.path.arm,[mfilename,'.arm.mat']),'ARmodel');
        end
    else
        data(nniz(fet.data),:) = fet.data(nniz(fet.data),:);
    end

     svout = cell([1,nargout-3]);
    [ys,fs,ts,svout{:}] = spec(str2func(mode),data,parspec);


    % Modify time stamps and spec; add padding (0's)
    ts = ts+(parspec.WinLength/2)/fet.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(fet.size(1)-round(parspec.WinLength/2),parspec.WinLength)/fet.sampleRate].*ssr)-[1,0];
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
        rhm.resample(newSR);
        temp_ts = Trial.xyz.copy;
        temp_ts.data = ts;
        temp_ts.resample(newSR);
        ts = temp_ts.data;
    end

    if strcmp(type,'raw'),
        rhm = rhm.data;
    end
    
    
    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));

end






