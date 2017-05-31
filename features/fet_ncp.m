function [ncp,fs,ts] = fet_ncp(Trial,varargin)
%function [ncp,fs,ts] = fet_ncp(Trial,varargin)
%[sampleRate,mode,chans] = DefaultArgs(varargin,{'xyz','',2});
%SampleRate is always xyz sampleRate at the moment
%
%
%

parspec = empty_spec;
xyz = Trial.load('xyz');
varargout = cell([1,nargout-1]);

[sampleRate,mode,chans,defspec,overwrite,newSR,type] = DefaultArgs(varargin,{'xyz','mta',2,def_spec_parm(xyz),false,[],'mta'});

fs = []; ts = [];

% load NCP channel from lfp file
fet = Trial.lfp.copy;
if ~isempty(chans),
    fet.filename = [Trial.name,'.lfp'];
    fet.load(Trial,chans);
end




%SampleRate is always xyz sampleRate at the moment
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    fet.resample(sampleRate);    
elseif Trial.xyz.sampleRate<120,
    xyz = Trial.load('xyz');
    fet.resample(xyz);
else
    xyz = Trial.load('xyz').resample(120);
    fet.resample(xyz);
end


switch mode
  case 'mta'
    ncp = MTADfet.encapsulate(Trial,...
                              fet.data,...
                              fet.sampleRate,...
                              'nasal cavity pressure',...
                              'ncp',...
                              'n');
  case 'raw'
    ncp = fet.data;
  otherwise
    dsf = fieldnames(defspec);
    for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
    data = zeros(fet.size);
    
    if 0%wsig,
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
    pad = round([ts(1),size(fet,1)./fet.sampleRate-ts(end)].*ssr)-[1,0];
    %pad = round([ts(1),mod(fet.size(1)-round(parspec.WinLength/2),parspec.WinLength)/fet.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    ncp = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

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
        ncp = ncp.data;
    end

    
    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));
end




