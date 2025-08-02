function [rhm,varargout] = fet_rbm(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate
%
% 
%

% DEFARGS ------------------------------------------------------------------------------------------
varargout = cell([1,nargout-1]);

Trial = MTATrial.validate(Trial);

parspec = empty_spec;
xyz = Trial.load('xyz');

defargs = struct('sampleRate',                            'xyz',                                 ...
                 'mode',                                  'mta',                                 ...
                 'wsig',                                  true,                                  ...
                 'defspec',                               def_spec_parm(xyz),                    ...
                 'overwrite',                             false,                                 ...
                 'newSR',                                 [],                                    ...
                 'type',                                  'mta'                                  ...
);

[sampleRate,mode,wsig,defspec,overwrite,newSR,type] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

fs = []; ts = [];

% if xyz sampling rat e is greater than 120 Hz then resample it to 120 Hz
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    xyz.resample(sampleRate);    
elseif sampleRate > 120, 
    xyz.resample(120); 
end

ang = create(MTADang,Trial,xyz);
ang.data(~nniz(xyz(:,1,1)),:,:,:) = 0;

fet = xyz.copy;

fet.data = [ang(:,'spine_lower','spine_middle',3)-ang(:,'spine_lower','spine_upper',3)];
fet.filter('RectFilter',3,5);
fet.data = diff(fet.data);
fet.filter('RectFilter',3,5);
fet.data = [0;diff(fet.data);0];
fet.data(~nniz(fet.data(:))) = 1;

switch mode
  case 'mta'
    rhm = MTADfet.encapsulate(Trial,...
                              fet.data,...
                              fet.sampleRate,...
                              'rhythmic body motion',...
                              'rbm',...
                              'b');
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






