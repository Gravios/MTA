function [rhm,varargout] = fet_spec(Trial,fet,varargin)
% [rhm,varargout] = fet_spec(Trial,fet,varargin)
%
%[mode,wsig,sampleRate,defspec,overwrite] = ...
%    DefaultArgs(varargin,{'raw',true,120,...
%                    struct('nFFT',2^9,'Fs',fet.sampleRate,...
%                           'WinLength',2^7,'nOverlap',2^7*.875,...
%                           'FreqRange',[1,40]),...
%                    true});
%
%
parspec = empty_spec;

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('mode',                    'mtchglong',                                         ...
                 'wsig',                    true,                                                ...
                 'sampleRate',              120,                                                 ...
                 'defspec',                 def_spec_parm(fet),                                  ...
                 'overwrite',               true                                                 ...
);
[mode,wsig,sampleRate,defspec,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

varargout = cell([1,nargout-1]);

if sampleRate<fet.sampleRate&&~strcmp(fet.ext,'lfp'),
    fet.resample(sampleRate); 
else
    sampleRate = fet.sampleRate;
end



switch mode

  case 'mta'
    rhm = MTADlfp('data',fet.data,...
                  'sampleRate',fet.sampleRate,...
                  'syncPeriods',fet.sync.copy,...
                  'syncOrigin',fet.origin);
  case 'raw'
    rhm = fet.data;
  otherwise % mode {mtcsdglong,mtchglong,...}
    dsf = fieldnames(defspec);
    for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end

    data = zeros(fet.size);
    
    if wsig,
        if ischar(wsig),
            try,load(fullfile(Trial.path.arm,wsig));end
        else
            try,load(fullfile(Trial.path.arm,[mfilename,fet.label,'_' fet.key '.arm.mat']));end
        end
        
        if ~exist('ARmodel','var'), overwrite = true; end
        
        if overwrite,
            [data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                                [],...
                                                                true);
            save(fullfile(Trial.path.arm,[mfilename,fet.label,'_' fet.key '.arm.mat']),'ARmodel');
        else
            data(nniz(fet.data),:) = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                      [],...
                                                      true,...
                                                      ARmodel);
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
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

    %ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));        
    ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

    if numel(svout)>0,
        for i = 1:numel(svout),
        svout{i} = cat(1,zeros([pad(1),size(svout{i},2),size(svout{i},3),size(svout{i},4)]),...
                         svout{i},...
                         zeros([pad(2),size(svout{i},2),size(svout{i},3),size(svout{i},4)]));
        end
    end

    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));
    
end


end

% END MAIN -----------------------------------------------------------------------------------------