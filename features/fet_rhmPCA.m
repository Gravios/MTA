function [rhm,varargout] = fet_rhmPCA(Trial,varargin)
% [rhm,fs,ts] = fet_rhmPCA(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate
%
% 
% fet_rhmPCA is for capturing the movements of the head using
% moving time embedded pca 
%

% DEFARGS ------------------------------------------------------------------------------------------
xyz = Trial.load('xyz');
%xyz = Trial.load('xyz','seh');
%xyz.data(isnan(xyz.data(:)))=0;
varargout = cell([1,nargout-1]);

parspec = struct('nFFT',2^9,...
                 'Fs',  xyz.sampleRate,...
                 'WinLength',2^8,...
                 'nOverlap',2^8*.875,...
                 'NW',3,...
                 'Detrend',[],...
                 'nTapers',[],...
                 'FreqRange',[.1,20]);



[sampleRate,mode,wsig,defspec,overwrite,type,markers] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'mta',1,parspec,false,'mta',{'head_back','head_left','head_front','head_right'}});
fs = []; ts = [];


flxyz = xyz.copy;
flxyz.filter('ButFilter',3,2,'low');

dflxyz = xyz.copy;
dflxyz.data = xyz(:,markers,:)-flxyz(:,markers,:);
%dflxyz.filter('ButFilter',3,20,'low');



dflxyz.resample(60);

nfet = size(dflxyz,2).*size(dflxyz,3);

embeddingWindow = 80;
sdx = dflxyz.copy;
sdx.data = dflxyz.segs([],embeddingWindow);
sdx.data = reshape(permute(sdx.data,[2,1,3,4]),size(sdx,2),embeddingWindow,nfet);

U = zeros([size(sdx,1),80,nfet]);
V = zeros([size(sdx,1),nfet,nfet]);
try
    for i = 1:size(sdx,1),
        [U(i,:,:),~,V(i,:,:)] = svd(sq(sdx(i,:,:)),0);
    end
end

Vmask = ones(size(V,1),nfet);
Vmask(1:2:end,:) = 0;

pv = V(:,:,1);
for i = 2:size(V,1)
    dpv = sqrt(sum([pv(i,:)-pv(i-1,:);pv(i,:)+pv(i-1,:)].^2,2));
    [~,dind] = min(dpv);
    if dind == 2,
        pv(i,:) = -pv(i,:);
    end
end

fet = MTADfet.encapsulate(Trial,[],dflxyz.sampleRate,'Slow Rhythmic Respiration','rhmPCA','s');
fet.data = zeros([size(dflxyz,1),1]);
for i = 1:size(V,1)
    dfx = dflxyz(i,:,:);
    fet.data(i) = dfx(:)'*sq(pv(i,:))';
end    

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


    rhm.resample(sampleRate);
    temp_ts = Trial.xyz.copy;
    temp_ts.data = ts;
    temp_ts.resample(sampleRate);
    ts = temp_ts.data;

    if strcmp(type,'raw'),
        rhm = rhm.data;
    end
    
    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));

end






