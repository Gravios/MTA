function bhv_spect(Trial,varargin)
% function bhv_spect(Trial,varargin)
% calculates the spectral segments around onsets/offsets of behaviors
% [behaviors,nffts,winlen,channels,frange] = DefaultArgs(varargin,{{'walk','rear'},2^11,2^10,65:96,[1,30]})
[behaviors,nffts,winlen,seglen,channels,frange] = DefaultArgs(varargin,{{'walk','rear'},2^10,2^9,8,65:96,[1,40]});

%% Compute Spectrograms around rear & walk - onset/offset
if ~isa(Trial,'MTATrial'),
    Trial = MTATrial(Trial);
end

Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;

Trial.lfp.load(Trial,channels);
wlfp = WhitenSignal(Trial.lfp.data,[],1);
nchan = numel(channels);

yr = {};
t = [];
f = [];


%if matlabpool('size')~=12,matlabpool open 12,end

for bhv = 1:length(behaviors),
    y{bhv}=[];
    for b = 1:2,
        rwlfp = GetSegs(wlfp,round(cell2mat({Trial.stc{behaviors{bhv},Trial.lfp.sampleRate}.data(:,b)})-(seglen/2).*Trial.lfp.sampleRate),seglen*Trial.lfp.sampleRate+winlen,0);
        for j = 1:nchan,
            for i = 1:size(rwlfp,2),
                [y{bhv}(:,:,i,b,j),f,t] = mtchglong(rwlfp(:,i,j),...
                                                     nffts,...
                                                     Trial.lfp.sampleRate,...
                                                     winlen,...
                                                     winlen*0.875,[],[],[],frange);
            end
        end
    end

end


t = t+diff(t(1:2))/2;

thpow_walk = sq(log10(nanmedian(y{1}(:,f<12&f>6,:,:,:),2)));
thpow_rear = sq(log10(nanmedian(y{2}(:,f<12&f>6,:,:,:),2)));
figure,imagesc(sq(nanmean(thpow_rear(:,:,1,:),2))'),axis xy

figure,imagesc(sq(nanmean(thpow_walk(:,:,1,:),2))'),axis xy


