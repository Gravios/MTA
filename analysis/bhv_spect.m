function bhv_spect(Trial,varargin)
% function bhv_spect(Trial,varargin)
% [behaviors,nffts,winlen,channels,frange] = DefaultArgs(varargin,{{'walk','rear'},2^11,2^10,65:96,[1,30]})
[behaviors,nffts,winlen,seglen,channels,frange] = DefaultArgs(varargin,{{'walk','rear'},2^11,2^10,8,65:96,[1,30]});

%% Compute Spectrograms around rear & walk - onset/offset
if ~isa(Trial,'MTATrial'),
    Trial = MTATrial(Trial,{{'lfp',channels}},'all');
end

if isempty(Trial.lfp),
    Trial = Trial.load_lfp(channels);
end

wlfp = WhitenSignal(Trial.lfp,[],1);
nchan = length(channels);

yr = {};
t = [];
f = [];


for bhv = 1:length(behaviors),
    y{bhv}=[];
    for b = 1:2,
        rwlfp = GetSegs(wlfp,...
                        round(((Trial.Bhv.getState(behaviors{bhv}).state(:,b)./Trial.xyzSampleRate)-round(seglen/2)).*Trial.lfpSampleRate),...
                        seglen*Trial.lfpSampleRate,0);
        for j = 1:nchan,
            for i = 1:size(rwlfp,2),
                [y{bhv}(:,:,i,b,j),f,t] = mtchglong(rwlfp(:,i,j),...
                                                     nffts,...
                                                     Trial.lfpSampleRate,...
                                                     winlen,...
                                                     winlen*0.875,[],[],[],frange);
            end
        end
    end

end

t = t+diff(t(1:2))/2;

thpow_rear = sq(median(yr(:,f<12&f>6,:,:,:),2));
thpow_walk = sq(median(yw(:,f<12&f>6,:,:,:),2));
