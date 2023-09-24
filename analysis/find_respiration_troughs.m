function  [amins,ncpSampleRate] = find_respiration_troughs(Trial,varargin)

overwrite = false;
if numel(varargin)==1
    overwrite = varargin{1};
end

filepath = fullfile(Trial.spath, [Trial.filebase,'.respiration_troughs.mat']);

if ~exist(filepath) | overwrite
ncp = Trial.load('lfp',Trial.meta.channelGroup.respiration);
spow = conv(sqrt(get(filter(copy(ncp),'ButFilter',4,[6,12],'bandpass'),'data').^2),...
            gausswin(128)./sum(gausswin(128)),'same');
rpow = conv(sqrt(get(filter(copy(ncp),'ButFilter',4,[2],'low'),'data').^2),...
            gausswin(128)./sum(gausswin(128)),'same');

ncpFiltHigh = filter(copy(ncp),'ButFilter',4,[16],'low');
ncpFiltLow = filter(copy(ncp),'ButFilter',4,[2],'low');

% COMPUTE ncp derivative
ncpFiltHighUnity = copy(ncpFiltHigh);
ncpFiltHighUnity.data = nunity(ncpFiltHighUnity);
ncpFiltDerUnity = copy(ncpFiltHigh);
ncpFiltDerUnity.data = nunity(ncpFiltHigh.data);
ncpFiltDerUnity.data = circshift(ncpFiltDerUnity.data,-1)-circshift(ncpFiltDerUnity.data,1);
[dmins,dmvals] = LocalMinima(abs(ncpFiltDerUnity.data),15,0.5);
dmvals = dmvals(spow(dmins)>rpow(dmins));
dmins = dmins(spow(dmins)>rpow(dmins));

% SELECT trough minimas
hfig = figure(50001),
plot(ncpFiltHighUnity(dmins),dmvals,'.')
waitforbuttonpress();
waitforbuttonpress();
waitforbuttonpress();    
[cpnts] = ClusterPP(hfig);

% SELECT trough points from the manual clustering above
cmins = dmins(cpnts==1);
cvals = dmvals(cpnts==1);
close(figure(50001));

% COMPUTE ncp derivative
ncpFiltLowUnity = copy(ncpFiltLow);
ncpFiltLowUnity.data = nunity(ncpFiltLowUnity);
ncpFiltDerUnity = copy(ncpFiltLow);
ncpFiltDerUnity.data = nunity(ncpFiltLow.data);
ncpFiltDerUnity.data = circshift(ncpFiltDerUnity.data,-1)-circshift(ncpFiltDerUnity.data,1);
[dmins,dmvals] = LocalMinima(abs(ncpFiltDerUnity.data),150,0.5);
dmvals = dmvals(spow(dmins)<rpow(dmins));
dmins = dmins(spow(dmins)<rpow(dmins));

% SELECT trough minimas
hfig = figure(50001),
plot(ncpFiltLowUnity(dmins),dmvals,'.')
waitforbuttonpress();
waitforbuttonpress();
waitforbuttonpress();    
[lpnts] = ClusterPP(hfig);

lmins = dmins(lpnts==1);
lvals = dmvals(lpnts==1);


amins = sort([cmins;lmins]);

binds = circshift(1./(diff(amins)./ncp.sampleRate)>16,1);
ncpFiltDerUnity = copy(ncpFiltHigh);
ncpFiltDerUnity.data = nunity(ncpFiltHigh.data);
ncpFiltDerUnity.data = circshift(ncpFiltDerUnity.data,-1)-circshift(ncpFiltDerUnity.data,1);


amins = amins(~binds);
amins = amins(ncpFiltHighUnity(amins)<-0.1);
ncpSampleRate = ncp.sampleRate;

save(filepath,'amins','ncpSampleRate');
else
    load(filepath);
end
