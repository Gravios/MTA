

Trial = MTATrial();
load(['/data/homes/gravio/data/analysis/' jg05-20120310 '/' jg05-20120310 ...
      '.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )

% $$$ Trial = MTATrial('jg04-20120130');
% $$$ load(['/data/homes/gravio/data/analysis/jg04-20120130/jg04-' ...
% $$$  '20120130.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )

% $$$ Trial = MTATrial('jg04-20120213');
% $$$ load(['/data/homes/gravio/data/analysis/jg04-20120213/jg04-' ...
% $$$  '20120213.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )

% $$$ Trial = MTATrial('er01-20110721');
% $$$ load(['/data/homes/gravio/data/analysis/er01-20110721/er01-' ...
% $$$  '20110721.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )

% $$$ pfw = MTAPlaceField(Trial,[],'walk');
% $$$ pfr = MTAPlaceField(Trial,[],'rear');


tlist = {'jg05-20120309','jg05-20120310','jg05-20120315','jg05-20120317'};
adata = [];

for i=1:numel(tlist),
Trial = MTATrial(tlist{i});
load(['/data/homes/gravio/data/analysis/' tlist{i} '/' tlist{i} ...
      '.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )

adata = cat(1,adata,data);
end
data = adata;

dlen = length(data);

pes = permute(reshape([data(:).pdd_edd_snr],auxdata.nbins,auxdata.ntrans,dlen),[1,3,2]);

fpes = zeros(auxdata.nbins,dlen,auxdata.ntrans);
for s = 1:4,
    fpes(:,:,s) = -Filter0(gausswin(5)./sum(gausswin(5)),pes(:,:,s));
end

[~,mpi] = max(fpes);
mpi = sq(mpi);
[~,spi] = sort(mpi);
spi = sq(spi);


figure
for s = 1:4,
subplot(2,2,s);
imagesc(-pes(:,spi(:,4),s)');
Lines(37,[],'k');
caxis([-1,1])
end

c = 'bgrk';
for o= 1:4;
figure,
for s = 1:4,
subplot(2,2,s);
plot(mean(-pes(21:35,:,o))-mean(-pes(39:45,:,o)), ...
     mean(-pes(21:35,:,s))-mean(-pes(39:45,:,s)),['.' c(s)])
xlim([-1,1])
ylim([-1,1])
line([-1,1],[-1,1]);
end
end


