Trial = MTATrial('jg05-20120317');

pfw = MTAPlaceField(Trial,[],'walk');
pfr = MTAPlaceField(Trial,[],'rear');


load(['/data/homes/gravio/data/analysis/jg05-20120317/jg05-' ...
 '20120317.cof.all.cufrccg.5000_64.pfc_walk.sur_walk.mat'] )


dlen = length(data);

prs = zeros(auxdata.nbins,dlen,auxdata.ntrans);
hrs = zeros(auxdata.nbins,dlen,auxdata.ntrans);
pks = zeros(auxdata.nbins,dlen,auxdata.ntrans);
hks = zeros(auxdata.nbins,dlen,auxdata.ntrans);
for u = 1:dlen
for p1 = 1:auxdata.ntrans,
for t = 1:auxdata.nbins,
[prs(t,u,p1),hrs(t,u,p1)] = ranksum(sq(data(u).pdd(t,p1,:)), ...
                                    sq(data(u).edd(t,p1,:)));
[hks(t,u,p1),pks(t,u,p1)] = kstest2(sq(data(u).pdd(t,p1,:)), ...
                                    sq(data(u).idd(t,p1,:)));

end
end
end

fpks = reshape(Filter0(gausswin(5)./sum(gausswin(5)),-log10(eps+pks)),[],dlen,auxdata.ntrans);
[~,mpi] = max(fpks);
mpi = sq(mpi);
[~,spi] = sort(mpi);
spi = sq(spi);

figure
for s = 1:4;
subplot(2,2,s)
imagesc(-log10(eps+pks(:,spi(:,s),s))');axis xy
end

s = 1;

[ni,xi] = hist(sq(data([data(:).clu]==u).idd(t,s,:)),100);
[ns,xs] = hist(sq(data([data(:).clu]==u).idd(t,5,:)),100);
[np,xp] = hist(sq(data([data(:).clu]==u).pdd(t,s,:)),100);
[ne,xe] = hist(sq(data([data(:).clu]==u).edd(t,s,:)),100);

figure,plot(xi,ni,'-b',xp,np,'-r',xs,ns,'-g',xe,ne,'-k')
Lines(data([data(:).clu]==u).idd_p95(t,s),[],'k'); 
Lines(data([data(:).clu]==u).idd_p05(t,s),[],'k'); 
Lines(data([data(:).clu]==u).pdd_p05(t,s),[],'g'); 
Lines(data([data(:).clu]==u).pdd_p95(t,s),[],'g'); 



figure,plot(-data([data(:).clu]==u).pdd_edd_snr(:,s)); 


s =2

figure,plot(auxdata.tbins,-log10(eps+prs(:,s)),'-b',auxdata.tbins,-log10(eps+pks(:,s)),'-r')
axis tight,ylim([0,0.000001])