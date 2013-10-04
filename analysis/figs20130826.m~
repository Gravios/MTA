Trial = MTATrial('jg05-20120310');

pfw = MTAPlaceField(Trial,[],'walk');
pfr = MTAPlaceField(Trial,[],'rear');


subplot(4,


t = 40;
u = 24;


prs = zeros(auxdata.nbins,auxdata.ntrans);
hrs = zeros(auxdata.nbins,auxdata.ntrans);
pks = zeros(auxdata.nbins,auxdata.ntrans);
hks = zeros(auxdata.nbins,auxdata.ntrans);
for p1 = 1:auxdata.ntrans,
for t = 1:auxdata.nbins,
[prs(t,p1),hrs(t,p1)] = ranksum(sq(data([data(:).clu]==u).pdd(t,p1,:)), ...
                    sq(data([data(:).clu]==u).edd(t,p1,:)));
[hks(t,p1),pks(t,p1)] = kstest2(sq(data([data(:).clu]==u).pdd(t,p1,:)), ...
                    sq(data([data(:).clu]==u).idd(t,p1,:)));

end
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