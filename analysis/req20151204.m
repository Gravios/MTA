%tsne mapping of spike on to 2D plane
cd /storage/gravio/data/project/general/Ed10-20140820/

%Spk = LoadSpk('Ed10-20140820.fet.1');
[Clu,nclu] = LoadClu('Ed10-20140820.clu.1');
Clu = Clu(ismember(Clu,1:36));
[Fet,nfet] = LoadFet('Ed10-20140820.fet.1',[],1:36);

Fet = Fet(1:20:end,mod(1:size(Fet,2),3)~=0);
Fet = Fet(:,1:end-1);
Clu = Clu(1:20:end);

c = jet(36);
msc = c(Clu,:);
mx = tsne(Fet,msc,2,8,100);

figure(1)
xlm = xlim(gca);
ylm = ylim(gca);


hfig = figure(2); hold on
plot(mx(Clu~=1,1),mx(Clu~=1,2),'.b');
plot(mx(Clu==1,1),mx(Clu==1,2),'.r');
xlim(xlm);
ylim(ylm);