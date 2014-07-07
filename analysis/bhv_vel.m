
Trial = MTATrial('jg05-20120310');
xyz = Trial.xyz.copy;
xyz.filter
vel = Trial.transformOrigin;
%figure,plot(hang.roll);
%figure,hist(hang.roll,1000);

gr = nniz(hang.roll);
wrol = nan(size(hang.roll));
wrol(gr) = WhitenSignal(hang.roll(gr));


[ys,fs,ts,phi,fsta] = mtchglong(wrol,2^8,xyz.sampleRate,2^7,2^7*.875,[],'linear',[],[1,20]);
ts = ts+(2^6)/xyz.sampleRate;
ssr = 1/diff(ts(1:2));
pad = round([ts(1),2^7/xyz.sampleRate].*ssr);
szy = size(ys);
ys = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),...
             'sampleRate',ssr);
ts = cat(1,zeros([pad(1),1]),ts,zeros([pad(2),1]));
         

s = 'r';
fhscore = zeros(ys.size(2));
for i = 1:ys.size(2),
    ystf = log10(ys(:,i)); ystf(~nniz(ystf))=[];
    yssf = log10(ys(Trial.stc{s},i)); yssf(~nniz(yssf))=[];
    bins = linspace(prctile(ystf,1),prctile(ystf,99),32);
    Nt = histc(ystf,bins); Ht=-sum(Nt(find(Nt))./sum(Nt).*log2(Nt(find(Nt))./sum(Nt)))+log2(abs(diff(bins(1:2))));
    Ns = histc(yssf,bins); Hs=-sum(Ns(find(Ns))./sum(Ns).*log2(Ns(find(Ns))./sum(Ns)))+log2(abs(diff(bins(1:2))));
    fhscore(i) = (Ht-Hs)/Ht;
end



gfet = median(ys(:,6<fs&fs<14),2);
gsfet = (median(ys(:,8<fs&fs<14),2)-median(ys(:,8>fs|14<fs),2))./(median(ys(:,8<fs&fs<14),2)+median(ys(:,8>fs|14<fs),2));
figure,
sp(1) = subplot(311); imagesc(ys.data'),axis xy
sp(2) = subplot(312); plot(gfet);
sp(3) = subplot(313); plot(gsfet);
linkaxes(sp,'x');


figure,
sp(1) = subplot(311); imagesc(log10(ys.data')),axis xy
ys.data = ys.data./repmat(sum(ys.data,2),[1,ys.size(2)]);
ys.data = log10(ys.data);
sp(2) = subplot(312); imagesc(ys.data'),axis xy
ys.data(nniz(ys),:) = ys(nniz(ys),:)-repmat(mean(ys(nniz(ys),:)),[sum(nniz(ys.data)),1]);
sp(3) = subplot(313); imagesc(ys.data'),axis xy
linkaxes(sp,'x');



[U,D,V] = svd(cov(ys(nniz(ys),:)));
figure,imagesc(fs,1:ys.size(2),V')

gsc = [];
for i = 1:5,
gsc(:,i) = sum(ys.data.*repmat(V(:,i)',ys.size(1),1),2);
end

figure,hist2([max(gsc(:,1:2),[],2),gsc(:,5)],100,100)
caxis([0,100])


[A,W]=fastica(ys(nniz(ys),:)','g','gauss','approach','symm','displayMode','off');
