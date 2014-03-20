MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
Trial = MTATrial('jg05-20120317');
fwin = gausswin(11)./sum(gausswin(11));
Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.xyz.filter(fwin);

angvel = [];
for i = 6:6:60,
    angvel(:,end+1) = circshift(circ_dist(Trial.ang(:,5,7,1),circshift(Trial.ang(:,5,7,1),i)),0);%,round(i/2));
end

figure,plot(angvel)
hold on
Lines(Trial.stc{'r'}.data(:),[],'r')
Lines(Trial.stc{'w'}.data(:),[],'k')

%figure,hist(angvel(:,2),100)
angvel(isinf(angvel))=0;
angvel(isnan(angvel))=0;

wav = WhitenSignal(angvel(:,2:2:10),[],1);

[ya,fa,ta] = mtchglong(wav,2^9,Trial.ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);

figure
ns = size(ya,3);
sp = [];
for i = 1:ns
sp(i) = subplot(ns,1,i);
imagesc(ta,fa,log10(ya(:,:,i,i))'),axis xy
end
linkaxes(sp,'xy');

yasw = GetSegs(log10(ya(:,:,5,5)),Trial.stc{'w',1/(diff(ta(1:2))/2)}.data(:,1)-10,21,0);
yasw(isinf(yasw))=nan;

figure(12332), imagesc(1:21,fa,sq(nanmean(yasw,2))'),axis xy
figure(12333), imagesc(1:21,fa,sq(nanstd(yasw,[],2))'),axis xy


avm = [max(angvel,[],2),min(angvel,[],2)];

figure,plot(diff(avm,1,2));
Lines(Trial.stc{'r'}.data(:),[],'r')
Lines(Trial.stc{'w'}.data(:),[],'k')
