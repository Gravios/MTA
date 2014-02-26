

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317');


v = Trial.vel;

v = Filter0(gausswin(21)./sum(gausswin(21)),v);




c =  xcorr(log10(v(v(:,1)~=0,1)),log10(v(v(:,1)~=0,1)),10000);
figure,plot(c./max(c))

hist2([log10(v(v(:,1)~=0,1)),circshift(log10(v(v(:,1)~=0,1)),800)],100,100);

