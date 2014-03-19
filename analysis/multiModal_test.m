
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317');
Trial.load('xyz')
Trial.xyz.filter(gausswin(61)./sum(gausswin(61)));

marker = 7;
smed = 'descend';
smed = 'ascend';

ws = 100;

xyzparts = GetSegs(sort(log10(Trial.xyz(Trial.xyz(:,marker,3)>0,marker,3)),smed),1:25:sum(Trial.xyz(:,marker,3)>0),ws,1);


if matlabpool('size')==0,matlabpool open 12,end

bns = round(linspace(1,size(xyzparts,2)-100,500));
mss = zeros(500,1);
sss = zeros(500,1);
rind = randi(100,100,10000);
ms = zeros(size(rind,2),1);

for s = 1:500,
parfor i = 1:size(rind,2)
xp = xyzparts(:,bns(s)+rind(:,i));
ms(i) = mean(xp(:));
end
mss(s) = mean(ms);
sss(s) = std(ms);
end

%figure,plot(ms)
% $$$ figure,plot(mean(xyzparts(:,1:end-1))',Filter0(gausswin(61)./sum(gausswin(61)),diff(Filter0(gausswin(61)./sum(gausswin(61)),ms))))
% $$$ 
% $$$ figure,
% $$$ hold on
% $$$ plot(mean(xyzparts(:,1:end-2))',-diff(Filter0(gausswin(61)./sum(gausswin(61)),diff(Filter0(gausswin(61)./sum(gausswin(61)),ms)))),'b')
% $$$ Lines([],0,'k')
figure,hist(xyzparts(:),1000)
figure,hist(ms(:),1000)

d=2;
figure,plot(mean(xyzparts(:,1:end-d))',abs(diff(ms(:,:),d)'))
hold on,plot(mean(xyzparts(:,1:end-d))',abs(Filter0(gausswin(61)./sum(gausswin(61)),diff(ms(:,:),d))'),'m')


figure,plot(xyzparts(100,bns(1:end)),(Filter0(gausswin(5)./sum(gausswin(5)),mss)))
hold on,plot(xyzparts(100,bns(1:end)),(Filter0(gausswin(5)./sum(gausswin(5)),mss+sss)),'r')

figure,plot(xyzparts(100,bns(1:end-2)),diff(Filter0(gausswin(5)./sum(gausswin(5)),mss),2))

figure,
subplot(211),plot(xyzparts(100,bns(1:end)),-log10(sss),'r')
xlim([.5,2.5])
subplot(212),hist(xyzparts(:),1000)
xlim([.5,2.5])


hold on,plot(xyzparts(100,bns(1:end)),diff(Filter0(gausswin(5)./sum(gausswin(5)),mss)),'r')

