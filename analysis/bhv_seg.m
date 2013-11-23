
t = MTATrial('jg05-20120317');

windows = [9,27,63,81,121];
dwindows = round([0,diff(windows)]);
dims = [1,2];
dsf = 4;

sl = t.xyz.size(1);
sl = sl-mod(sl,dsf);
fx = zeros(sl/dsf,t.xyz.model.N,numel(windows));
for n = 1:numel(windows)
    tfx = zeros(windows(n),sl/dsf,t.xyz.model.N,numel(dims));
    for m =1:t.xyz.model.N;
        tfx(:,:,m,:) = GetSegs(Filter0(gausswin(9)./sum(gausswin(9)),sq(t.xyz(1:sl,m,dims))),1:dsf:sl,windows(n),0);
    end
    fx(:,:,n) = sq(median(sqrt(sum((tfx-repmat(tfx(1,:,:,:),windows(n),1)).^2,4))));
end

dims = [1,2];
sl = t.xyz.size(1);
sl = sl-mod(sl,dsf);
fcv = zeros(4,4,sl/dsf,numel(windows));
for n = 1:numel(windows)
    tfcv = zeros(windows(n),sl/dsf,4,numel(dims));
    for m =1:4;
        tfcv(:,:,m,:) = GetSegs(Filter0(gausswin(9)./sum(gausswin(9)),sq(t.xyz(1:sl,m,dims))),1:dsf:sl,windows(n),0);
    end
    ttfcv = sqrt(sum((tfcv-repmat(tfcv(1,:,:,:),windows(n),1)).^2,4));
    for i = 1:sl/dsf,
        fcv(:,:,i,n) = cov(sq(ttfcv(:,i,:)));
    end
end

dfcv = abs(diff(fcv,1,2));
tsc = diff(sq(std(dfcv)),2,1);
score = sq(abs(tsc));

rx = zeros(sl/dsf,t.xyz.model.N,numel(windows));
for n = 1:numel(windows)
    trx = zeros(windows(n),sl/dsf,t.xyz.model.N,numel(dims));
    for m =1:t.xyz.model.N;
        trx(:,:,m,:) = GetSegs(Filter0(gausswin(9)./sum(gausswin(9)),sq(t.xyz(1:sl,m,dims))),1:dsf:sl,windows(n),0);
    end
    rx(:,:,n) = sq(median(sqrt(sum((trx-repmat(trx(windows(n),:,:,:),windows(n),1)).^2,4))));
end


dx = zeros(sl/dsf,t.xyz.model.N,numel(windows));
for n = 1:numel(windows)
    tdx = zeros(windows(n),sl/dsf,t.xyz.model.N,numel(dims));
    for m =1:t.xyz.model.N;
        tdx(:,:,m,:) = GetSegs(Filter0(gausswin(9)./sum(gausswin(9)),sq(t.xyz(1:sl,m,dims))),1:dsf:sl,windows(n),0);
    end
    dx(:,:,n) = sq(sum(sqrt(sum(diff(tdx).^2,4))));
end

mfx = MTADxyz([],[],fx,t.xyz.sampleRate/dsf);
mdmfx = MTADxyz([],[],diff(mfx.data,1,3),t.xyz.sampleRate/dsf);

mrx = MTADxyz([],[],rx,t.xyz.sampleRate/dsf);
mdmrx = MTADxyz([],[],diff(mrx.data,1,3),t.xyz.sampleRate/dsf);

mdx = MTADxyz([],[],dx,t.xyz.sampleRate/dsf);
mdmdx = MTADxyz([],[],diff(mdx.data,1,3),t.xyz.sampleRate/dsf);



t.stc{'r',mdmdx.sampleRate};
t.stc{'w',mdmdx.sampleRate};
t.stc{'t',mdmdx.sampleRate};

m = 7;
w = 3;
figure,hist2([log10(abs(mfx(t.stc{'w'}.data,1,w)+.01)),log10(abs(mfx(t.stc{'w'}.data,m,w)+.01))],100,100)
Lines([],[1,2],'w');line([0,2],[0,2],'color',[1,1,1]);
figure,hist2([log10(abs(mfx(t.stc{'r'}.data,1,w)+.01)),log10(abs(mfx(t.stc{'r'}.data,m,w)+.01))],100,100)
Lines([],[1,2],'w');line([0,2],[0,2],'color',[1,1,1]);
figure,hist2([log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),1,w)+.01)),log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+.01))],100,100)
Lines([],[1,2],'w');line([0,2],[0,2],'color',[1,1,1]);

figure,hist2([log10(abs(fx(:,m,w)+.01)),log10(abs(dx(:,m,w)+.01))],100,100)

figure,hist2([log10(abs(mfx(t.stc{'w'},1,w)+.1)),log10(abs(mdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(t.stc{'r'},1,w)+.1)),log10(abs(mdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),1,w)+.1)),log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])

figure,hist2([log10(abs(mfx(t.stc{'w'},m,w)+.1)),log10(abs(mdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(t.stc{'r'},m,w)+.1)),log10(abs(mdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+.1)),log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])

o = 1;
figure,hist2([log10(abs(mfx(t.stc{'w'},o,w)+.1)),log10(abs(mdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(t.stc{'r'},o,w)+.1)),log10(abs(mdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])

figure,hist2([log10(abs(mfx(:,o,w)+.1)),log10(abs(mdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])


m = 4;
o = 7;
figure,hist2([log10(abs(mfx(t.stc{'w'},o,w)+.1)),log10(abs(mdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(t.stc{'r'},o,w)+.1)),log10(abs(mdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mfx(:,o,w)+.1)),log10(abs(mdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])



m = 1;
o = 7;
w = 2;
figure,hist2([log10(abs(mdx(t.stc{'w'},o,w)+.1)),log10(abs(mdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdx(t.stc{'r'},o,w)+.1)),log10(abs(mdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdx(:,o,w)+.1)),log10(abs(mdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])


m = 1;
o = 7;
w = 2;
figure,hist2([log10(abs(mdmdx(t.stc{'w'},o,w)+.1)),log10(abs(mdmdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmdx(t.stc{'r'},o,w)+.1)),log10(abs(mdmdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdmdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,100]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmdx(:,o,w)+.1)),log10(abs(mdmdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])

m = 1;
o = 7;
figure,hist2([log10(abs(mdmfx(t.stc{'w'},o,w)+.1)),log10(abs(mdmdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmfx(t.stc{'r'},o,w)+.1)),log10(abs(mdmdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdmdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
figure,hist2([log10(abs(mdmfx(:,o,w)+.1)),log10(abs(mdmdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])


m = 2;
o = 7;
figure,
for w = 1:4,
subplot2(4,4,1,w),hist2([log10(abs(mdmfx(t.stc{'w'},o,w)+.1)),log10(abs(mdmdx(t.stc{'w'},m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
subplot2(4,4,2,w),hist2([log10(abs(mdmfx(t.stc{'r'},o,w)+.1)),log10(abs(mdmdx(t.stc{'r'},m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
subplot2(4,4,3,w),hist2([log10(abs(mdmfx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),o,w)+.1)),log10(abs(mdmdx(cat(1,t.stc{'walk'}.data,t.stc{'rear'}.data),m,w)+1))],100,100),caxis([0,50]),line([0,2],[0,2],'color',[1,1,1])
subplot2(4,4,4,w),hist2([log10(abs(mdmfx(:,o,w)+.1)),log10(abs(mdmdx(:,m,w)+1))],100,100),caxis([0,500]),line([0,2],[0,2],'color',[1,1,1])
end


ind = log10(mdx(:,7,1)+1)>.5&log10(mdx(:,1,1)+1)<1.5;
%ind = t.stc{'r'};
%ind = cat(1,t.stc{'w'}.data,t.stc{'r'}.data);
figure,
for i = 2:5,
    subplot(1,4,i-1),
    hist2(log10([mdx(ind,1,1),mdx(ind,7,i)]+1),100,100),
    caxis([0,20]),
    line([0,4],[0,4],'color',[1,1,1]),
end


mlims = [-.25,2.5;0,2.5;2.5,2.5;-.25,0];

w=2;
m=1;
o=7;
figure,
subplot(221),hist2(cat(1,[clip(log10(abs(mdmfx(:,m,w)+1)),-0.25,5),log10(abs(mdmdx(:,o,w)+1))],mlims),100,100);caxis([0,150])
subplot(222),hist2(cat(1,[clip(log10(abs(mdmfx(t.stc{'r'},m,w)+1)),-0.25,5),log10(abs(mdmdx(t.stc{'r'},o,w)+1))],mlims),100,100);caxis([0,50])
subplot(223),hist2(cat(1,[clip(log10(abs(mdmfx(t.stc{'w'},m,w)+1)),-0.25,5),log10(abs(mdmdx(t.stc{'w'},o,w)+1))],mlims),100,100);caxis([0,50])


figure,
subplot(221),hist2(cat(1,[clip(log10(abs(mfx(:,m,w)+1)),-0.25,5),log10(abs(mdx(:,o,w)+1))],mlims),100,100);caxis([0,150])
subplot(222),hist2(cat(1,[clip(log10(abs(mfx(t.stc{'r'},m,w)+1)),-0.25,5),log10(abs(mdx(t.stc{'r'},o,w)+1))],mlims),100,100);caxis([0,50])
subplot(223),hist2(cat(1,[clip(log10(abs(mfx(t.stc{'w'},m,w)+1)),-0.25,5),log10(abs(mdx(t.stc{'w'},o,w)+1))],mlims),100,100);caxis([0,50])


figure,
subplot(221),hist2(cat(1,[clip(log10(abs(mfx(:,m,w)+1)),-0.25,5),log10(abs(mfx(:,o,w)+1))],mlims),100,100);caxis([0,150])
subplot(222),hist2(cat(1,[clip(log10(abs(mfx(t.stc{'r'},m,w)+1)),-0.25,5),log10(abs(mfx(t.stc{'r'},o,w)+1))],mlims),100,100);caxis([0,50])
subplot(223),hist2(cat(1,[clip(log10(abs(mfx(t.stc{'w'},m,w)+1)),-0.25,5),log10(abs(mfx(t.stc{'w'},o,w)+1))],mlims),100,100);caxis([0,50])


figure,
subplot(221),hist2(cat(1,[clip(log10(abs(mdx(:,m,w)+1)),-0.25,5),log10(abs(mdx(:,o,w)+1))],mlims),100,100);caxis([0,150])
subplot(222),hist2(cat(1,[clip(log10(abs(mdx(t.stc{'r'},m,w)+1)),-0.25,5),log10(abs(mdx(t.stc{'r'},o,w)+1))],mlims),100,100);caxis([0,50])
subplot(223),hist2(cat(1,[clip(log10(abs(mdx(t.stc{'w'},m,w)+1)),-0.25,5),log10(abs(mdx(t.stc{'w'},o,w)+1))],mlims),100,100);caxis([0,50])




mz = MTADxyz([],[],t.xyz(round(linspace(1,t.xyz.size(1),mdx.size(1))),7,3),t.xyz.sampleRate/dsf);
figure,hold on
plot3(log10(abs(mdx(:,m,w)+1)),log10(abs(mdx(:,o,w)+1)),mz.data,'.k')
plot3(log10(abs(mdx(t.stc{'w'},m,w)+1)),log10(abs(mdx(t.stc{'w'},o,w)+1)),mz(t.stc{'w'}),'.')
plot3(log10(abs(mdx(t.stc{'r'},m,w)+1)),log10(abs(mdx(t.stc{'r'},o,w)+1)),mz(t.stc{'r'}),'.r')




mz = MTADxyz([],[],t.xyz(round(linspace(1,t.xyz.size(1),mdx.size(1))),7,3),t.xyz.sampleRate/dsf);
figure,hold on
%plot3(log10(abs(mfx(:,m,w)+1)),log10(abs(mfx(:,o,w)+1)),mz.data,'.k')
plot3(log10(abs(mfx(t.stc{'w'},m,w)+1)),log10(abs(mfx(t.stc{'w'},o,w)+1)),mz(t.stc{'w'}),'.')
plot3(log10(abs(mfx(t.stc{'r'},m,w)+1)),log10(abs(mfx(t.stc{'r'},o,w)+1)),mz(t.stc{'r'}),'.r')







figure,hist2(cat(1,[clip(log10(abs(mdmfx(:,m,w)+1)),-0.25,5),log10(abs(mdmdx(:,o,w)+1))],mlims),100,100);caxis([0,150])
figure,hist2([log10(abs(mdmdx(:,m,w)+1)),log10(abs(mdmdx(:,o,w)+1))],100,100);caxis([0,150])
figure,hist2([log10(abs(mdmfx(:,m,w)+1)),log10(abs(mdmfx(:,o,w)+1))],100,100);caxis([0,150])



figure,hist2([log10(abs(mdmfx(:,m,w)+1)),log10(abs(mdmfx(:,o,w)+1))],100,100);caxis([0,150])
figure,hist2([log10(abs(mfx(:,m,w)+1)),log10(abs(mfx(:,o,w)+1))],100,100);caxis([0,150])
figure,hist2([log10(abs(mrx(:,m,w)+1)),log10(abs(mrx(:,o,w)+1))],100,100);caxis([0,150])
figure,hist2([log10(abs(mrx(:,m,w)+1)),log10(abs(mfx(:,o,w)+1))],100,100);caxis([0,150])


figure,hist2([log10(abs(mfx(:,o,w)+1)),log10(abs(mrx(:,o,w)+1))],100,100);caxis([0,150])
w = 2;
i=1;
figure
subplot(3,3,i),hist2([log10(abs(mfx(:,m,w)+1)),log10(abs(mfx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(:,m,w)+1)),log10(abs(mdx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(:,m,w)+1)),log10(abs(mfx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;

subplot(3,3,i),hist2([log10(abs(mfx(t.stc{'w'},m,w)+1)),log10(abs(mfx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(t.stc{'w'},m,w)+1)),log10(abs(mdx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(t.stc{'w'},m,w)+1)),log10(abs(mfx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;

subplot(3,3,i),hist2([log10(abs(mfx(t.stc{'r'},m,w)+1)),log10(abs(mfx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(t.stc{'r'},m,w)+1)),log10(abs(mdx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdx(t.stc{'r'},m,w)+1)),log10(abs(mfx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
ForAllSubplots('xlim([0,2.5]),ylim([0,2.5])')


i=1;
figure
subplot(3,3,i),hist2([log10(abs(mdmfx(:,m,w)+1)),log10(abs(mdmfx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(:,m,w)+1)),log10(abs(mdmdx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(:,m,w)+1)),log10(abs(mdmfx(:,o,w)+1))],100,100);caxis([0,150]),i=i+1;

subplot(3,3,i),hist2([log10(abs(mdmfx(t.stc{'w'},m,w)+1)),log10(abs(mdmfx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(t.stc{'w'},m,w)+1)),log10(abs(mdmdx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(t.stc{'w'},m,w)+1)),log10(abs(mdmfx(t.stc{'w'},o,w)+1))],50,50);caxis([0,50]),i=i+1;

subplot(3,3,i),hist2([log10(abs(mdmfx(t.stc{'r'},m,w)+1)),log10(abs(mdmfx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(t.stc{'r'},m,w)+1)),log10(abs(mdmdx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
subplot(3,3,i),hist2([log10(abs(mdmdx(t.stc{'r'},m,w)+1)),log10(abs(mdmfx(t.stc{'r'},o,w)+1))],50,50);caxis([0,50]),i=i+1;
ForAllSubplots('xlim([-.5,2.5]),ylim([-.5,2.5])')



m = 2;
o = 7;
w = 2;
sp = [];
mind = abs(mdmdx(:,m,w)~=0)&abs(mdmfx(:,m,w)~=0);
%pc = princomp([log10(abs(mdmfx(mind,o,w))),log10(abs(mdmdx(mind,m,w)))]);
[A,W] = fastica([log10(abs(mdmfx(mind,o,w))),log10(abs(mdmdx(mind,m,w)))]');
tcor = [log10(abs(mdmfx(:,o,w)+.1)),log10(abs(mdmdx(:,m,w)+1))]*A';
figure,hist2(tcor(mind,:),100,100),caxis([0,500])

tci = tcor(:,2)>.1598;
tcil = find(~tci);
tcig = find(tci);

figure,plot(mdmfx(:,1,1),'.')
hold on,plot(tcil,mdmfx(~tci,1,1),'.k')

scor = tcor(tci,:);
[SA,SW] = fastica(tcor(tci,:)');
scor = scor*inv(SA');
figure,hist2(scor,100,100),caxis([0,300])

mci = scor(:,1)>2.066;
mcil = find(~mci);
mcig = find(mci);


hcor = scor(mci,:);
[HA,HW] = fastica(hcor');
thcor = hcor*HA';
figure,hist2(thcor,100,100),caxis([0,100])

hci = thcor(:,2)<-2.14;
hcil = find(~hci);
hcig = find(hci);

dmx = [mdmdx(:,1,w),mdmdx(:,7,w)];
sdmx = GetSegs(dmx,1:size(dmx,1),21,0);

rdmx = zeros(size(dmx,1),1);
pdmx = zeros(size(dmx,1),1);
for i = 1:size(dmx),
     [r,p]= corr(sq(sdmx(:,i,:)));
     rdmx(i) = r(1,2);
     pdmx(i) = p(1,2);
end

figure,hold on,
plot(mdmfx(:,1,w))
plot(mdmfx(:,7,w))
plot(mdmrx(:,1,w),'m')
plot(mdmrx(:,7,w),'m')
plot(mdmdx(:,7,w),'g')
plot(mdmdx(:,1,w),'g')
plot(tcil,mdmdx(tcil,1,w),'.k')
plot(tcig(mcil),mdmdx(tcig(mcil),1,w),'.m')
plot(tcig(mcig(hcil)),mdmdx(tcig(mcig(hcil)),1,w),'.g')
plot(tcig(mcig(hcig)),mdmdx(tcig(mcig(hcig)),1,w),'.c')
Lines(t.stc{'w',mdmfx.sampleRate}(:)-diff(windows([w,w+1]))/2,[],'k');
Lines(t.stc{'r',mdmfx.sampleRate}(:)-diff(windows([w,w+1]))/2,[],'r');
plot(log10(abs(log10(abs(mean(mdmdx(:,1:3,w),2).^2-mdmdx(:,7,w).*mean(mdmdx(:,1:3,w),2))+eps).*mean([mean(mdmdx(:,1:3,w),2)],2))+eps)*50,'.y')
plot(cv*35,'.r')
%plot(log10(abs(mean(mdmdx(:,1:2,w),2).^2-mdmdx(:,7,w).*mean(mdmdx(:,1:2,w),2)+eps)),'.r')
Lines([],100,'k')

m = [1,2];
h = [5,6,7,8];
cv = log10(abs(mean(mdmdx(:,m,w),2).^2-mean(mdmdx(:,h,w),2).*mean(mdmdx(:,m,w),2)+eps)) +log10(abs(mean(mdmdx(:,h,w),2)-mean(mdmdx(:,m,w),2))./abs(mean(mdmdx(:,h,w),2)+mean(mdmdx(:,m,w),2))+eps);
plot(cv,'r')


figure,hold on,
plot(mdmfx(:,1,w)+mdmrx(:,1,w),'b')
plot(mdmfx(:,7,w)+mdmrx(:,7,w),'r')
plot(log10(abs(log10(abs(mean(mdmdx(:,1:2,w),2).^2-mdmdx(:,7,w).*mean(mdmdx(:,1:2,w),2))+eps).*mean([mean(mdmdx(:,1:2,w),2)],2))+eps)*50,'.y')
Lines(t.stc{'w',mdmfx.sampleRate}(:)-diff(windows([w,w+1]))/2,[],'k');
Lines(t.stc{'r',mdmfx.sampleRate}(:)-diff(windows([w,w+1]))/2,[],'r');
Lines([],100,'k')



plot([1:size(rdmx,1)],rdmx.*20,'.y')
plot([1:size(rdmx,1)],-log10(pdmx).*5,'.r')

figure,
sp(1) = subplot(3,1,1);plot(log10(abs(mdmfx(:,o,w)+.1)),log10(abs(mdmdx(:,m,w)+1)),'.','MarkerSize',5),line([0,4],[0,4],'color',[0,0,0]);
sp(2) = subplot(3,1,2);plot(log10(abs(mdmfx(t.stc{'w'},o,w)+.1)),log10(abs(mdmdx(t.stc{'w'},m,w)+1)),'.g','MarkerSize',5),line([0,4],[0,4],'color',[0,0,0]);
sp(3) = subplot(3,1,3);plot(log10(abs(mdmfx(t.stc{'r'},o,w)+.1)),log10(abs(mdmdx(t.stc{'r'},m,w)+1)),'.r','MarkerSize',5),line([0,4],[0,4],'color',[0,0,0]);
linkaxes(sp,'xy')


m = 7;
figure, hold on
plot(mfx(:,m,2))
plot(mrx(:,m,2),'r')
plot(mean([mrx(:,m,2),mfx(:,m,2)],2),'g')
plot(diff([mrx(:,m,2),mfx(:,m,2)],1,2))



ds = [diff([mrx(:,1,2),mfx(:,1,2)],1,2),diff([mrx(:,7,2),mfx(:,7,2)],1,2)];
wds = WhitenSignal(ds,[],1);
[yo,fo,to,phi,fst] = mtchglong(wds,2^6,mrx.sampleRate,2^5,2^5-1,[],[],[],[0.1,15]);


figure,imagesc(to,fo,log10(yo(:,:,1,1)')),axis xy,caxis([-2.57,0])
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,2)),[],'r');
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,1)),[],'k');

sp = [];
figure,
sp(1)=subplot(211);imagesc(to,fo,log10(yo(:,:,1,1)')),axis xy,
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,2)),[],'r');
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,1)),[],'k');
sp(2)=subplot(212);imagesc(to,fo,log10(yo(:,:,2,2)')),axis xy,
linkaxes(sp,'xy')
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,2)),[],'r');
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,1)),[],'k');
    

figure,imagescnan({to,fo,phi(:,:,1,2)'},[],1,1),axis xy
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,2)),[],'r');
Lines(to(t.stc{'w',mdmfx.sampleRate}(:,1)),[],'k');

