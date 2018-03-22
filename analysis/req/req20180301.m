

% Nope, shouldn't be touching this .... 

Trial = MTATrial.validate('Ed10-20140816.cof.all');

ncp = Trial.load('lfp',66);
ofb = Trial.load('lfp',33:64);

ncp.filter('RectFilter',5,5);


[mins,mval] = LocalMinima(ncp.data,100,-1000);

%scycles = zeros([numel(mins)-1,1000]);
for index = 1:numel(mins)-1,
    scycles(index,:) = interp1(1:size(ncp,1),ncp.data,linspace(mins(index),mins(index+1),1000));
end

save('/storage/gravio/data/project/general/analysis/req20180301.mat','scycles');
load('/storage/gravio/data/project/general/analysis/req20180301.mat');

figure,imagesc(scycles)
ucycles = unity(scycles')';
figure,imagesc(ucycles)

tucycles = ucycles;

%[LU,LR,FSr,VT] = erpPCA(scycles(1:2:end,:),6);
[LU,LR,FSr,VT] = erpPCA(ucycles(1:2:end,:),6);


figure,plot(VT(:,4))
figure,plot(LR(:,1:6))

figure,plot3(FSr(:,1),FSr(:,2),FSr(:,3),'.');

figure,hist(,linspace(0,30,30))
figure,plot(mean([circshift(mins,-1)-mins,mins-circshift(mins,1)],2)./1250)
sfreq = clip(1250./mean([circshift(mins,-1)-mins,mins-circshift(mins,1)],2),0.01,100);

figure,hist(sfreq,linspace(0,30,30))

eds = linspace(0,15,30);
sci = discretize(sfreq,eds);
cs = jet(30);

%tmap = tsne([FSr(1:2:end,1:6),sfreq(1:2:index)],cs(sci(1:2:end),:),2,3,80);
perps = [25,50,75,100];
tmap = {};
for p=perps,
    figure,
    tmap{end+1} = tsne([FSr(:,1:6)],cs(sci(1:2:end),:),2,4,p);
end

hfig = gcf();

hfig = figure();
plot(tmap{3}(:,1),tmap{3}(:,2),'.','MarkerSize',1)
%scatter(tmap{3}(:,1),tmap{3}(:,2),5,cs(sci(1:2:end),:));
cpnts =  ClusterPP(hfig);

whos cpnts
whos ucycles
tucycles = ucycles(1:2:end,:);
figure();
hold('on');
cids = unique(cpnts);
cluColors = parula(numel(cids));
for t = cids
    plot(mean(tucycles(cpnts==t,:)),'Color',cluColors(t+1,:))
end

figure,
for t = cids,   
subplot(2,6,t+1),
plot(tucycles(cpnts==t,:)')
end



defspec = struct('nFFT',2^8,'Fs',ofb.sampleRate,...
                 'WinLength',2^7,'nOverlap',2^7*.875,...
                 'FreqRange',[20,150]);

for c = 1:32,
    tofb = ofb.copy('empty');
    tofb.data = ofb.data(:,c);
    [tys,fs,ts] = fet_spec(Trial,tofb,'mtchglong',true,[],defspec);
    if c==1, 
        ys = tys;
    else
        ys.data = cat(3,ys.data,tys.data);
    end
end

tmpMins = mins(1:2:end);

nys = ys.copy;
nys.data = nunity(nys.data);

tt = [-40:40]./ys.sampleRate;


t = 6;
ysegs = ys.segs(round(tmpMins(cpnts==t)./1250*ys.sampleRate)'-40,80);

%ind = {':',':',18,':'};
ind = {':',':',':',2};
%scaleFun = @(x) log10(x);
scaleFun = @(x) x;

figure,
subplot(311);
imagesc(tt,fs,sq(mean(scaleFun(ysegs(ind{:})),2,'omitnan'))'),axis('xy');
subplot(312);
imagesc(tt,fs,[sq(mean(scaleFun(ysegs(ind{:})),2,'omitnan'))...
               ./sq(std(scaleFun(ysegs(ind{:})),[],2,'omitnan'))]'),axis('xy');
subplot(313);
plot(mean(tucycles(cpnts==t,:)),'Color',cluColors(t+1,:))


tsfreq = sfreq(1:2:end);
[~,sfi] = sort(tsfreq(cpnts==t));



figure,imagesc(ysegs(:,sfi,18,2)')
%caxis([-1,3])

nofb = ofb.copy();
nofb.data = nunity(nofb);

osegs = nofb.segs(tmpMins(cpnts==t)-250,500);

tsfreq = sfreq(1:2:end);
[~,sfi] = sort(tsfreq(cpnts==t));

figure,imagesc(osegs(:,sfi,2)');

