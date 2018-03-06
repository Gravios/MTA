

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
    tmap{end+1} = tsne([FSr(:,1:6)],cs(sci(1:2:end),:),2,2,perps);
end

cpnts =  ClusterPP;