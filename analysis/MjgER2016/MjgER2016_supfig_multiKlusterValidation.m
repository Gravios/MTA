

Trial = MTATrial.validate('ER06-20130612.cof.all');
xyz = preproc_xyz(Trial,'trb');
tper = Trial.stc{'t',xyz.sampleRate};

cd(Trial.spath);
[resOld,cluOld,mapOld] = LoadCluRes('ER06-20130612',[1,2]);
[resOld,rind] = SelectPeriods(resOld,round(Trial.sync.data.*Trial.sampleRate),'d',1,0);
cluOld = cluOld(rind);
resOld = round(resOld./Trial.sampleRate.*xyz.sampleRate);
resOld = resOld-round(Trial.sync(1).*xyz.sampleRate);
cluOld(resOld<1|resOld>size(xyz,1)) = [];
resOld(resOld<1|resOld>size(xyz,1)) = [];
[resOld,rind] = SelectPeriods(resOld,tper.data,'d',1,0);
cluOld = cluOld(rind);


cd(fullfile(Trial.spath,'klustering'));
[resNew,cluNew,mapNew] = LoadCluRes('ER06-20130612',[1,2]);
[resNew,rind] = SelectPeriods(resNew,round(Trial.sync.data.*Trial.sampleRate),'d',1,0);
cluNew = cluNew(rind);
resNew = round(resNew./Trial.sampleRate.*xyz.sampleRate);
resNew = resNew-round(Trial.sync(1).*xyz.sampleRate);
cluNew(resNew<1|resNew>size(xyz,1)) = [];
resNew(resNew<1|resNew>size(xyz,1)) = [];
[resNew,rind] = SelectPeriods(resNew,tper.data,'d',1,0);
cluNew = cluNew(rind);


mccg = zeros([numel(unique(cluOld)),numel(unique(cluNew)),2,2]);
for o = unique(cluOld)',
    for n = unique(cluNew)',
        uccg = CCG([resOld(cluOld==o);resNew(cluNew==n)],... 
                   [ones([sum(cluOld==o),1]);2*ones([sum(cluNew==n),1])],...
                   32,6,32000,[1,2],'count');        
        mccg(o,n,:,:) = permute(uccg(7,:,:),[4,1,2,3]);
    end
end


figure();
imagesc(log10(mccg(:,:,1,2)'+1));
axis('xy');


xbinInds = discretize(xyz(:,'hcom',1),linspace(-500,500,21));
ybinInds = discretize(xyz(:,'hcom',2),linspace(-500,500,21));

xyocc = accumarray([SelectPeriods(xbinInds,tper.data,'c',1,1),SelectPeriods(ybinInds,tper.data,'c',1,1)],1,[20,20],@sum);

figure,
for o = 1:92;
[~,nind] = max(mccg(o,:,1,2));
clf();
subplot(121);
imagesc(linspace(-500,500,21),...
        linspace(-500,500,21),...
        (accumarray([xbinInds(resOld(cluOld==o)),ybinInds(resOld(cluOld==o))],1,[20,20],@sum)./xyocc.*xyz.sampleRate)');
axis('xy');
caxis([0,8]);
title(num2str(o));
subplot(122);
imagesc(linspace(-500,500,21),...
        linspace(-500,500,21),...
        (accumarray([xbinInds(resNew(cluNew==nind)),ybinInds(resNew(cluNew==nind))],1,[20,20],@sum)./xyocc.*xyz.sampleRate)');
axis('xy');
caxis([0,8]);
title(num2str(nind));
waitforbuttonpress();
end