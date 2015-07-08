Trial = MTATrial('jg05-20120317');

xyz = Trial.xyz.copy;
xyz.load(Trial);

[rhm,fs] = fet_rhm(Trial,[],'wspectral');
% $$$ rhmn = 1./(log10(rhm.data)./repmat(nanmean(clip(log10(rhm.data),-10,6)),[rhm.size(1),1]));
% $$$ rhmn = rhmn./repmat(nanmean(rhmn,2),[1,rhm.size(2)]);

rhmp = max(log10(rhm.data(:,fs>5&fs<13)),[],2);
rhmp = Filter0(gtwin(.2,rhm.sampleRate),rhmp);

%rhmp = mean(rhmn(:,15:30),2);
%rhmp = max(rhmn(:,15:30),[],2);
%rhmp = hold on,plot(log10(max(rhm.data(:,15:30),[],2)),'r')
%rhmp = log10(rhmp);
%

bound_lims = [-500,500;-500,500;-.15,.15];

xyr = MTADxyz('data',[xyz(:,7,1),xyz(:,7,2),rhmp],'sampleRate',xyz.sampleRate);
n
pfs = MTAApfs(Trial,[],'theta',true,'bhvfettest41',[30,30,.2],[1.2,1.2,2.5],'type','xyr','xyz',xyr,'bound_lims',bound_lims);
pfs.parameters.type = 'xyz';

xyz.filter(gtwin(.25,xyz.sampleRate));
vel = xyz.vel(7);
bound_limv = [-500,500;-500,500;-2,2];
xyv = MTADxyz('data',[xyz(:,7,1),xyz(:,7,2),log10(vel.data)],'sampleRate',xyz.sampleRate);

pfv = MTAApfs(Trial,[],'theta',true,'bhvfettest42',[30,30,.2],[1.2,1.2,2.5],'type','xyr','xyz',xyv,'bound_lims',bound_limv);
pfv.parameters.type = 'xyz';



cs = [];
figH = figure(32884),
%for i =1:95,


i =21;
clf
sp(1) = subplot(131);
spfs = reshape(pfs.data.rateMap(:,i),pfs.adata.binSizes');
spfs = spfs(:,:,20:30);
s = slice(spfs,[20],[16],[1:2:10]);set(s,'EdgeColor','none')
sp(2) = subplot(132);
spfv = reshape(pfv.data.rateMap(:,i),pfv.adata.binSizes');
spfv = spfv(:,:,10:20);
s = slice(spfv,[20],[16],[1:2:10]);set(s,'EdgeColor','none')
linkprop(sp,{'cameraposition','cameraupvector'});
subplot(133),plot(spfs(:),spfv(:),'.')
set(figH,'position',[159,351,1380,420])
p_rs = ranksum(spfs(:),spfv(:));
[c_sp,c_p] = corr([spfs(nniz(spfs(:))&nniz(spfv(:))),spfv(nniz(spfs(:))&nniz(spfv(:)))],'type','spearman');
disp(['ranksum p: ' num2str(p_rs)])
disp(['spearman corr: ' num2str(c_sp(2)),'  spearman pval: ' num2str(c_p(2))])
cs(i,:) = [c_sp(2),c_p(2)];
pause(.2)
%end




xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(1.25,xyz.sampleRate));

tag = 'bhvfettest48';
bvel = xyz.vel(1,[1,2]);
bound_limc = [-6,-2;-.5,1.5];
sv = MTADxyz('data',[rhmp,log10(bvel.data)],'sampleRate',xyz.sampleRate);
pfc = MTAApfs(Trial,[],'theta',true,tag,[.2,.1],[1.2,1.2,2.5],'type','xy','xyz',sv,'bound_lims',bound_limc);
pfc.parameters.type = 'xy';


hvel = xyz.vel(7,[1,2]);
bound_limc = [-6,-2;-.5,1.5];
sv = MTADxyz('data',[rhmp,log10(hvel.data)],'sampleRate',xyz.sampleRate);
pfh = MTAApfs(Trial,[],'theta',true,tag,[.2,.1],[1.2,1.2,2.5],'type','xy','xyz',sv,'bound_lims',bound_limc);
pfh.parameters.type = 'xy';

figure,
for i= 1:95
subplot(121),pfc.plot(i,[],1);,title(num2str(i))
subplot(122);pfh.plot(i,[],1);
waitforbuttonpress
end






xyz.load(Trial);
xyz.filter(gtwin(1.25,xyz.sampleRate));


vel = xyz.vel(7);
bound_limx = [-500,500;-6,-2.5;-.5,2];
xsv = MTADxyz('data',[xyz(:,7,1),rhmp,log10(vel.data)],'sampleRate',xyz.sampleRate);
pfsv = MTAApfs(Trial,[],'theta',true,'bhvfettest54',[30,.2,.1],[1.2,1.8,1.8],'type','xyr','xyz',xsv,'bound_lims',bound_limx);
pfsv.parameters.type = 'xyz';

xy = MTADxyz('data',[xyz(:,7,1),xyz(:,7,2),xyz(:,7,3)],'sampleRate',xyz.sampleRate);
pft = MTAApfs(Trial,[],'theta',false,[],[30,30,30],[1.2,1.2,1.2],'type','xy','xyz',xy,'bound_lims',[-500,500;-500,500;0,300]);
pft.parameters.type = 'xyz';



i = 14;
spfs = reshape(pfsv.data.rateMap(:,i),pfsv.adata.binSizes');
sp(1) = subplot(121);
s = slice(spfs,[15],[23],[1:3:25]);set(s,'EdgeColor','none')
sp(2) = subplot(122);
cla
pft.plot(i);
%linkprop(sp,{'cameraposition','cameraupvector'});



xyz.load(Trial);
xyz.filter(gtwin(1.25,xyz.sampleRate));


vel = xyz.vel(7);
bound_limx = [-500,500;-6,-2.5;-.5,2];
xsv = MTADxyz('data',[xyz(:,7,1),rhmp,log10(vel.data)],'sampleRate',xyz.sampleRate);
pfsv = MTAApfs(Trial,[],'theta',false,'bhvfettest57',[30,.2,.1],[1.2,1.8,1.8],'type','xyr','xyz',xsv,'bound_lims',bound_limx);
pfsv.parameters.type = 'xyz';

vel = xyz.vel(1);
bound_limx = [-500,500;-6,-2.5;-.5,2];
xsv = MTADxyz('data',[xyz(:,7,1),rhmp,log10(vel.data)],'sampleRate',xyz.sampleRate);
pfsvb = MTAApfs(Trial,[],'theta',false,'bhvfettest56',[30,.2,.1],[1.2,1.8,1.8],'type','xyr','xyz',xsv,'bound_lims',bound_limx);
pfsvb.parameters.type = 'xyz';

xy = MTADxyz('data',[xyz(:,7,1),xyz(:,7,2),xyz(:,7,3)],'sampleRate',xyz.sampleRate);
pft = MTAApfs(Trial,[],'theta',false,[],[30,30,30],[1.2,1.2,1.2],'type','xy','xyz',xy,'bound_lims',[-500,500;-500,500;0,300]);
pft.parameters.type = 'xyz';


figure,
i = 23;
sp(1) = subplot(131);
spfs = reshape(pfsv.data.rateMap(:,i),pfsv.adata.binSizes');
s = slice(spfs,[15],[23],[1:2:25]);set(s,'EdgeColor','none')
colorbar
sp(2) = subplot(132);
spfsb = reshape(pfsvb.data.rateMap(:,i),pfsvb.adata.binSizes');
s = slice(spfsb,[15],[23],[1:2:25]);set(s,'EdgeColor','none')
colorbar
sp(3) = subplot(133);
cla
pft.plot(i);
linkprop(sp(1:2),{'cameraposition','cameraupvector'});
