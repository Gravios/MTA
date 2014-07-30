Trial = MTATrial('jg05-20120317');

xyz = Trial.xyz.copy;
xyz.load(Trial);

rhm = fet_rhm(Trial,[],'wspectral');
rhmn = 1./(log10(rhm.data)./repmat(nanmean(clip(log10(rhm.data),-10,6)),[rhm.size(1),1]));
rhmn = rhmn./repmat(nanmean(rhmn,2),[1,rhm.size(2)]);

rhmp = max(log10(rhm.data(:,15:30)),[],2);
rhmp = mean(rhmn(:,15:30),2);
%rhmp = max(rhmn(:,15:30),[],2);
%rhmp = hold on,plot(log10(max(rhm.data(:,15:30),[],2)),'r')
rhmp = log10(rhmp);
rhmp = Filter0(gausswin(31)./sum(gausswin(31)),rhmp);

bound_lims = [-500,500;-500,500;-.15,.15];
bound_lims = [-500,500;-500,500;-5,-2];
xyr = MTADxyz('data',[xyz(:,7,1),xyz(:,7,2),rhmp],'sampleRate',xyz.sampleRate);

pfs = MTAApfs(Trial,1:16,'walk',true,'bhvfettest39',[30,30,.2],[1.2,1.2,2.5],'type','xyr','xyz',xyr,'bound_lims',bound_lims);
pfs.parameters.type = 'xyz';

