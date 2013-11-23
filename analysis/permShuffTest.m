




psigw = 1./sum(permute(repmat(permute(pfknnmrw,circshift([1:ndims(pfknnmrw)+1]',1)'),niter,1),circshift([1:ndims(pfknnmrw)+1]',-1)')>pfknnmrs,3);
psigw = 1./sum(permute(repmat(permute(pfknnmrw,circshift([1:ndims(pfknnmrw)+1]',1)'),niter,1),circshift([1:ndims(pfknnmrw)+1]',-1)')>repmat(max(max(pfknnmrs)),50,50),3);
psigw(isinf(psigw))=nan;
figure,
imagescnan({xbins,ybins,psigw},[],[],1,[0,0,0]);,axis xy 
figure,imagesc(psigw)


psigr = 1./sum(permute(repmat(permute(pfknnmrr,circshift([1:ndims(pfknnmrr)+1]',1)'),niter,1),circshift([1:ndims(pfknnmrr)+1]',-1)')>pfknnmrs,3);
psigr = 1./sum(permute(repmat(permute(pfknnmrr,circshift([1:ndims(pfknnmrr)+1]',1)'),niter,1),circshift([1:ndims(pfknnmrr)+1]',-1)')>repmat(max(max(pfknnmrs)),50,50),3);
psigr(isinf(psigr))=nan;
figure,
imagescnan({xbins,ybins,psigr},[],[],1,[0,0,0]);,axis xy 
figure,imagesc(psigr)



psubw  =  pfknnmrw-mean(pfknnmrs,3);
psubr  =  pfknnmrr-mean(pfknnmrs,3);

figure
subplot2(2,3,1,1)
imagescnan({xbins,ybins,psubr},[],[],1,[0,0,0]);,axis xy 
subplot2(2,3,2,1)
imagescnan({xbins,ybins,psubw},[],[],1,[0,0,0]);,axis xy 
subplot2(2,3,1,2)
imagescnan({xbins,ybins,pfknnmrr},[],[],1,[0,0,0]);,axis xy 
subplot2(2,3,2,2)
imagescnan({xbins,ybins,pfknnmrw},[],[],1,[0,0,0]);,axis xy 
subplot2(2,3,1,3)
pfr.plot(57,1)
subplot2(2,3,2,3)
pfw.plot(57,1)