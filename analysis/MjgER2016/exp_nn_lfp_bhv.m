Trial = MTATrial.validate('jg05-20120309.cof.all');
Trial.load('stc','NN0317R');

states = {'loc','rear','pause','groom','sit'};

chans = [65:1:96];

Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);





figure,
chan =4;
sp(1) = subplot(611);
imagesc(th,fh,log10(bsxfun(@rdivide,abs(yhd.data(:,:,chan)),sum(abs(yhd.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-3,-.7])
sp(2) = subplot(612);
imagesc(tl,fl,log10(bsxfun(@rdivide,abs(yld.data(:,:,chan)),sum(abs(yld.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-2,-1.6])

chan =14;
sp(3) = subplot(613);
imagesc(th,fh,log10(bsxfun(@rdivide,abs(yhd.data(:,:,chan)),sum(abs(yhd.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-3,-.7])
sp(4) = subplot(614);
imagesc(tl,fl,log10(bsxfun(@rdivide,abs(yld.data(:,:,chan)),sum(abs(yld.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-2,-1.6])

chan =20;
sp(5) = subplot(615);
imagesc(th,fh,log10(bsxfun(@rdivide,abs(yhd.data(:,:,chan)),sum(abs(yhd.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-3,-.7])
sp(6) = subplot(616);
imagesc(tl,fl,log10(bsxfun(@rdivide,abs(yld.data(:,:,chan)),sum(abs(yld.data(:,:,chan)),2)))')
axis xy,colormap jet,caxis([-2,-1.6])

linkaxes(sp,'x');



figure,
chan =4;
sp(1) = subplot(611);
imagesc(th,fh,log10(abs(yhd.data(:,:,chan)))')
axis xy,colormap jet,caxis([1,4])
sp(2) = subplot(612);
imagesc(tl,fl,log10(abs(yld.data(:,:,chan)))')
axis xy,colormap jet,caxis([-.2,.6])

chan =14;
sp(3) = subplot(613);
imagesc(th,fh,log10(abs(yhd.data(:,:,chan)))')
axis xy,colormap jet,caxis([1.5,4])
sp(4) = subplot(614);
imagesc(tl,fl,log10(abs(yld.data(:,:,chan)))')
axis xy,colormap jet,caxis([.2,.5])

chan =20;
sp(5) = subplot(615);
imagesc(th,fh,log10(abs(yhd.data(:,:,chan)))')
axis xy,colormap jet,caxis([1.5,4])
sp(6) = subplot(616);
imagesc(tl,fl,log10(abs(yld.data(:,:,chan)))')
axis xy,colormap jet,caxis([0.2,.5])

linkaxes(sp,'x');


yldt = yld.copy;
yhdt = yhd.copy;


data = reshape(cat(2,yld(:,2:2:end,1:2:end),yhd(:,2:2:end,1:2:end)),size(yhd,1),[]);


fet = MTADfet.encapsulate(Trial,data,yhd.sampleRate,'PSD depth lfp','lfpPSD','y');


[Stc,d_state,Model_Information,n_state] = bhv_nn(Trial,true,states,Trial.stc,fet,'nNeurons',200);