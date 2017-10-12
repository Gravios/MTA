% req20170925 ----------------------------------------------------
%
%  Status: active
%  Type: Utility
%  Final_Forms: NA
%  Description: Modify xyz data for wei wei. Create virtual marker
%               which best captures the rhythmic head motion 
%               feature.
%  Bugs: NA


addpath /storage/weiwei/data/VRnew/  

which VRobjLoadDataExp

Trial = MTATrial.validate('Ed01-20140707.cof.all');
[rhmo,fs,ts] = fet_rhm_old(Trial,[],'mtchglong');
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
[ncp,fs,ts] = fet_ncp(Trial,[],'mtchglong');


figure,
subplot(311);
imagesc(ts,fs,log10(rhm.data)');
axis('xy');
colormap('jet');
subplot(312);
imagesc(ts,fs,log10(rhmo.data)');
axis('xy');
colormap('jet');
subplot(313);
imagesc(ts,fs,log10(ncp.data)');
axis('xy');
colormap('jet');
linkaxes(findobj(gcf,'Type','Axes'),'xy');