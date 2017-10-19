% req20170925 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: new rhm feature base on proper filtering
%
%  Bugs: NA


%addpath /storage/weiwei/matlab/NewFolder/
%nrhm = DeConductance(log10(rhm(:,10:100)),[1,91],[]);

%Trial = MTATrial.validate('Ed01-20140707.cof.all');

Trial = MTATrial.validate('Ed05-20140529.ont.all');
[rhmo,fs,ts] = fet_rhm_old(Trial,[],'mtchglong');
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
[ncp,fs,ts] = fet_ncp(Trial,[],'mtchglong');


figure,
subplot(311);
imagesc(ts,fs,log10(bsxfun(@rdivide,rhm.data,nansum(rhm.data')')'));
%imagesc(ts,fs,log10(rhm.data)');
%caxis([-6,-4]);
grid('on');
%imagesc(ts,fs(6:50),nrhm');grid('on');
axis('xy');
colormap('jet');
subplot(312);
imagesc(ts,fs,log10(rhmo.data)');
caxis([-6,-3]);
grid('on');
%imagesc(ts,fs,log10(rhm.data)');caxis([-9,-4]);grid('on');
axis('xy');
colormap('jet');
subplot(313);
imagesc(ts,fs,log10(ncp.data)');caxis([3,6]);grid('on');
axis('xy');
colormap('jet');
linkaxes(findobj(gcf,'Type','Axes'),'xy');

OwnDir  = '/storage/gravio/nextcloud/';
FigDir  = 'Shared/Behavior Paper/Figures/Suplementary';
FigName = 'comparison_of_versions_of_rhm';
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
