Trial = MTATrial('jg05-20120310');

Trial.load('stc','nn0317_PP');

[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong',true,'overwrite',true);

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,40,'low');

ang = create(MTADang,Trial,xyz);

tx = [0:(xyz.size(1)-1)]/xyz.sampleRate;
tx = MTADxyz('data',tx,'sampleRate',xyz.sampleRate);

th = tx(Trial.stc{'hwalk'});
tl = tx(Trial.stc{'lwalk'});

hper = Trial.stc{'hwalk'};hper.resample(1);
lper = Trial.stc{'lwalk'};lper.resample(1);
rper = Trial.stc{'rear'};rper.resample(1);


% Figure Props
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8);
hfig = figure(2016070601);
clf;
set(hfig,'units','centimeters')
set(hfig,'Position',[0,0,15,8])
set(hfig,'PaperPositionMode','auto');

sp = [];

sp(1) = subplot(2,1,1);

% Plot Head Pitch
plot(tx.data,ang(:,5,7,2),'k'),hold on,
Lines([],-.47,'r')
ylabel('Head Pitch (rad)')

% Plot the First Period of each behavior so the legend only
% displays important elements
t = hper.data(1,:)';
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'b');

t = lper.data(1,:)',
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'g');

t = rper.data(1,:)',
psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
ptx = tx(psts);
ptx = [ptx(1),ptx,ptx(end)];
fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'r');
legend('Head Pitch','Angle Threshold','High Walk','Low Walk','Rearing');

for t = hper.data(:,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(1,psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'b');
end


for t = lper.data(:,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(1,psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'g');
end


for t = rper.data(:,:)',
    psts = find(tx>=t(1),1,'first'):find(tx>=t(2),1,'first');
    ptx = tx(psts);
    ptx = [ptx(1),ptx,ptx(end)];
    fill(ptx,[-.47,ang(psts,5,7,2)',-.47],'r');
end


sp(2) = subplot(2,1,2);
imagesc([0:(rhm.size(1)-1)]/rhm.sampleRate,fs,log10(rhm.data)');
axis xy;
caxis([-5,-2.8]);
linkaxes(sp,'x');
ylabel('RHM Frequencey Hz');
xlabel('Time (s)');
colormap('jet');

print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_4',...
                     'F4_sup_pRHM_pHANG.eps'))



% RHM Distrb
Trial = MTATrial('jg05-20120310');
fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
afig = hgload(fpath);
caxis([.2,1.5])
xlim([-.85,.85])

caxis([.4,1.6])
xlim([1.7,2.2])

caxis([1,2])
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_4',...
                     'F4_sup_pRHM_distrb.eps'))


Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev2_Ed');
bhv_mean_coherence(Trial);


xlim ([.4,2.4])
caxis([.2,1.6])

xlim ([-1.2,1])
caxis([.1,.9])

xlim ([-1.2,1])
caxis([.2,1.6])


%V
xlim ([.1,.55])
caxis([.1,1.1])

xlim ([.1,.55])
caxis([.2,1.8])


%
xlim ([-.6,.6])
caxis([.1,.9])

xlim ([-.6,.6])
caxis([.1,1.3])

%
xlim ([-2.2,1.7])
caxis([0,.9])
caxis([0,2])

print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_4',...
                     'F4_sup_bhv_mean_cohere_Ed01.eps'))





Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev2_Ed');
[rhm] = fet_rhm(Trial,[],'mta');
[ncp] = fet_ncp(Trial,[],'mta',2);

figure,
plot(nunity(ncp.data))
hold on
plot(nunity(rhm.data))

legend({'Nasal Pressure Sensor','Rhythmic Head Motion'})

print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_4',...
                     'F4_sup_rhm_ncp_trace_Ed01.eps'))
