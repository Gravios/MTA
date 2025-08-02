% req20190530 forbidden analysis

MjgER2016_load_data();
pfs = cf(@(t,u) pfs_2d_states(t,u), Trials,units);

Trial = Trials{20};


fet = fet_bref(Trial);


figure();
plot(fet(:,1:5));
figure();
plot(fet(:,6:10));
legend({'sl','pr','sm','su','hc'})

figure();
plot(fet(:,17:2:25));
legend({'sl','pr','sm','su','hc'})

dhl = MTADfet.encapsulate(Trial,[0,0,0,0,0;diff(RectFilter(fet(:,16:2:24),3,5))],fet.sampleRate,'dlhm','dlhm','d');
figure,plot(dhl.data);


[ys,fs,ts,pho,fst] = fet_spec(Trial,dhl,'flagCrossSpec',true);

numChan = size(ys,3);
figure(1);
for s = 1:numChan,
    subplot(numChan,1,s);
    imagesc(ts,fs,nunity(ys(:,:,s,s))');
    caxis([-0,0.2]);
    %imagesc(ts,fs,log10(ys(:,:,s,s)'));    
    %caxis([-8-s,-5+s]);
    axis('xy');
end

figure(2);
for s = 2:numChan,
    subplot(numChan+1,1,s);
    imagesc(ts,fs,pho(:,:,1,s)');
    axis('xy');
    colormap('hsv');
    caxis([-pi,pi]);
    %    ylim([0,20]);
end
subplot(numChan+1,1,1);
    imagesc(ts,fs,circ_dist(pho(:,:,1,2),pho(:,:,2,3))');
    axis('xy');
    colormap('hsv');
    caxis([-pi,pi]);
    %        ylim([0,20]);
subplot(numChan+1,1,numChan+1);
plotSTC(Trial.stc,1);a

linkaxes([findobj(figure(1),'Type','Axes');findobj(figure(2),'Type','Axes')],'x');