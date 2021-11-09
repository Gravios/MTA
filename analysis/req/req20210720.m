% req20210720

t = 20
Trial = Trials{t};
unitset = units{t};
sampleRate = 30;
halfSpkWindow = 0.15;
ufrWindow = 0.15;
rateOffset = 1e-3;

xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);

fet =  fet_HB_pitchB(Trial,sampleRate);


spk = Trial.load('spk',             ... object name
                 sampleRate,        ... final samplerate of spike times,...
                 '',                ... state
                 unitset,           ... collection of units
                 'deburst');          % Spike filtering mode
                 
ufr = Trial.load('ufr',             ... object name
                 xyz,               ... reference object
                 spk,               ... collection of spike times
                 unitset,           ... collection of units
                 ufrWindow,         ... window length (seconds)
                 'boxcar',          ... window shape
                 true);               % overwrite

pfs = compute_xyhb_ratemaps( Trial, unitset);

ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
mask = ds.mask;

ratemaps = collate_ratemaps(pfs,unitset,mask,rateOffset);

index = 9360;

xyzi = discretize(sq(xyz(:,'hcom',[1,2])),linspace(-500,500,21));
feti(:,1) = discretize(fet(:,1),linspace(-2,0.8,29));
feti(:,2) = discretize(fet(:,2),linspace(-0.8,2,29));


binGrid = cell([1,numel(pfs.adata.bins)]); 
[binGrid{:}] = ndgrid(pfs.adata.bins{:});
gbinm = nan([size(ratemaps,1),numel(pfs.adata.bins)]);
for d = 1:numel(pfs.adata.bins)
    gbinm(:,d) = cat(2,binGrid{d}(:));
end

sum(gbinm.*repmat(posterior(:),[1,numel(pfs.adata.bins)]),'omitnan')

smoothingWeights = [800.^2,800.^2, 1.2.^2, 1.2.^2];
    
smoothingWeights = diag(smoothingWeights);    


indices = 8900:10560;

figdir = create_directory('/storage/share/Projects/BehaviorPlaceCode/decoding/posteriorVideo');
imgdir = create_directory(...
    fullfile(figdir,[Trial.filebase,'_sr_',num2str(sampleRate),'_frames',num2str(indices([1,end]),'_%d_to_%d')]));

[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'1080p','landscape',[],100,100,10,10);


%%%<<< ADD pcolor of log10(posterior) slice in position space 
[yind, yOffSet, xind, xOffSet] = deal( 1, -100, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2,                                ...
                              fig.subplot.height*2],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
pax_pos_log = pcolor(pfs.adata.bins{1}-25,pfs.adata.bins{2}-25,log10(posterior(:,:,mind(3),mind(4)))');
set(pax_pos_log,'EdgeColor','none');
ln_pos_com_bay_log = plot(0,0,'*r');
ln_pos_sax_log_log = plot(0,0,'*k');
ln_pos_com_log_log = plot(0,0,'*g');
ln_pos_real_log = plot(0,0,'oc');
caxis([-10,-1]);    
xlim([-500,500]);
ylim([-500,500]);
colorbar(sax(end),'EastOutside');
daspect(sax(end),[1,1,1]);
%%%>>>

%%%<<< ADD pcolor of log10(posterior) slice in behavior space 
[yind, yOffSet, xind, xOffSet] = deal( 1, -100, 3, 0);
% CREATE subplot axes
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2,                                ...
                              fig.subplot.height*2],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
pax_bhv_log = pcolor(pfs.adata.bins{3}-0.05,pfs.adata.bins{4}-0.05,log10(sq(posterior(mind(1),mind(2),:,:)))');
set(pax_bhv_log,'EdgeColor','none');
ln_bhv_com_bay_log = plot(0,0,'*r');
ln_bhv_sax_log_log = plot(0,0,'*k');
ln_bhv_com_log_log = plot(0,0,'*g');
ln_bhv_real_log = plot(0,0,'oc');
caxis([-10,-1]);    
xlim([-1.8,0.8]);
ylim([-0.8,1.8]);
colorbar(sax(end),'EastOutside');
daspect(sax(end),[1,1,1]);
%%%>>>

%%%<<< ADD pcolor of posterior slice in position space 
[yind, yOffSet, xind, xOffSet] = deal( 3, -100, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2,                                ...
                              fig.subplot.height*2],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
pax_pos_bay = pcolor(pfs.adata.bins{1}-25,pfs.adata.bins{2}-25,posterior(:,:,mind(3),mind(4))');
set(pax_pos_bay,'EdgeColor','none');
ln_pos_com_bay_bay = plot(0,0,'*r');
ln_pos_sax_log_bay = plot(0,0,'*k');
ln_pos_com_log_bay = plot(0,0,'*g');
ln_pos_real_bay = plot(0,0,'oc');
caxis([0.0000001,0.01]);    
xlim([-500,500]);
ylim([-500,500]);
colorbar(sax(end),'EastOutside');
daspect(sax(end),[1,1,1]);
%%%>>>

%%%<<< ADD pcolor of posterior slice in behavior space 
[yind, yOffSet, xind, xOffSet] = deal( 3, -100, 3, 0);
% CREATE subplot axes
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2,                                ...
                              fig.subplot.height*2],                              ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
pax_bhv_bay = pcolor(pfs.adata.bins{3}-0.05,pfs.adata.bins{4}-0.05,sq(posterior(mind(1),mind(2),:,:))');
set(pax_bhv_bay,'EdgeColor','none');
ln_bhv_com_bay_bay = plot(0,0,'*r');
ln_bhv_sax_log_bay = plot(0,0,'*k');
ln_bhv_com_log_bay = plot(0,0,'*g');
ln_bhv_real_bay    = plot(0,0,'oc');
caxis([0.0000001,0.01]);    
xlim([-1.8,0.8]);
ylim([-0.8,1.8]);
colorbar(sax(end),'EastOutside');
daspect(sax(end),[1,1,1]);
%%%>>>

%%%<<< ADD animated line for x position
[yind, yOffSet, xind, xOffSet] = deal( 1, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
xlim(indices([1,end]));
ylim([-500,500]);
ln_x_real_bay = animatedline(sax(end));
ln_x_real_bay.LineStyle = '-';
ln_x_real_bay.Color = 'c';
ln_x_com_bay = animatedline(sax(end));
ln_x_com_bay.LineStyle = ':';
ln_x_com_bay.Color = 'r';
ln_x_com_log = animatedline(sax(end));
ln_x_com_log.LineStyle = ':';
ln_x_com_log.Color = 'r';
ln_x_sax_log = animatedline(sax(end));
ln_x_sax_log.LineStyle = ':';
ln_x_sax_log.Color = 'k';
ylabel('mm');
sax(end).XTickLabel = [];
%%%>>>

%%%<<< ADD animated line for y position
[yind, yOffSet, xind, xOffSet] = deal( 2, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
xlim(indices([1,end]));
ylim([-500,500]);
ln_y_real_bay = animatedline(sax(end));
ln_y_real_bay.LineStyle = '-';
ln_y_real_bay.Color = 'c';
ln_y_com_bay = animatedline(sax(end));
ln_y_com_bay.LineStyle = ':';
ln_y_com_bay.Color = 'r';
ln_y_com_log = animatedline(sax(end));
ln_y_com_log.LineStyle = ':';
ln_y_com_log.Color = 'r';
ln_y_sax_log = animatedline(sax(end));
ln_y_sax_log.LineStyle = ':';
ln_y_sax_log.Color = 'k';
ylabel('mm');
sax(end).XTickLabel = [];
%%%>>>

%%%<<< ADD animated lines for head pitch
[yind, yOffSet, xind, xOffSet] = deal( 3, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
xlim(indices([1,end]));
ylim([-1.8,0.8]);
ln_h_real_bay = animatedline(sax(end));
ln_h_real_bay.LineStyle = '-';
ln_h_real_bay.Color = 'c';
ln_h_com_bay = animatedline(sax(end));
ln_h_com_bay.LineStyle = ':';
ln_h_com_bay.Color = 'r';
ln_h_com_log = animatedline(sax(end));
ln_h_com_log.LineStyle = ':';
ln_h_com_log.Color = 'r';
ln_h_sax_log = animatedline(sax(end));
ln_h_sax_log.LineStyle = ':';
ln_h_sax_log.Color = 'k';
ylabel('rad');
sax(end).XTickLabel = [];
%%%>>>

%%%<<< ADD animated lines for body pitch
[yind, yOffSet, xind, xOffSet] = deal( 4, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
xlim(indices([1,end]));
ylim([-0.8,1.8]);
ln_b_real_bay = animatedline(sax(end));
ln_b_real_bay.LineStyle = '-';
ln_b_real_bay.Color = 'c';
ln_b_com_bay = animatedline(sax(end));
ln_b_com_bay.LineStyle = ':';
ln_b_com_bay.Color = 'r';
ln_b_com_log = animatedline(sax(end));
ln_b_com_log.LineStyle = ':';
ln_b_com_log.Color = 'r';
ln_b_sax_log = animatedline(sax(end));
ln_b_sax_log.LineStyle = ':';
ln_b_sax_log.Color = 'k';
ylabel('rad');
sax(end).XTickLabel = [];
%%%>>>


%%%<<< ADD animated lines for unit count
[yind, yOffSet, xind, xOffSet] = deal( 5, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
xlim(indices([1,end]));
ylim([0,25]);
ln_ucnt = animatedline(sax(end));
ln_ucnt.LineStyle = '-';
ln_ucnt.Color = 'k';
ylabel({'unit','count'});
sax(end).XTickLabel = [];
%%%>>>

%%%<<< ADD states timeseries
[yind, yOffSet, xind, xOffSet] = deal( 6, 1, 5, 50);
sax(end+1) = axes('Units',fig.subplot.units,                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*9,                              ...
                              fig.subplot.height],                             ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
plotSTC(Trial.stc,sampleRate,'text',{'theta','lpause','lloc','hpause','hloc','rear'},'kbbggr');
xlim(indices([1,end]));

%%%>>>


frameCount = 1;
for index = indices;

    
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar( Trial, index, sampleRate, ufr, ratemaps, halfSpkWindow);

    if isempty(posterior)
        imwrite(frame2im(getframe(hfig)),fullfile(imgdir,['posterior_',num2str(frameCount,'%08.f'),'.jpeg']), 'jpeg');
        frameCount = frameCount+1;
        continue;
    end
    
    lpost = log10(posterior(:));
    lpost = lpost+8;
    lpost(lpost<=0) = eps;
    lpost = lpost./sum(lpost,'omitnan');
    lpos = sum(gbinm.*repmat(lpost(:),[1,numel(pfs.adata.bins)]),'omitnan');

    cpos = sum(gbinm.*repmat(posterior(:),[1,numel(pfs.adata.bins)]),'omitnan');

    posterior = reshape(posterior,pfs.adata.binSizes');
    mind = LocalMinimaN(-posterior,0,100);
    
    wbinm = bsxfun(@minus,gbinm,gbinm(sub2ind(pfs.adata.binSizes',mind(1),mind(2),mind(3),mind(4)),:));                
    weights = exp(multiprod(-wbinm,multiprod(inv(smoothingWeights),wbinm,[1,2],[2]),[2],[2]));
    weights = weights./sum(weights,'omitnan');
    weights = bsxfun(@times,weights,repmat(lpost(:),[1,numel(pfs.adata.bins)]));
    weights = weights./sum(weights,'omitnan');
    wpos = sum(gbinm.*weights,'omitnan');
    
    hfig.CurrentAxes = sax(1);
    pax_pos_log.CData = log10(posterior(:,:,mind(3),mind(4)))';
    ln_pos_real_log.XData = xyz(index,'hcom',1);
    ln_pos_real_log.YData = xyz(index,'hcom',2);
    ln_pos_com_bay_log.XData = cpos(1);    ln_pos_com_bay_log.YData = cpos(2);
    ln_pos_sax_log_log.XData = wpos(1);    ln_pos_sax_log_log.YData = wpos(2);
    ln_pos_com_log_log.XData = lpos(1);    ln_pos_com_bay_log.YData = lpos(2);

    hfig.CurrentAxes = sax(2);    
    pax_bhv_log.CData = log10(sq(posterior(mind(1),mind(2),:,:)))';
    ln_bhv_real_log.XData = fet(index,1);
    ln_bhv_real_log.YData = fet(index,2);
    ln_bhv_com_bay_log.XData = cpos(3);    ln_bhv_com_bay_log.YData = cpos(4);
    ln_bhv_sax_log_log.XData = wpos(3);    ln_bhv_sax_log_log.YData = wpos(4);
    ln_bhv_com_log_log.XData = lpos(3);    ln_bhv_com_log_log.YData = lpos(4);    
    
    hfig.CurrentAxes = sax(3);
    pax_pos_bay.CData = posterior(:,:,mind(3),mind(4))';
    ln_pos_real_bay.XData = xyz(index,'hcom',1);
    ln_pos_real_bay.YData = xyz(index,'hcom',2);
    ln_pos_com_bay_bay.XData = cpos(1);    ln_pos_com_bay_bay.YData = cpos(2);    
    ln_pos_sax_log_bay.XData = wpos(1);    ln_pos_sax_log_bay.YData = wpos(2);
    ln_pos_com_log_bay.XData = lpos(1);    ln_pos_com_log_bay.YData = lpos(2);

    hfig.CurrentAxes = sax(4);    
    pax_bhv_bay.CData = sq(posterior(mind(1),mind(2),:,:))';
    ln_bhv_real_bay.XData = fet(index,1);
    ln_bhv_real_bay.YData = fet(index,2);
    ln_bhv_com_bay_bay.XData = cpos(3);    ln_bhv_com_bay_bay.YData = cpos(4);
    ln_bhv_sax_log_bay.XData = wpos(3);    ln_bhv_sax_log_bay.YData = wpos(4);
    ln_bhv_com_log_bay.XData = lpos(3);    ln_bhv_com_log_bay.YData = lpos(4);
        
    hfig.CurrentAxes = sax(5);    
    ln_x_real_bay.addpoints(index,xyz(index,'hcom',1));
    ln_x_com_bay.addpoints(index,cpos(1));
    ln_x_com_log.addpoints(index,lpos(1));
    ln_x_sax_log.addpoints(index,wpos(1));

    hfig.CurrentAxes = sax(6);    
    ln_y_real_bay.addpoints(index,xyz(index,'hcom',2));
    ln_y_com_bay.addpoints(index,cpos(2));
    ln_y_com_log.addpoints(index,lpos(2));
    ln_y_sax_log.addpoints(index,wpos(2));

    hfig.CurrentAxes = sax(7);    
    ln_h_real_bay.addpoints(index,fet(index,1));
    ln_h_com_bay.addpoints(index,cpos(3));
    ln_h_com_log.addpoints(index,lpos(3));
    ln_h_sax_log.addpoints(index,wpos(3));

    hfig.CurrentAxes = sax(8);    
    ln_b_real_bay.addpoints(index,fet(index,2));
    ln_b_com_bay.addpoints(index,cpos(4));
    ln_b_com_log.addpoints(index,lpos(4));
    ln_b_sax_log.addpoints(index,wpos(4));

    hfig.CurrentAxes = sax(9);    
    ln_ucnt.addpoints(index,unitCount);

    
    drawnow('update');
    imwrite(frame2im(getframe(hfig)),fullfile(imgdir,['posterior_',num2str(frameCount,'%08.f'),'.jpeg']), 'jpeg');
    frameCount = frameCount+1;
end

%    IMWRITE(A,FILENAME,FMT) writes the image A to the file specified by
cwd = pwd();
cd(imgdir);
system('ffmpeg -r 1/1 -i posterior_%08d.jpeg -c:v libx264 -vf fps=30 -pix_fmt yuv1080p out.mp4')
cd(cwd);


rmap = ratemaps(:,2);
rmap = reshape(rmap,pfs.adata.binSizes');


rmap  = pfs.plot(unitset(2),1,false,[],false,0.25,false);


figure();
set(pcolor(rmap(:,:,mind(3),mind(4)))','EdgeColor','none');


figure();
%set(pcolor(mask(:,:,mind(3),mind(4)))','EdgeColor','none');
set(pcolor(sq(mask(mind(1),mind(2),:,:)))','EdgeColor','none');
ffmpeg -framerate 1 -pattern_type glob -i '*.jpeg'  -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

ffmpeg -f image2 -r 60 -i path/filename%03d.jpg -vcodec libx264 -crf 18  -pix_fmt yuv420p test.mp4

ffmpeg -r 30 -f image2 -s 1920x1080 -i posterior_%08d.png -vcodec libx264 -crf 18  -pix_fmt yuv420p test.mp4

mthd = 'esax';
mthd = 'ecom';
figure,
for sts = 1:6,
    subplot2(6,2,sts,1);
    ind = dca.stcm(:,1)==1 & dca.stcm(:,sts)==sts;
    hist2([dca.ucnt(ind),...
           dca.(mthd)(ind,3)],...
          linspace(0,30,31),...
          linspace(-1.2,1.2,30),'yprob');
    title(states{sts});
    caxis([0,0.1])

    subplot2(6,2,sts,2);
    hist2([dca.ucnt(ind),...
           dca.(mthd)(ind,4)],...
          linspace(0,30,31),...
          linspace(-1.2,1.2,30),'yprob');
    title(states{sts});
    caxis([0,0.1])
end



mthd = 'esax';
figure();
ind = dca.stcm(:,1)==1 & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6);
for uc = 1:25,
subplot2(5,5,);
hist2([sqrt(sum(dca.(mthd)(ind,1:2).^2,2)),sqrt(sum(dca.(mthd)(ind,3:4).^2,2))],30, ...
      30,[],'mud');


mthd = 'com';
figure,
hist2(dca.(mthd)(ind,1:2),linspace(-500,500,30),linspace(-500,500,30),[],'mud')