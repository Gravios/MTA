


%% COMPUTE pfd HPITCHxBPITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dspch = pch.copy(); 
dspch.resample(5);
dsxyz = xyz.copy();
dsxyz.resample(5);

hfig = figure(666002);
hfig.Units = 'centimeters';
%hfig.Position = [0.5,0.5,35,25];
hfig.Position = [0.5,0.5,16,12];
hfig.PaperPositionMode = 'auto';
ny = 12;
hax = gobjects([1,4]);
for u = 1:numel(units{tind}), 
    clf();    
    maxPfsRate = max([pft.maxRate(units{tind}(u)),pfd{tind,pfindex}.maxRate(units{tind}(u),'mazeMaskFlag',false)]);
    
    % PLOT placefield rate map
    hax(1) = subplot(221);  hold('on');  plot(pft,units{tind}(u),'mean',true,maxPfsRate,false,0.99);
    plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta Place Field, unit:',num2str(units{tind}(u))]);
    
    % PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(2) = subplot(222);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
    colorbar();    
    plot(dspch(drzState{u},2),...
         dspch(drzState{u},1),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('body pitch (rad)');  ylim([-pi/2,pi/2]);
    title('RateMap');

    % PLOT placefield rate map
    hax(3) = subplot(223);  
    hold('on');  
    plot(pft,units{tind}(u),'snr',true,[],false,0.99);
    plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);
    
    % PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(4) = subplot(224);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'snr',true,5,false,0.85,false);
    plot(dspch(drzState{u},2),...
         dspch(drzState{u},1),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('body pitch (rad)');  ylim([-pi/2,pi/2]);
    title('SNR Map');
    
    af(@(h) set(h,'Units','centimeters'),            hax);    
    af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
    af(@(h) set(h.Title,'Units','pixels'),           hax);
    af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);
    % SAVE figure
    drawnow();
    figName = ['rateMap_BHPITCHxBPITCH','_',Trial.filebase,'_unit-',num2str(units{tind}(u))];
    print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
end%for u





%% VISUALIZE erpPCA_HPITCHxBPITCH eigenvectors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reshape_eigen_vector = @(V,pfd) fliplr(rot90(reshape(V(:,1),pfd{1}.adata.binSizes')',-1));

hfig = figure(666003);clf();
hfig.Units = 'centimeters';
%hfig.Position = [0.5,0.5,35,25];
hfig.Position = [0.5,0.5,20,6];
hfig.PaperPositionMode = 'auto';

%dS = diag(S);
nV = 5;
hax = gobjects([1,5]);
for i = 1:nV,
    hax(i) = subplot(1,nV,i);
    fpc = nan([zdims(1),1]);
    fpc(validDimsInds) = V(:,i);
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc,pfd))},[-0.5,3.5],'linear',false,[0,0,0],1,1);       % PRINT eigenvectors
    colorbar();                                           
    %title(sprintf('PC%i Var:%3.2f',i,dS(i)));             % APPEND rank and variance as title
    axis('xy');
    axis('tight');
    hold('on');
    for s = 1:numel(states),                              % OVERLAY state Contours
        copyobj(H{s},gca);
    end
    caxis([-0.5,3.5])
end
af(@(h) set(h,'Units','centimeters'),            hax);    
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
ForAllSubplots('daspect([1,1,1])');

figName = ['erpPCA_HPITCHxBPITCH','_',sessionListName];
print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));




%% VISUALIZE pfd HPITCHxVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dspch = pch.copy();
dspch.resample(5);
dsxyz = xyz.copy();
dsxyz.resample(5);        
dsvxy = dsxyz.vel(['spine_lower'],[1,2]);
dsvxy.data(dsvxy.data<1e-3) = 1e-3;
dsvxy.data = log10(dsvxy.data);
dsfet = dspch.copy('empty');
dsfet.label = 'fet_VP';       
dsfet.data = [dspch(:,2),dsvxy(:)];


hfig = figure(666003);
hfig.Units = 'centimeters';
%hfig.Position = [0.5,0.5,35,25];
hfig.Position = [0.5,0.5,16,12];
hfig.PaperPositionMode = 'auto';
ny = 12;
hax = gobjects([1,4]);
for u = 1:numel(units{tind}), 
    clf();    
    maxPfsRate = max([pft.maxRate(units{tind}(u)),pfd{tind,pfindex}.maxRate(units{tind}(u),'mazeMaskFlag',false)]);
    
    % PLOT placefield rate map
    hax(1) = subplot(221);  hold('on');  plot(pft,units{tind}(u),'mean',true,maxPfsRate,false,0.99);
    plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta Place Field, unit:',num2str(units{tind}(u))]);
    
    % PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(2) = subplot(222);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
    colorbar();    
    plot(dsfet(drzState{u},1),...
         dsfet(drzState{u},2),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
    ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
    title('RateMap');

    % PLOT placefield rate map
    hax(3) = subplot(223);  
    hold('on');  
    plot(pft,units{tind}(u),'snr',true,[],false,0.99);
    plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);
    
    % PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(4) = subplot(224);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'snr',true,5,false,0.85,false);
    plot(dsfet(drzState{u},1),...
         dsfet(drzState{u},2),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
    ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
    title('SNR Map');
    
    af(@(h) set(h,'Units','centimeters'),            hax);    
    af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
    af(@(h) set(h.Title,'Units','pixels'),           hax);
    af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);
    % SAVE figure
    drawnow();
    figName = ['rateMap_BHPITCHxVEL','_',Trial.filebase,'_unit-',num2str(units{tind}(u))];
    print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
end%for u







%% erpPCA_HPITCHxVEL
% VISUALIZE vector X sample space
figure,imagesc(rmaps(nind,sind)');

% VISUALIZE individual sample reshaped vector space
figure,
for u = 1:size(rmaps,2), clf();
    imagesc(reshape(rmaps(:,u),cellfun(@numel,bins))');
    axis('xy');
    waitforbuttonpress();
end
