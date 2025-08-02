% VISUALIZATION script for req20180123.m



hfig = figure(666003);
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,16,12];
hfig.PaperPositionMode = 'auto';
ny = 12;
hax = gobjects([1,4]);
for u = 1:numel(units{tind}), 
    clf();    
    maxPfsRate = max([pft.maxRate(units{tind}(u),false),pfd{tind,pfindex}.maxRate(units{tind}(u),false)]);
    
% PLOT theta placefield rate map            
% PLOT theta placefield SNR map            
    hax(1) = subplot(221);  hold('on');  plot(pft,units{tind}(u),'mean',true,maxPfsRate,false,0.99);
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

    hax(3) = subplot(223);  
    hold('on');  
    plot(pft,units{tind}(u),'snr',true,[],false,0.99);
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

    
% PLOT rate map HPITCH x BSPEED | DRZ[-0.5,0.5]
% PLOT SNR  map HPITCH x BSPEED | DRZ[-0.5,0.5]
    hax(2) = subplot(222);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
    colorbar();    
    xlabel('head-body pitch (rad)');   xlim([-2,pi/2]);
    ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
    title('RateMap');

    hax(4) = subplot(224);  
    hold('on');  
    plot(pfd{tind,pfindex},units{tind}(u),'snr',true,5,false,0.85,false);
    xlabel('head-body pitch (rad)');   xlim([-2,pi/2]);
    ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
    title('SNR Map');
    
% SET figure parameters
    af(@(h) set(h,'Units','centimeters'),            hax);    
    af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax);
    af(@(h) set(h.Title,'Units','pixels'),           hax);
    af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);

% SAVE figure
    drawnow();
    figName = ['rateMap_',tags{pfindex},'_v',version,'_',Trial.filebase,'_unit-',num2str(units{tind}(u))];    
    print(hfig,'-depsc2',fullfile(figDir,analDir,[figName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,analDir,[figName,'.png']));
end%for u

