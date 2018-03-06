

analDir = ['parts-v',version];
create_directory(fullfile(figDir,analDir));

for tind = 1:numTrials,
    hfig = figure(666004);
    hfig.Units = 'centimeters';
    hfig.Position = [0.5,0.5,20,12];
    hfig.PaperPositionMode = 'auto';
    ny = 12;
    hax = gobjects([1,6]);
    for u = 1:numel(units{tind}), 
        clf();    
        maxPfsRate = max([maxRate(pft{tind}  ,units{tind}(u),false,'mean'),...
                          maxRate(pfd{tind,1},units{tind}(u),false,'mean'),...
                          maxRate(pfd{tind,2},units{tind}(u),false,'mean')]);

% PLOT placefield RATE map
% PLOT placefield SNR map        
        hax(1) = subplot(231);  hold('on');  plot(pft{tind},units{tind}(u),'mean',true,maxPfsRate,false,0.99);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

        hax(2) = subplot(234);  
        hold('on');  
        plot(pft{tind},units{tind}(u),'snr',true,[],false,0.99);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

        
% PLOT Rate map HPITCH x BPITCH | DRZ [-0.5,0.5]
% PLOT SNR  map HPITCH x BPITCH | DRZ [-0.5,0.5]        
        hax(3) = subplot(232);  
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body pitch (rad)');  ylim([-2,pi/2]);
        title('RateMap');

        hax(4) = subplot(235);  
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-2,2]);
        ylabel('body pitch (rad)');  ylim([-2,2]);
        title('SNR Map');

% PLOT Rate map HPITCH x BSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map HPITCH x BSPEED | DRZ [-0.5,0.5]        
        hax(5) = subplot(233);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);                
        title('RateMap');

        hax(6) = subplot(236);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);
        title('SNR Map');
        
        
% FORMAT figure
        af(@(h) set(h,'Units','centimeters'),            hax);    
        af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax(1:2));
        af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax(3:end));        
        af(@(h) set(h.Title,'Units','pixels'),           hax);
        af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);

% SAVE figure
        drawnow();
        figName = ['rateMap_BHPITCHxBPITCHxBSPEED_v',version,'_',Trials{tind}.filebase,'_unit-',num2str(units{tind}(u))];
        print(hfig,'-depsc2',fullfile(figDir,analDir,[figName,'.eps']));        
        print(hfig,'-dpng',  fullfile(figDir,analDir,[figName,'.png']));
    end%for u
end%for tind



