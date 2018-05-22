
% vars:
%    Trials
%    units
%    figDir
%    


analDir = ['parts-v',version];
create_directory(fullfile(figDir,analDir));

if ~exist('pft','var'), 
    pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);
end

interpPar = struct('bins',{{linspace(-500,500,100),linspace(-500,500,100)}},             ...
                   'nanMaskThreshold', 0,                                                ...
                   'methodNanMap',     'linear',                                         ...
                   'methodRateMap',    'linear');



ny = 2;
nx = 7;


for tind = 1:numTrials,
    Trial = Trials{tind};
    
    [accg,tbins] = autoccg(Trial);
    
    hfig = figure(666004);
    hfig.Units = 'centimeters';
    hfig.Position = [0.5,0.5,34,10];
    hfig.PaperPositionMode = 'auto';


    
    for u = 1:numel(units{tind}), 
        hax = gobjects([1,0]);        
        clf();    
        maxPfsRate = max([maxRate(pft{tind},units{tind}(u),false,'prctile99'),...
                          maxRate(pfd{tind,1},units{tind}(u),false,'prctile99'),...                          
                          maxRate(pfd{tind,2},units{tind}(u),false,'prctile99'),...
                          maxRate(pfd{tind,3},units{tind}(u),false,'prctile99'),...                          
                          maxRate(pfd{tind,4},units{tind}(u),false,'prctile99')]);

% PLOT Unit Auto Correlogram
        hax(end+1) = subplot(ny,nx,1);
        bar(tbins,accg(:,units{tind}(u)));axis tight;
        
        
% PLOT placefield RATE map
% PLOT placefield SNR map        
        hax(2) = subplot(ny,nx,2);  
        hold('on');  
        plot(pft{tind},units{tind}(u),'mean',true,maxPfsRate,true,0.99,false,interpPar);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

        hax(end+1) = subplot(ny,nx,9);  
        hold('on');  
        plot(pft{tind},units{tind}(u),'snr',true,[],false,0.99);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

        
% PLOT Rate map HPITCH x BPITCH | DRZ [-0.5,0.5]
% PLOT SNR  map HPITCH x BPITCH | DRZ [-0.5,0.5]        
        hax(end+1) = subplot(ny,nx,3);
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
        ylabel('body pitch (rad)');  ylim([-pi/2,2]);
        title('RateMap');

        hax(end+1) = subplot(ny,nx,10);  
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
        ylabel('body pitch (rad)');  ylim([-pi/2,2]);
        title('SNR Map');

% PLOT Rate map HPITCH x BSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map HPITCH x BSPEED | DRZ [-0.5,0.5]        
        hax(end+1) = subplot(ny,nx,4);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);                
        title('RateMap');

        hax(end+1) = subplot(ny,nx,11);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);
        title('SNR Map');

% PLOT Rate map BPITCH x BSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map BPITCH x BSPEED | DRZ [-0.5,0.5]        
        hax(end+1) = subplot(ny,nx,5);  
        hold('on');  
        plot(pfd{tind,3},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);                
        title('RateMap');

        hax(end+1) = subplot(ny,nx,12);  
        hold('on');  
        plot(pfd{tind,3},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);
        title('SNR Map');

% PLOT Rate map BPITCH x HSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map BPITCH x HSPEED | DRZ [-0.5,0.5]        
        hax(end+1) = subplot(ny,nx,6);  
        hold('on');  
        plot(pfd{tind,4},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('head speed log(cm/s)');  ylim([-2,2]);                
        title('RateMap');

        hax(end+1) = subplot(ny,nx,13);  
        hold('on');  
        plot(pfd{tind,4},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('head speed log(cm/s)');  ylim([-2,2]);
        title('SNR Map');

% PLOT Rate map BPITCH x HSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map BPITCH x HSPEED | DRZ [-0.5,0.5]        
        hax(end+1) = subplot(ny,nx,7);  
        hold('on');  
        plot(pfd{tind,5},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head pitch (rad)');  xlim([-2,pi/2]);
        ylabel('RHM Pow 6-12Hz (A.U.)');  ylim([-9,-3]);                
        title('RateMap');

        hax(end+1) = subplot(ny,nx,14);  
        hold('on');  
        plot(pfd{tind,5},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head pitch (rad)');  xlim([-2,pi/2]);        
        ylabel('RHM Pow 6-12Hz (A.U.)');  ylim([-9,-3]);                        
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



