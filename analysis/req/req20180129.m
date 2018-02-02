% req20180129 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: M: spatial self-representation
%               D: continuation of req20171215.m
%  Bugs: NA



Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial.load('stc','msnn_ppsvd_raux');
units = select_placefields(Trial);
xyz = preproc_xyz(Trial,'trb');
shifts = [-1.5:0.125:1.5];
marker = 'hcom';
states = {'theta-groom-sit','loc&theta','pause&theta'};

%ufr = Trial.load('ufr',xyz,'theta-groom-sit',units,1.25,true);
%pft = pfs_2d_theta(Trial);
%drz = compute_drz(Trial,units,pft);
% $$$ figure();
% $$$ plot(drz(:,units(u)),ufr(:,units(u)),'.');


pftt = cell([1,numel(shifts)]);
pftp = cell([1,numel(shifts)]);
pftl = cell([1,numel(shifts)]);

for s = 1:numel(shifts),  pftt{s} = pfs_2d_state_traj(Trial,shifts(s),marker,'theta-groom-sit','overwrite',true);end 
for s = 1:numel(shifts),  pftp{s} = pfs_2d_state_traj(Trial,shifts(s),marker,'pause&theta','overwrite',true);end 
for s = 1:numel(shifts),  pftl{s} = pfs_2d_state_traj(Trial,shifts(s),marker,'loc&theta','overwrite',true);end 

mrtt = cell2mat(cf(@(p) p.maxRate(), pftt));
mrtp = cell2mat(cf(@(p) p.maxRate(), pftp));
mrtl = cell2mat(cf(@(p) p.maxRate(), pftl));

srtt = cell2mat(cf(@(p) p.maxRate('mode','std'), pftt));
srtp = cell2mat(cf(@(p) p.maxRate('mode','std'), pftp));
srtl = cell2mat(cf(@(p) p.maxRate('mode','std'), pftl));


FigDir = create_directory(fullfile('/storage/gravio/figures/placefields/pfs_traj_shift/',Trial.filebase));
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
ny = 6;
hax = gobjects([ny,numel(pftl)]);
for u = 1:numel(units)
    clf();    
    unit = units(u);
    [mrt,mrit] = max(mrtt(u,:));
    [mrl,mril] = max(mrtl(u,:));
    [mrp,mrip] = max(mrtp(u,:));
    
    [srt,srit] = max(srtt(u,:));
    [srl,sril] = max(srtl(u,:));
    [srp,srip] = max(srtp(u,:));

    
    mrm = max([mrl,mrp,mrt]);
    srm = max([srl,srp,srt]);    
    
% THETA 
    for p = 1:numel(shifts),
        hax(1,p) = subplot2(ny,numel(pftt),1,p);
        plot(pftt{p},unit,'mean',[],mrm);
        title(num2str(shifts(p)));
        if p == 1,    ylabel({['s:loc&theta'],['u:',num2str(unit)],['rate:',sprintf('%2.2f',mrt)]});  end
        clear_axes_ticks(hax(1,p));
        xlabel(sprintf('%2.2f',mrtt(u,p)));
    end
    for p = 1:numel(shifts),
        hax(2,p) = subplot2(ny,numel(pftt),2,p);
        plot(pftt{p},unit,'std',[],srm);
        title(num2str(shifts(p)));                          
        if p == 1,    ylabel( {['s:loc&theta'],['u:',num2str(unit)],['std:',sprintf('%2.2f',srt)]});  end
        clear_axes_ticks(hax(2,p));
        xlabel(sprintf('%2.2f',srtt(u,p)));
    end
% LOC 
    for p = 1:numel(shifts),
        hax(3,p) = subplot2(ny,numel(pftl),3,p);
        plot(pftl{p},unit,'mean',[],mrm);
        title(num2str(shifts(p)));                          
        if p == 1,    ylabel({['s:loc&theta'],['u:',num2str(unit)],['rate:',sprintf('%2.2f',mrl)]});  end
        clear_axes_ticks(hax(3,p));
        xlabel(sprintf('%2.2f',mrtl(u,p)));
    end
    for p = 1:numel(shifts),
        hax(4,p) = subplot2(ny,numel(pftl),4,p);
        plot(pftl{p},unit,'std',[],srm);
        title(num2str(shifts(p)));                          
        if p == 1,    ylabel({['s:loc&theta'],['u:',num2str(unit)],['std:',sprintf('%2.2f',srl)]});  end
        clear_axes_ticks(hax(4,p));
        xlabel(sprintf('%2.2f',srtl(u,p)));
    end
% PAUSE     
    for p = 1:numel(shifts),
        hax(5,p) = subplot2(ny,numel(pftp),5,p);
        plot(pftp{p},unit,'mean',[],mrm);
        if p == 1,    ylabel({['s:pause&theta'],['u:',num2str(unit)],['rate:',sprintf('%2.2f',srp)]});  end
        clear_axes_ticks(hax(5,p));
        xlabel(sprintf('%2.2f',mrtp(u,p)));
    end
    for p = 1:numel(shifts),
        hax(6,p) = subplot2(ny,numel(pftp),6,p);
        plot(pftp{p},unit,'std',[],srm);
        if p == 1,    ylabel({['s:pause&theta'],['u:',num2str(unit)],['std:',sprintf('%2.2f',srp)]});  end
        clear_axes_ticks(hax(6,p));
        xlabel(sprintf('%2.2f',srtp(u,p)));
    end
    
    hax(1,mrit).XColor = [1,0,0];    
    hax(1,mrit).YColor = [1,0,0];        
    hax(2,mrit).XColor = [1,0,0];    
    hax(2,mrit).YColor = [1,0,0];        
    hax(3,mril).XColor = [1,0,0];    
    hax(3,mril).YColor = [1,0,0];        
    hax(4,mril).XColor = [1,0,0];    
    hax(4,mril).YColor = [1,0,0];        
    hax(5,mrip).XColor = [1,0,0];    
    hax(5,mrip).YColor = [1,0,0];        
    hax(6,mrip).XColor = [1,0,0];    
    hax(6,mrip).YColor = [1,0,0];        

    %hfig.Position = [0.5,0.5,50,20];
    af(@(h) set(h,'Units','centimeters'),  hax);
    af(@(h) set(h,'LineWidth',1),     hax);    
    af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]),  hax);
    FigName = ['pfs_traj_shift_',Trial.filebase,'_',marker,'_unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));
end

sil = cell2mat(cf(@(p) p.data.si', pftl));
sip = cell2mat(cf(@(p) p.data.si', pftp));

u = 84;
figure();
plot([shifts',shifts'],[sil(u,:)',sip(u,:)']);
%plot([shifts',shifts'],[mrtl(u,:)',mrtp(u,:)']);
title(['unit: ' num2str(u)]);
xlabel('time shift (s)');
ylabel('max pfs rate Hz');

pf = {pftt,pftl,pftp};

[~,mrttsind] = sort(mrtt(:,round(numel(shifts)/2)));

msnr = nan([numel(units),numel(shifts),3]);
for s = 1:3,
for u = mrttsind',%1:numel(units)
    unit = units(u);
    mrate = plot(pf{s}{round(numel(shifts)/2)},unit,'mean');        
    snr   = plot( pf{s}{round(numel(shifts)/2)},unit,'snr');    
    ind = mrate(:)./prctile(mrate(:),95)>0.5;
    for p = 1:numel(shifts),
        snr   = plot( pf{s}{p},unit,'snr');    
        msnr(u,p,s) = mean(snr(ind));
    end
end
end



hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
for s = 1:3,
subplot(1,3,s);
imagesc(shifts,1:numel(units),bsxfun(@rdivide,msnr(mrttsind,:,s),max(msnr(mrttsind,:,s),[],2)));
caxis([0.5,1]);
axis('xy');
xlabel('time shift (s)');
if s == 1, ylabel('units'); end
title([states{s},': normalized pfs patch snr']);
end
suptitle({[Trial.filebase,': ',marker],'Placefields with time shifted trajectories'});
hfig.Position = [0.5,0.5,28,13];
FigName = ['pfs_traj_shift_meanSNR_',Trial.filebase,'_',marker];
print(hfig,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(FigDir,[FigName,'.png']));
