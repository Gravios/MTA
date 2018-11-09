
% Theta phase restricted egocentric placefield ratemaps


MjgER2016_load_data();

sampleRate = 250;
pfsState = 'theta-groom-sit-rear';
spkMode = 'deburst';
binPhzs = linspace(-pi,pi,13);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;

for tind = 17:23; % jg05-20120312.cof.all

    Trial = Trials{tind};
    unitSubset = units{tind};


    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    lfp = load(Trial,'lfp',sessionList(tind).thetaRef);
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    

    fxyz = filter(copy(xyz),'ButFilter',3,20,'low');
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    
    pft = pfs_2d_theta(Trial,unitSubset);

    spk = Trial.spk.copy;
    spk.create(Trial,sampleRate,pfsState,unitSubset,spkMode);
    
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
    
    pfTemp = Trial;

    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.units        = unitSubset;
    pargs.tag          = 'egofield';
    pargs.binDims      = [20 20];
    pargs.SmoothingWeights = [1.8 1.8];
    pargs.halfsample   = false;
    pargs.numIter      = 1;   
    pargs.boundaryLimits = [-400,400;-400,400];
    pargs.states       = '';
    pargs.overwrite    = true;
    pargs.autoSaveFlag = false;    
    electrode = 0;
    
    for p = 1:numel(binPhzc)
        pargs.tag          = ['egofield_theta_phase_',num2str(p)];
        for u = 1:numel(unitSubset),
            if u==1 | electrode~=spk.map(spk.map(:,1)==unitSubset(u),2), % update phase state            
                pargs.spk = copy(spk);
                electrode = spk.map(spk.map(:,1)==unitSubset(u),2);
                pargs.states = copy(thetaState);
                pargs.states.label = ['thetaCA3_',num2str(p)];
                pargs.states.data( phz(:,electrode)<binPhzs(p) | phz(:,electrode)>binPhzs(p+1)) = 0;
                cast(pargs.states,'TimePeriods');
                resInd = WithinRanges(pargs.spk.res,pargs.states.data);
                pargs.spk.res = pargs.spk.res(resInd);
                pargs.spk.clu = pargs.spk.clu(resInd);
            end 
            
            [mxr,mxp] = pft.maxRate(unitSubset(u));
            pfsCenterHR = MTADfet.encapsulate(Trial,                                               ...
                                              multiprod(bsxfun(@minus,mxp,sq(xyz(:,'hcom',[1,2]))),...
                                                        hvec,2,[2,3]),                             ...
                                              sampleRate,                                          ...
                                              'placefield_center_referenced_to_head',              ...
                                              'pfsCenterHR',                                       ...
                                              'p'                                                  ...
                                              );
            pargs.xyzp = pfsCenterHR;
            pargs.units  = unitSubset(u);
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
            if u==1,
                pfTemp.purge_savefile();
                pfTemp.save();        
            end    
        end
        pfTemp.save();
        pfTemp = Trial;
    end
end


tind = 18;
Trial = Trials{tind};
unitSubset = units{tind};
tprpf = {};
for p = 1:numel(binPhzc),
    tprpf{p} = MTAApfs(Trial,[],[],[],['egofield_theta_phase_',num2str(p)]);
end
pft = pfs_2d_theta(Trial,unitSubset);

figure,
hax = tight_subplot(1,numel(binPhzc)+1,0.01,0.2,0.1);
for u = unitSubset,
    for p = 1:numel(binPhzc),            
        axes(hax(p));
        cla();
        plot(tprpf{p},u,1,'',[0,15],false);
        title({['unit: ' num2str(u)],['phzBin: ' num2str((binPhzc(p))/pi*180)]});        
    end
    axes(hax(end));
    h = gca();
    cla();
    plot(pft,u,'mean','colorbar',[],true);    
    h.Position(3) = hax(end-1).Position(3);
    waitforbuttonpress();
end


tprpf = {};
for tind = 17:23;
    Trial = Trials{tind};
    unitSubset = units{tind};
    for p = 1:numel(binPhzc),
        tprpf{tind-16,p} = MTAApfs(Trial,[],[],[],['egofield_theta_phase_',num2str(p)]);
    end
end


mpfsr = repmat({[]},[1,numel(binPhzc)]);
for tind = 1:size(tprpf,1),
    for p = 1:numel(binPhzc),
    mpfsr{p} = cat(2,mpfsr{p},tprpf{tind,p}.data.rateMap);
    end
end


cf(@(t) t.load('nq'), Trials(17:23));
edist = cf(@(t,u) t.nq.eDist(u), Trials(17:23),units(17:23));
edist = cat(1,edist{:});


mpfsrm = cf(@(m,p,e) reshape(mean(m(:,max(m)'>1&e>20),2,'omitnan'),p.adata.binSizes'),mpfsr,tprpf(1,:),repmat({edist},[1,numel(mpfsr)]));

figure,

clf();
hax = tight_subplot(1,numel(binPhzc),0.01,0.2,0.1);
for p = 1:numel(binPhzc),            
    axes(hax(p));
    cla();
    imagesc(tprpf{p}.adata.bins{:},mpfsrm{p}');
    title({['phzBin: ' num2str((binPhzc(p))/pi*180)]});        
    caxis([0,7.5]);
    grid('on');
    clear_axes_labels(gca);    
    lax = Lines(0,[],'r');
    if p == 1,
        hax(p).XAxisLocation = 'bottom';        
        xlabel('longitudinal (back/front)');
        ylabel('lateral (left/right)');
    end
    axis('xy');    
end
ForAllSubplots('h = gca();h.Units=''centimeters'';h.Position(4) = h.Position(3);');
ForAllSubplots('h = gca();h.GridColor=[1,1,1];');

fax = axes();
fax.Position = [0,0,1,1];
fax.Visible = 'off';
xlim([0,1]);
ylim([0,1]);
tax = text(fax,0.5,0.95,'Mean firing rate of place field in egocentric coordinates partitioned by theta phase. Subject:jg05','HorizontalAlignment','center');

hax(end).Units = 'normalized';
line(hax(end).Position(1)+hax(end).Position(3).*[0.75,1],...
     [0.1,0.1],'LineWidth',1,'Color',[0,0,0]);
text(hax(end).Position(1)+hax(end).Position(3).*0.75,...
     0.05,'10 cm');
cax = colorbar(hax(end));
cax.Position(1) = cax.Position(1)+0.03;

print(gcf,'-dpng',fullfile('/storage/share/Projects/BehaviorPlaceCode/phase_precession',...
                             'ego_pp_meanFiringRate.png'));
print(gcf,'-depsc2',fullfile('/storage/share/Projects/BehaviorPlaceCode/phase_precession',...
                             'ego_pp_meanFiringRate.eps'));



