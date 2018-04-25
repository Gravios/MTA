% placefield_summary
% 
% description: plot a variety of measures for units with place fields 
%
% Mod:20171211: reorganize the signature of multiple functions to include units
%               as the second input option.
%

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';
generateFigureParts = false;
marker = 'nose';

if generateFigureParts,
    FigDir = create_directory('/storage/gravio/figures/analysis/parts/placefields');
else,
    FigDir = create_directory('/storage/gravio/figures/analysis/placefields'); 
end        
        
% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);


% $$$ Sessions  = cf(@(t)  MTASession.validate(t.filebase),   Trials);
% $$$ cf(@(s)  create(s.spk,s),   Sessions);
% $$$ cf(@(s)  save(s),           Sessions);
% $$$ clear('Sessions');
% $$$ Trials  = af(@(t)  MTATrial.validate(t),   sessionList);


states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
             'theta-groom-sit'};

numStates = numel(states);



interpPar = struct('bins',{{linspace(-500,500,100),linspace(-500,500,100)}},             ...
                   'nanMaskThreshold', 0,                                                ...
                   'methodNanMap',     'linear',                                         ...
                   'methodRateMap',    'linear');


pfdVersion = '7';

% PROCESS session data
% $$$ cf(@(t)  compute_neuron_quality(t,[],[],[],true), Trials);
% $$$ units  =  cf(@(t)  select_placefields(t,30,true),  Trials);
% $$$ cf(@(t)  pfs_2d_theta(t,'overwrite',true),  Trials);
% $$$ cf(@(t)  pfs_2d_states(t,'states',{'hpause&theta'},'overwrite',true),  Trials);
% $$$ units  =  cf(@(t)  select_placefields(t,30,false),  Trials);
% $$$ ind = 21:23;
% $$$ cf(@(t,u)  pfs_2d_theta(t,u,false,true),  Trials(ind),units(ind));
% $$$ cf(@(t,u)  pfs_2d_states(t,u,[],[],false,'',true),  Trials(ind),units(ind));
% $$$ cf(@(t)    req20180123_ver5(t,[],'7',true), Trials(ind));
% $$$ cf(@(t)  compute_pfstats_bs(t,'overwrite',true),  Trials);
% $$$ cf(@(t)  MjgER2016_drzfields(t,[],true), Trials);
% $$$ cf(@(t)  pfd_2d_pp(t), Trials);
% $$$ pft = {};
% $$$ for t = 1:numel(Trials),
% $$$      pft{t} = pfs_2d_theta(Trials{t},[],false,false,1);
% $$$ end



% FOR each Trial -------------------------------------------------------------

hfig = figure(666001);
hfig.Position = [1, 1, 1269, 681];
hfig.PaperPositionMode = 'auto';
ny = 12;

for t = 1:numel(Trials),
    clf();

% LOAD Trial 
% CREATE Trial directory in specified location
% PRINT intent of trial processing
    Trial = Trials{t};    
    Trial.load('nq');
    create_directory(fullfile(FigDir,Trial.filebase));
    disp(['Processing Trial: ' Trial.filebase]);

% COPY stc
% LOAD spike waveforms
% LOAD marker positions
% COMPUTE head pitch
% COREGISTER head pitch to reference session
    stc = Trial.stc.copy();
    spkw = Trial.spk.copy();
    spkw.load_spk(Trial);
    xyz = Trial.load('xyz');
    pch = fet_HB_pitch(Trial);
    map_to_reference_session(pch,Trial,pitchReferenceTrial);    

% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
    lfp = Trial.load('lfp',sessionList(t).thetaRef);
    lfp.resample(xyz);
    phz = lfp.phase([6,12]);

% SELECT high quality units with place fields
% COMPUTE placefields' statistics
% LOAD placefields for theta state
% LOAD spikes 
    units = select_placefields(Trial);
    if isempty(units), continue;end;
    pft = pfs_2d_theta(Trial,units,false,false);    
    %pfstats = compute_pfstats_bs(Trial,units);
    pfstats = [];
    spk = Trial.spk.copy();
    spk.create(Trial,xyz.sampleRate,'theta',units,'deburst');    
    
% LOAD placefields and subsampled estimate
    pfs = pfs_2d_states(Trial,units);
    pfs{end+1} = pft;
    
% COMPUTE behavioral state occupancy maps
    stsocc = cell(size(states));
    for sts = 1:numStates,    
        [stsocc(sts),stsoccBins] = xyocc(Trial,stc{states{sts}});
    end
    maxocc = max(cellfun(@(r) max(r(:)),stsocc));

% COMPUTE drz and ddz
    drz = compute_drz(Trial,units,pft,pfstats,[],marker);
    ddz = compute_ddz(Trial,units,pft,pfstats,[],marker);
    
% COMPUTE unit auto correlogram
    [accg,tbins] = autoccg(Trial);
    
% COMPUTE place fields and subsampled estimate
% $$$     for sts = 1:numStates,
% $$$         [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5,units,4,pft);
% $$$         sper{sts}{1}(sper{sts}{1}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
% $$$         sper{sts}{2}(sper{sts}{2}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
% $$$     end
    
% LOAD DRZ fields
    %dfs = MjgER2016_drzfields(Trial,units,false);%,pfstats);
    %dfst = {'height','rhm','Bpitch','HBpitch'};
    dfs = req20180123_ver5(Trial,[],pfdVersion);
    dfst = {'HPITCHxBPITCH','HPITCHxBSPEED','BPITCHxBSPEED','BPITCHxHSPEED','HPITCHxRHM'};
    
% COMPUTE phase precession
    drzp = drz;    drzp(drzp<0)=nan;
    ddzp = ddz;    ddzp(ddzp<0)=nan;
    drzn = drz;    drzn(drzn>0)=nan;
    ddzn = ddz;    ddzn(ddzn>0)=nan;    
    P  = {};   phzStats  = {};   Rmax  = {};
    PP = {};   phzStatsP = {};   RmaxP = {};
    PN = {};   phzStatsN = {};   RmaxN = {};
    
    
    for s = 1:numStates,
        spkpp = Trial.spk.copy();
        spkpp.create(Trial,xyz.sampleRate,states{s},units,'deburst');
        [P{s},  phzStats{s}, Rmax{s}, drzHCnt{s}] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spkpp,units);
        [PP{s},phzStatsP{s},RmaxP{s},drzHCntP{s}] = MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spkpp,units);
        [PN{s},phzStatsN{s},RmaxN{s},drzHCntN{s}] = MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spkpp,units);
    end
    
    
% $$$     mCom = pfstats.peakPatchCOM(:,:,ismember(pfstats.cluMap,units),:);
% $$$     mCom(mCom==0) = nan;
% $$$     mCom = sq(mean(pfstats.peakPatchCOM,2,'omitnan'));
    
    
    for u = 1:numel(units),  tic        
        %figure(hfig);
        clf();
        sp = gobjects([1,0]);
        unit = units(u);
        
% $$$         mpfsRate = max(cell2mat(cf(@(p,u) max(p.maxRate(u)),...
% $$$                                    pfkbs,repmat({unit},[1,numStates]))));
        maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                      [pfs,dfs],repmat({units(u)},[1,numel(pfs)+numel(dfs)]))));

        pfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean',0.5)),pfs,...
                                  repmat({unit},[1,numStates])));
        dfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean',0.5)),dfs,...
                                  repmat({unit},[1,numel(dfs)])));
        
        %mpfsRate = max([pfsMaxRates]);        
        %mpfsRate = max([dfsMaxRates]);        
        %mpfdRate = max([dfsMaxRates]);
        mpfsRate = maxPfsRate;        
        mpfdRate = maxPfsRate;
        
        pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({unit},[1,numStates])));
        
        
        
% $$$         mccgRate = max(cell2mat(cf(@(c,u) max(max(c.ccg(:,u==c.cluMap(:,1),:,1))),...
% $$$                                    bhvccg,repmat({unit},[1,numel(bhvccg)]))));
% $$$         if mccgRate<=0,mccgRate=1;end % default to 1 if 0
        
        uResw = spkw.res(spkw.clu==unit);
        uSpkw = spkw.spk(spkw.clu==unit,:,:);
        [~,sInd] = SelectPeriods(uResw,[stc{states{end},1}],'d',1,0);
        if numel(sInd)<=1, mspkt = zeros([size(uSpkw,2),size(uSpkw,3)]);
        else,              mspkt = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
        end

        FigInfo = uicontrol('Parent',hfig,...
                            'Style','text',...
                            'String',{['Unit: ',num2str(unit)],...
                                       Trial.filebase,...
                                      ['stcMode: ',Trial.stc.mode],...
                                      ['eDist:   ',num2str(Trial.nq.eDist(unit))],...
                                      ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],...
                                      ['SNR:     ',num2str(Trial.nq.SNR(unit))],...
                                      ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],...
                                      ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]...
                                     },...
                            'Units','centimeters',...
                            'Position',[2,18.5,6,3.5]);

        
        % row1col1 - auto ccg
        sp(end+1) = subplot2(ny,numStates+2,3:4,1);
        bar(tbins,accg(:,unit));axis tight;


        % row1col1 - pitch hl loc distrib
        sp(end+1) = subplot2(ny,numStates+2,5:6,1);hold('on');
        ind = [stc{'lloc'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,3),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'c';
        hax.EdgeColor = 'c';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;
        ind = [stc{'hloc'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,3),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'r';
        hax.EdgeColor = 'r';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;

        % row1col1 - pitch hl pause distrib
        sp(end+1) = subplot2(ny,numStates+2,7:8,1);hold('on');
        ind = [stc{'lpause'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,3),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'c';
        hax.EdgeColor = 'c';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;
        ind = [stc{'hpause'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,3),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'r';
        hax.EdgeColor = 'r';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;
         
        
        
        % row9-12col1 - ave spike waveform
        sp(end+1) = subplot2(ny,numStates+2,[9:12],1);
        hold('on');            
        [~,sInd] = SelectPeriods(uResw,[stc{states{end},1}],'d',1,0);
        if numel(sInd)>5,
            mspk = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
            plot(mspk,'b');
            sspk = sq(std(uSpkw(sInd,:,:)))';
            plot(mspk-sspk,'r');
            plot(mspk+sspk,'r');
        end
        xlim([0,52]);
        ylim([0,9000]);

        
        
        for s = 1:numStates,
% GET current state
% SELECT spikes within current state 
% SKIP if spike count is less than 50
            state = [stc{states{s},spk.sampleRate}];
            res = spk(unit);
            res = res(WithinRanges(res,state.data));
            res(abs(ddz(res,u))>250) = [];
            if numel(res) >10,
                res(res>xyz.size(1))=[];
                drzspk = drz(res,u);
                ddzspk = ddz(res,u);
                phzspk = phz(res,spk.map(spk.map(:,1)==units(u),2));
                gind = ~isnan(drzspk)&~isnan(phzspk);                
            else
                res = [];
                drzspk=[];
                ddzspk=[];
                phzspk=[];
                gind=[];
            end
            
% PLACEFIELDS MTAAknnpfs
            sp(end+1) = subplot2(ny,numStates+2,[1,2],s+1);
            plot(pfs{s},unit,1,'none',mpfsRate,true,[],false,interpPar,@jet);
            set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
            title(sprintf('Max Rate: %3.2f',pfsMaxRates(s)));

            
% PLACEFIELDS MTAApfs
            sp(end+1) = subplot2(ny,numStates+2,[3,4],s+1);
            plot(pfs{s},unit,'mean','none',mpfsRate,true,0.5,false,interpPar,@jet);
            set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
            title(sprintf('Max Rate: %3.2f',pfsMaxRatesMean(s)));
            
% TRANSITION Positions
            sp(end+1) = subplot2(ny,numStates+2,[5,6],s+1);
            hold('on');
            imagesc(stsoccBins{1},stsoccBins{2},stsocc{s}');
            axis('xy');
            set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});                
            caxis([0,maxocc])
            title(statesCcg{s});                            

% TRAJECTORIES Positions
% $$$             if unit==units(1),
% $$$                 subplot2(ny,numStates+2,[3,4],s+1);
% $$$                 hold('on');
% $$$                 rper = stc{states{s},xyz.sampleRate};
% $$$                 for r = 1:size(rper,1),
% $$$                     scatter(xyz(rper(r,:),7,1),xyz(rper(r,:),7,2),5,jet(diff(rper(r,:))+1),'filled')
% $$$                 end
% $$$                 xlim([-500,500]);ylim([-500,500]);
% $$$                 set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
% $$$             end


% PHASE PRECESSION 
% Plot phase drz relationship
            sp(end+1) = subplot2(ny,numStates+2,[7,8],s+1); hold('on');
            if sum(gind)>10,
                hist2([[drzspk(gind);drzspk(gind)],...
                      [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360]],...
                      linspace(-1,1,25),...
                      linspace(-180,540,30));
                xlim([-1,1]);
                ylim([-180,540]);
                
                plot([-1,1],circ_rad2ang(2*pi*P{s}(u,1,1)*[-1,1]+P{s}(u,1,2)),'-m','LineWidth',1)
                plot([-1,1],circ_rad2ang(2*pi*P{s}(u,1,1)*[-1,1]+P{s}(u,1,2))+360,'-m','LineWidth',1)

                plot([0,1],circ_rad2ang(2*pi*PP{s}(u,1,1)*[0,1]+PP{s}(u,1,2)),'-r','LineWidth',1)
                plot([0,1],circ_rad2ang(2*pi*PP{s}(u,1,1)*[0,1]+PP{s}(u,1,2))+360,'-r','LineWidth',1)

                plot([-1,0],circ_rad2ang(2*pi*PN{s}(u,1,1)*[-1,0]+PN{s}(u,1,2)),'-g','LineWidth',1)
                plot([-1,0],circ_rad2ang(2*pi*PN{s}(u,1,1)*[-1,0]+PN{s}(u,1,2))+360,'-g','LineWidth',1)
                
            end

            sp(end+1) = subplot2(ny,numStates+2,[9,10],s+1); hold('on');
            if sum(gind)>10,
                plot([drzspk(gind);drzspk(gind)],...
                     [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
                     '.','MarkerSize',1);
                xlim([-1,1]);
                ylim([-180,540]);
                
                plot([-1,1],circ_rad2ang(2*pi*P{s}(u,1,1)*[-1,1]+P{s}(u,1,2)),'-m','LineWidth',1)
                plot([-1,1],circ_rad2ang(2*pi*P{s}(u,1,1)*[-1,1]+P{s}(u,1,2))+360,'-m','LineWidth',1)

                plot([0,1],circ_rad2ang(2*pi*PP{s}(u,1,1)*[0,1]+PP{s}(u,1,2)),'-r','LineWidth',1)
                plot([0,1],circ_rad2ang(2*pi*PP{s}(u,1,1)*[0,1]+PP{s}(u,1,2))+360,'-r','LineWidth',1)

                plot([-1,0],circ_rad2ang(2*pi*PN{s}(u,1,1)*[-1,0]+PN{s}(u,1,2)),'-g','LineWidth',1)
                plot([-1,0],circ_rad2ang(2*pi*PN{s}(u,1,1)*[-1,0]+PN{s}(u,1,2))+360,'-g','LineWidth',1)
                xlim([-1,1]),
                ylim([-180,540])
            end

            

% TRANITION triggered histogram onset
% TRANITION triggered histogram offset
% $$$             sp(end+1) = subplot2(ny,numStates+2,[11,12],s+1); hold('on');            
% $$$             plot(bhvccg{sts}.tbin,RectFilter(sq(bhvccg{s}.ccg(:,unit==bhvccg{s}.cluMap(:,1),:,1))))
% $$$              %plot(bhvccg{s}.tbin,RectFilter(sq(bhvccg{s}.ccg(:,unit==bhvccg{s}.cluMap(:,1),:))));
% $$$             ylim([0,mccgRate]);

% WAVEFORM of unit
% $$$             sp(end+1) = subplot2(ny,numStates+2,[9:12],s+1);
% $$$             hold('on');            
% $$$             [~,sInd] = SelectPeriods(uResw,[stc{states{s},1}],'d',1,0);
% $$$             if numel(sInd)>5,
% $$$                 mspk = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
% $$$                 plot(bsxfun(@plus,mspk-mspkt,fliplr(linspace(250,2000,size( uSpkw,2)))),'b');
% $$$             end
% $$$             xlim([0,52]);
% $$$             ylim([0,2250]);
% $$$             if s~=1,set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});end
            
        end

        for s = 1:numel(dfs),
            sp(end+1) = subplot2(ny,numStates+2,[s*2:s*2+1]-1,numStates+2);
            dfs{s}.plot(unit,'mean',true,mpfdRate,false,0.85,false,[],@jet);
            hax = colorbar();
            hax.Position(1) = hax.Position(1) + 0.05;
            title(dfst{s});
        end
        

        % FORMAT for figure parts
        if generateFigureParts,
            af(@(h) set(h,'Units','centimeters'), sp);
            af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), sp);
        end

        pause(0.01);
        
        FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit)];
        print(hfig,'-depsc2',fullfile(FigDir,Trial.filebase,[FigName,'.eps']));        
        print(hfig,'-dpng',  fullfile(FigDir,Trial.filebase,[FigName,'.png']));
toc
        delete(sp);        

    end%for units
end%for Trials


