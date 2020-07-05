% PLACEFIELD
% _summary
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
marker = 'hcom';

if generateFigureParts,
    FigDir = create_directory('/storage/gravio/figures/analysis/parts/placefields');
else,
    FigDir = create_directory('/storage/gravio/figures/analysis/placefields'); 
end        
        
% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)    MTATrial.validate(t),   sessionList);
units   = cf(@(t)    select_placefields(t),  Trials); 
units   = cf(@(t,u)  remove_bad_units(t,u),  Trials,units);
%units = req20180123_remove_bad_units(units);

%unitsInt = cf(@(T)  select_units(T,'int'),  Trials); 

%units = cf(@(T)  select_placefields(T),  Trials); 

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


pfdVersion = '9';

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
pause(0.1);
hfig.Position = [1, 1, 2482, 1274];
ny = 15;
sampleRate = 250;

for t = 1:numel(Trials),
    clf();

% LOAD Trial 
% CREATE Trial directory in specified location
% PRINT intent of trial processing
    Trial = Trials{t};    
    Trial.load('nq');
    unitSubset = units{t};
    %unitSubset = [units{t},unitsInt{t}];    
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
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);    
    pch = fet_HB_pitchB(Trial,sampleRate);

% CREATE lowpass filtered xyz object
% COMPUTE basis vector aligned to the head
    fxyz = filter(copy(xyz),'ButFilter',3,20,'low');    
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);

% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
    lfp = Trial.load('lfp',sessionList(t).thetaRefGeneral);
    phz = lfp.phase([5,13]);    
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);


% SELECT high quality units with place fields
% COMPUTE placefields' statistics
% LOAD placefields for theta state
% LOAD spikes 
    
    if isempty(unitSubset), continue;end;
    pft = pfs_2d_theta(Trial,unitSubset,false,'pfsArgsOverride',struct('numIter',1,'halfsample',false),'overwrite',true);    
    %pfstats = compute_pfstats_bs(Trial,unitSubset);
    pfstats = [];
    spk = Trial.spk.copy();
    spk.create(Trial,xyz.sampleRate,'theta',unitSubset,'');    
    
% LOAD placefields and subsampled estimate
    pfs = pfs_2d_states(Trial,unitSubset,[],states(1:end-1),'pfsArgsOverride',struct('numIter',1,'halfsample',false),'overwrite',true);
    pfs{end+1} = pft;
    
% COMPUTE behavioral state occupancy maps
    stsocc = cell(size(states));
    for sts = 1:numStates,    
        [stsocc(sts),stsoccBins] = xyocc(Trial,stc{states{sts}});
    end
    maxocc = max(cellfun(@(r) max(r(:)),stsocc));

% COMPUTE drz and ddz
    drz = compute_drz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
    ddz = compute_ddz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
    hrz = compute_hrz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
    ghz = compute_ghz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
    gdz = compute_gdz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);    
    
% COMPUTE unit auto correlogram
    [accg,tbins] = autoccg(Trial);
    
% COMPUTE place fields and subsampled estimate
% $$$     for sts = 1:numStates,
% $$$         [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5,unitSubset,4,pft);
% $$$         sper{sts}{1}(sper{sts}{1}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
% $$$         sper{sts}{2}(sper{sts}{2}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
% $$$     end

% LOAD DRZ fields
    %dfs = MjgER2016_drzfields(Trial,unitSubset,false);%,pfstats);
    %dfst = {'height','rhm','Bpitch','HBpitch'};
    %dfs = req20180123_ver5(Trial,[],pfdVersion);
    %dfst = {'HPITCHxBPITCH','HPITCHxBSPEED','BPITCHxBSPEED','BPITCHxHSPEED','HPITCHxRHM'};
    dfs = {compute_bhv_ratemaps(Trial,unitSubset,[],[],pft,[],0.8,250,'overwrite',false)};
    dfst = {'HPITCHxBPITCH'};

% COMPUTE phase precession
    drzp = drz;    drzp(drzp<0)=nan;
    ddzp = ddz;    ddzp(ddzp<0)=nan;    
    drzn = drz;    drzn(drzn>0)=nan;
    ddzn = ddz;    ddzn(ddzn>0)=nan;    
    P  = {};   phzStats  = {};   Rmax  = {};
    PP = {};   phzStatsP = {};   RmaxP = {};
    PN = {};   phzStatsN = {};   RmaxN = {};
    rho= {};

    hrzp = hrz;    hrzp(hrzp<0)=nan;
    hdzp = ddz;    hdzp(hdzp<0)=nan;    
    hrzn = hrz;    hrzn(hrzn>0)=nan;
    hdzn = ddz;    hdzn(hdzn>0)=nan;        
    P_H  = {};   phzStats_H  = {};   Rmax_H  = {};
    PP_H = {};   phzStatsP_H = {};   RmaxP_H = {};
    PN_H = {};   phzStatsN_H = {};   RmaxN_H = {};
    rho_H= {};

    gdzp = gdz;    gdzp(gdzp<0)=nan;
    gdzn = gdz;    gdzn(gdzn>0)=nan;
    P_HDG  = {};   phzStats_HDG  = {};   Rmax_HDG  = {};
    PP_HDG = {};   phzStatsP_HDG = {};   RmaxP_HDG = {};
    PN_HDG = {};   phzStatsN_HDG = {};   RmaxN_HDG = {};
    rho_HDG= {};

    ghzp = ghz;    ghzp(ghzp<0)=nan;
    ghzn = ghz;    ghzn(ghzn>0)=nan;
    P_HHG  = {};   phzStats_HHG  = {};   Rmax_HHG  = {};
    PP_HHG = {};   phzStatsP_HHG = {};   RmaxP_HHG = {};
    PN_HHG = {};   phzStatsN_HHG = {};   RmaxN_HHG = {};
    rho_HHG= {};
    
    for s = 1:numStates,
        spkpp = Trial.spk.copy();
        spkpp.create(Trial,xyz.sampleRate,states{s},unitSubset,'deburst');

        % DRZ phase precession
        [P{s},  phzStats{s},  Rmax{s},  rho{s}] = ...
            MjgER2016_phasePrecession(Trial,drz,ddz,phz,spkpp,unitSubset,[],[],[],['placefield_summary_drz',DataHash(states{s})]);
        [PP{s}, phzStatsP{s},RmaxP{s}, rhoP{s}] = ...
            MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spkpp,unitSubset,[],[],[],['placefield_summary_drzp',DataHash(states{s})]);
        [PN{s}, phzStatsN{s},RmaxN{s}, rhoN{s}] = ...
            MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spkpp,unitSubset,[],[],[],['placefield_summary_drzn',DataHash(states{s})]);

        % HRZ phase precession
        [P_H{s},   phzStats_H{s}, Rmax_H{s}, rho_H{s}] = ...
            MjgER2016_phasePrecession(Trial,hrz,ddz,phz,spkpp,unitSubset,[],[],[],['placefield_summary_hrz',DataHash(states{s})]);
        [PP_H{s}, phzStatsP_H{s},RmaxP_H{s}, rhoP_H{s}] = ...
            MjgER2016_phasePrecession(Trial,hrzp,hdzp,phz,spkpp,unitSubset,[],[],[],['placefield_summary_hrzp',DataHash(states{s})]);
        [PN_H{s}, phzStatsN_H{s},RmaxN_H{s}, rhoN_H{s}] = ...
            MjgER2016_phasePrecession(Trial,hrzn,hdzn,phz,spkpp,unitSubset,[],[],[],['placefield_summary_hrzn',DataHash(states{s})]);
        
        % GHZ phase precession
        [P_HHG{s},   phzStats_HHG{s}, Rmax_HHG{s}, rho_HHG{s}] = ...
            MjgER2016_phasePrecession(Trial,ghz,ddz,phz,spkpp,unitSubset,[],[],[],['placefield_summary_ghz',DataHash(states{s})]);
        [PP_HHG{s}, phzStatsP_HHG{s},RmaxP_HHG{s}, rhoP_HHG{s}] = ...
            MjgER2016_phasePrecession(Trial,ghzp,ddzp,phz,spkpp,unitSubset,[],[],[],['placefield_summary_ghzp',DataHash(states{s})]);
        [PN_HHG{s}, phzStatsN_HHG{s},RmaxN_HHG{s}, rhoN_HHG{s}] = ...
            MjgER2016_phasePrecession(Trial,ghzn,ddzn,phz,spkpp,unitSubset,[],[],[],['placefield_summary_ghzn',DataHash(states{s})]);
        % GDZ phase precession
        [P_HDG{s},   phzStats_HDG{s}, Rmax_HDG{s}, rho_HDG{s}] = ...
            MjgER2016_phasePrecession(Trial,gdz,ddz,phz,spkpp,unitSubset,[],[],[],['placefield_summary_gdz',DataHash(states{s})]);
        [PP_HDG{s}, phzStatsP_HDG{s},RmaxP_HDG{s}, rhoP_HDG{s}] = ...
            MjgER2016_phasePrecession(Trial,gdzp,ddzp,phz,spkpp,unitSubset,[],[],[],['placefield_summary_gdzp',DataHash(states{s})]);
        [PN_HDG{s}, phzStatsN_HDG{s},RmaxN_HDG{s}, rhoN_HDG{s}] = ...
            MjgER2016_phasePrecession(Trial,gdzn,ddzn,phz,spkpp,unitSubset,[],[],[],['placefield_summary_gdn',DataHash(states{s})]);
    
    end

    
% $$$     mCom = pfstats.peakPatchCOM(:,:,ismember(pfstats.cluMap,unitSubset),:);
% $$$     mCom(mCom==0) = nan;
% $$$     mCom = sq(mean(pfstats.peakPatchCOM,2,'omitnan'));
    
    
    for u = 1:numel(unitSubset),  tic        
        %figure(hfig);
        clf();
        hfig.PaperPositionMode = 'auto';        
        sp = gobjects([1,0]);
        unit = unitSubset(u);
        
% $$$         mpfsRate = max(cell2mat(cf(@(p,u) max(p.maxRate(u)),...
% $$$                                    pfkbs,repmat({unit},[1,numStates]))));
        maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                      [pfs,dfs],repmat({unit},[1,numel(pfs)+numel(dfs)]))));

        pfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean',0.5)),pfs,...
                                  repmat({unit},[1,numStates])));
        dfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean',0.5)),dfs,...
                                  repmat({unit},[1,numel(dfs)])));
        
        if all(pfsMaxRates==0),
            continue;
        end
        
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
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'c';
        hax.EdgeColor = 'c';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;
        ind = [stc{'hloc'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'r';
        hax.EdgeColor = 'r';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;

        % row1col1 - pitch hl pause distrib
        sp(end+1) = subplot2(ny,numStates+2,7:8,1);hold('on');
        ind = [stc{'lpause'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
        hax.FaceColor = 'c';
        hax.EdgeColor = 'c';
        hax.FaceAlpha = 0.5;
        hax.EdgeAlpha = 0.5;
        ind = [stc{'hpause'}];
        hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
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
            res(abs(ddz(res,u))>350) = [];
            if numel(res) >10,
                res(res>xyz.size(1))=[];
                drzspk = drz(res,u);
                hrzspk = hrz(res,u);                
                ghzspk = ghz(res,u);
                gdzspk = gdz(res,u);
                ddzspk = ddz(res,u);
                phzspk = phz(res,1);
                gind = ~isnan(drzspk)&~isnan(phzspk)&~isnan(ghzspk)&~isnan(gdzspk);
            else
                res = [];
                drzspk=[];
                ddzspk=[];
                phzspk=[];
                ghzspk=[];
                gdzspk=[];
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
% $$$             if unit==unitSubset(1),
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
                plot([ghzspk(gind);ghzspk(gind)],...
                      [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
                      '.','MarkerSize',4);
                xlim([-1,1]);
                ylim([-180,540]);
                
                plot([-1,1],circ_rad2ang(P_HHG{s}(u,1,1)*[-1,1]+P_HHG{s}(u,2,1)),'-m','LineWidth',1)
                plot([-1,1],circ_rad2ang(P_HHG{s}(u,1,1)*[-1,1]+P_HHG{s}(u,2,1))+360,'-m','LineWidth',1)

                plot([0,1],circ_rad2ang(PP_HHG{s}(u,1,1)*[0,1]+PP_HHG{s}(u,2,1)),'-r','LineWidth',1)
                plot([0,1],circ_rad2ang(PP_HHG{s}(u,1,1)*[0,1]+PP_HHG{s}(u,2,1))+360,'-r','LineWidth',1)

                plot([-1,0],circ_rad2ang(PN_HHG{s}(u,1,1)*[-1,0]+PN_HHG{s}(u,2,1)),'-g','LineWidth',1)
                plot([-1,0],circ_rad2ang(PN_HHG{s}(u,1,1)*[-1,0]+PN_HHG{s}(u,2,1))+360,'-g','LineWidth',1)

                title(['GHZ rho: ',num2str(round(rho_HHG{s}(u),2))]);                
            end
            sp(end).XTickLabel = {};

            sp(end+1) = subplot2(ny,numStates+2,[9,10],s+1); hold('on');
            if sum(gind)>10,
                plot([hrzspk(gind);hrzspk(gind)],...
                     [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
                     '.','MarkerSize',4);
                xlim([-1,1]);
                ylim([-180,540]);
                
                plot([-1,1],circ_rad2ang(P_H{s}(u,1,1)*[-1,1]+P_H{s}(u,2,1)),'-m','LineWidth',1)
                plot([-1,1],circ_rad2ang(P_H{s}(u,1,1)*[-1,1]+P_H{s}(u,2,1))+360,'-m','LineWidth',1)

                plot([0,1],circ_rad2ang(PP_H{s}(u,1,1)*[0,1]+PP_H{s}(u,2,1)),'-r','LineWidth',1)
                plot([0,1],circ_rad2ang(PP_H{s}(u,1,1)*[0,1]+PP_H{s}(u,2,1))+360,'-r','LineWidth',1)

                plot([-1,0],circ_rad2ang(PN_H{s}(u,1,1)*[-1,0]+PN_H{s}(u,2,1)),'-g','LineWidth',1)
                plot([-1,0],circ_rad2ang(PN_H{s}(u,1,1)*[-1,0]+PN_H{s}(u,2,1))+360,'-g','LineWidth',1)
                xlim([-1,1]),
                ylim([-180,540])
                title(['HRZ rho: ',num2str(round(rho_H{s}(u),2))]);
            end

            
% GET placefield center
% COMPUTE position of the placefield center relative to head basis
            [mxr,mxp] = pft.maxRate(unit);
            pfsCenterHR = MTADfet.encapsulate(Trial,                                               ...
                                              multiprod(bsxfun(@minus,mxp,sq(xyz(:,'nose',[1,2]))),...
                                                        hvec,2,[2,3]),                             ...
                                              sampleRate,                                          ...
                                              'placefield_center_referenced_to_head',              ...
                                              'pfsCenterHR',                                       ...
                                              'p'                                                  ...
            );
            
            
            
            %res = res(WithinRanges(res,get([stc{'x+p&t'}],'data'))&cind(res));
            
% MEAN spike theta phase at location of place field center relative to head (Wait ... what?)
            sp(end+1) = subplot2(ny,numStates+2,[11,12],s+1);
            if sum(gind)>10&sum(sqrt(sum(pfsCenterHR(res,:).^2,2))<350)>10,
                scatter(pfsCenterHR(res,1),pfsCenterHR(res,2),5,...
                        phz(res,1),'filled');
                colormap(sp(end),'hsv');
                set(gca,'Color',[0.75,0.75,0.75]);
                ylim([-300,300]);
                xlim([-300,300]);
                ylabel('frontal');
                ylabel('lateral');
            end

% MEAN spike theta phase at location of place field center relative to head (Wait ... what?)
            sp(end+1) = subplot2(ny,numStates+2,[13,14],s+1);
            if sum(gind)>10&sum(sqrt(sum(pfsCenterHR(res,:).^2,2))<350)>10,
                bins = linspace(-300,300,15);
                hind = discretize(pfsCenterHR(res,:),bins);
                mtph = accumarray(hind(nniz(hind),:),                                ... subs
                                  phz(res(nniz(hind)),1),... vals
                                  repmat(numel(bins),[1,2]),                         ... size
                                   @circ_mean);                                         % func
                mtph(mtph==0) = nan;
                imagescnan({bins,bins,mtph'},[-pi,pi],'circular',true,'colorMap',@hsv);
                axis('xy');
            end

            
% TRANITION triggered histogram onset
% TRANITION triggered histogram offset
% $$$             sp(end+1) = subplot2(ny,numStates+2,[15],s+1); hold('on');            
% $$$             plot(bhvccg{sts}.tbin,RectFilter(sq(bhvccg{s}.ccg(:,unit==bhvccg{s}.cluMap(:,1),:,1))))
% $$$             sp(end+1) = subplot2(ny,numStates+2,[16],s+1); hold('on');            
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
            hax = colorbar(sp(end));
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
        %print(hfig,'-depsc2',fullfile(FigDir,Trial.filebase,[FigName,'.eps']));        
        print(hfig,'-dpng',  fullfile(FigDir,Trial.filebase,[FigName,'.png']));
toc
        delete(sp);        

    end%for unitSubset
end%for Trials


