function report_placefield_summary(Trial,varargin)
% function report_placefield_summary
% 
% description: plot a variety of measures for units with place fields 
%
% Mod:20171211: reorganize the signature of multiple functions to include units
%               as the second input option.
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('figDir',                '/storage/gravio/figures/analysis/placefields_new',    ...
                 'pitchReferenceTrial',   'Ed05-20140529.ont.all',                               ...
                 'marker',                'hcom',                                                ...
                 'overwrite',             false                                                  ...
);
[figDir,pitchReferenceTrial,marker,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

unitsPyr = Trial.spk.get_unit_set(Trial,'placecells');
%unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');
unitsInt = [];

states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};

numStates = numel(states);


pfdVersion = '9';


hfig = figure(666001);
pause(0.1);
hfig.Position = [1, 1, 2482, 1274];
ny = 15;
sampleRate = 250;

Trial.load('nq');

unitSubset = [unitsPyr];

create_directory(fullfile(figDir,Trial.filebase));
disp(['[INFO] MTA:analysis:report_placefield_summary: Processing Trial - ' Trial.filebase]);


% LOAD the state collection
stc = Trial.stc.copy();
% LOAD the spike waveforms
spkw = Trial.spk.copy();
spkw.load_spk(Trial);
% LOAD the marker positions
xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
% COMPUTE the head pitch and
% COREGISTER the current head pitch to a reference session
% USING the expected pitch given the head and body speed 
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
phz = load_theta_phase(Trial,xyz);

if isempty(unitSubset)
    return;
end

pft = pfs_2d_theta(Trial, unitSubset);    
%pfstats = compute_pfstats_bs(Trial,unitSubset);
pfstats = [];
spk = Trial.spk.copy();
spk.create(Trial,              ...  
           xyz.sampleRate,     ...
           'theta',            ...
           unitSubset,         ...
           '');    

% LOAD placefields and subsampled estimate
pfs = pfs_2d_states(Trial, unitSubset, [], states, 'overwrite', false);

% COMPUTE behavioral state occupancy maps
stsocc = cell(size(states));
for sts = 1:numStates,    
    [stsocc(sts),stsoccBins] = xyocc(Trial,stc{states{sts}});
end
maxocc = max(cellfun(@(r) max(r(:)),stsocc));

% COMPUTE drz and ddz
ddz = compute_ddz(Trial,unitSubset,pft,[],[],marker,[],xyz,[],[],sampleRate);
hrz = compute_hrz(Trial,unitSubset,pft,[],[],marker,[],xyz,[],[],sampleRate);
ghz = compute_ghz(Trial,unitSubset,pft,[],[],marker,[],xyz,[],[],sampleRate);


% COMPUTE unit auto correlogram
[accg,tbins] = autoccg(Trial);

% LOAD DRZ fields
dfs = {compute_bhv_ratemaps(Trial,unitSubset,[],[],pft,[],0.8,250,'overwrite',overwrite,'purge',true)};
dfst = {'HPITCHxBPITCH'};
mask_bhv =load('/storage/gravio/data/project/general/analysis/pfsHB_mask.mat');

% COMPUTE phase precession
ghzp = ghz;    ghzp(ghzp<0)=nan;
ghzn = ghz;    ghzn(ghzn>0)=nan;
P_HHG  = {};   phzStats_HHG  = {};   Rmax_HHG  = {};
PP_HHG = {};   phzStatsP_HHG = {};   RmaxP_HHG = {};
PN_HHG = {};   phzStatsN_HHG = {};   RmaxN_HHG = {};
rho_HHG= {};

for s = 1:numStates,
    spkpp = Trial.spk.copy();
    spkpp.create(Trial,xyz.sampleRate,states{s},unitSubset,'deburst');
    % GHZ phase precession
    [P_HHG{s},   phzStats_HHG{s}, Rmax_HHG{s}, rho_HHG{s}] = ...
        MjgER2016_phasePrecession(                           ...
            Trial,                                           ...
            ghz,                                             ...
            ddz,                                             ...
            phz,                                             ...
            spkpp,                                           ...
            unitSubset,                                      ...
            [],[],[],                                        ...
            ['placefield_summary_ghz',DataHash(states{s})],  ...
            true);
end


for u = 1:numel(unitSubset),  tic        
    clf();
    hfig.PaperPositionMode = 'auto';        
    sp = gobjects([1,0]);
    unit = unitSubset(u);
    maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                 [pfs,dfs],repmat({unit},[1,numel(pfs)+numel(dfs)]))));

    pfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean',0.5)),pfs,...
                              repmat({unit},[1,numStates])));
    dfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean',0.5)),dfs,...
                              repmat({unit},[1,numel(dfs)])));
    
    if all(pfsMaxRates==0),
        continue;
    end
    
    mpfsRate = maxPfsRate;        
    mpfdRate = maxPfsRate;
    
    pfsMaxRatesMean = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean')),pfs,repmat({unit},[1,numStates])));

% GET spke waveform of the current unit
    uResw = spkw.res(spkw.clu==unit);
    uSpkw = spkw.spk(spkw.clu==unit,:,:);
    [~,sInd] = SelectPeriods(uResw,[stc{states{end},1}],'d',1,0);
    if numel(sInd)<=1, mspkt = zeros([size(uSpkw,2),size(uSpkw,3)]);
    else,              mspkt = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
    end

% PLOT unit quality information
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
                        'Position',[8,22.5,6,3.5]);

    
% PLOT auto-cross-correleogram
    sp(end+1) = subplot2(ny,numStates+2,6:7,1);
    bar(tbins,accg(:,unit));axis tight;

    
% PLOT ave spike waveform
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
    sp(end).YTick = [];
    sp(end).XTick = [];

    
    
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
            ghzspk = ghz(res,u);
            ddzspk = ddz(res,u);
            phzspk = phz(res,1);
            gind = ~isnan(ddzspk)&~isnan(phzspk)&~isnan(ghzspk);
        else
            res = [];
            ddzspk=[];
            phzspk=[];
            ghzspk=[];
            gind=[];
        end
        
        
        % PLACEFIELDS MTAApfs
        sp(end+1) = subplot2(ny,numStates+2,[1,2],s+1);
        sp(end).Position(4) = 0.14;
        plot(pfs{s},unit,'mean','none',mpfsRate,true,0.5,false,[],@jet);
        set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
        title(sprintf('Max Rate: %3.2f',pfsMaxRatesMean(s)));
        
        % OCCUPANCY Positions
        sp(end+1) = subplot2(ny,numStates+2,[4,5],s+1);
        sp(end).Position(4) = 0.14;
        hold('on');
        imagesc(stsoccBins{1},stsoccBins{2},log10(stsocc{s}'));
        axis('xy');
        set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});                
        caxis([-1,log10(maxocc)]);
        circle(0,0,450,'g')
        colormap('jet');
        title(states{s});                            

        % PHASE PRECESSION 
        sp(end+1) = subplot2(ny,numStates+2,[7,8],s+1); hold('on');
        sp(end).Position(4) = 0.14;
        if sum(gind)>10,
            plot([ghzspk(gind);ghzspk(gind)],...
                 [circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360],...
                 '.','MarkerSize',4);
            xlim([-1,1]);
            ylim([-180,540]);
            plot([-1,1],circ_rad2ang(P_HHG{s}(u,1,1)*[-1,1]+P_HHG{s}(u,2,1)),'-m','LineWidth',1);
            plot([-1,1],circ_rad2ang(P_HHG{s}(u,1,1)*[-1,1]+P_HHG{s}(u,2,1))+360,'-m','LineWidth',1);
            title(['GHZ rho: ',num2str(round(rho_HHG{s}(u),2))]);                
        end
        sp(end).XTickLabel = {};

        % GET placefield center
        % COMPUTE position of the placefield center relative to head basis
        [mxr,mxp] = pft.maxRate(unit);
        pfsCenterHR = MTADfet.encapsulate( ...
            Trial,                                               ...
            multiprod(bsxfun(@minus,mxp,sq(xyz(:,'nose',[1,2]))),...
                      hvec,2,[2,3]),                             ...
            sampleRate,                                          ...
            'placefield_center_referenced_to_head',              ...
            'pfsCenterHR',                                       ...
            'p'                                                  ...
        );
        
        % MEAN spike theta phase at location of place field center relative to head
        sp(end+1) = subplot2(ny,numStates+2,[11,12],s+1);
        sp(end).Position(4) = 0.14;
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
    end

    for s = 1:numel(dfs),
        sp(end+1) = subplot2(ny,numStates+2,[s*2:s*2+1]-1,numStates+2);
        sp(end).Position(4) = 0.14;
        dfs{s}.plot(unit,1,'colorbar',[0,mpfdRate],true,[],false,[],@jet,mask_bhv.mask);
        hax = colorbar(sp(end));
        hax.Position(1) = hax.Position(1) + 0.05;
        title(dfst);
    end

    pause(0.01);
    
    FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit,'%04.f')];
    print(hfig,'-dpng',  fullfile(figDir,Trial.filebase,[FigName,'.png']));
    toc
    delete(sp);        
end%for unitSubset

%---------------------------------------------------------------------------------------------------
