function report_interneuron_summary(Trial,varargin)
% function report_placefield_summary
% 
% description: plot a variety of measures for units with place fields 
%
% Mod:20171211: reorganize the signature of multiple functions to include units
%               as the second input option.
%
% INFO :
%    spike info
%    spike waveform
%    place ratemaps given state
%    Theta phase preference

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('figDir',                '/storage/gravio/figures/analysis/interneurons',       ...
                 'pitchReferenceTrial',   'Ed05-20140529.ont.all',                               ...
                 'marker',                'hcom',                                                ...
                 'thetaRef',              1,                                                     ...
                 'phzCorrection',         0,                                                     ...
                 'overwrite',             false                                                  ...
);
[figDir,pitchReferenceTrial,marker,thetaRef,phzCorrection,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------


unitsInt   = select_units(Trial,'int');
unitsInt   = remove_bad_units(Trial,unitsInt);

%unitsInt = select_units(Trial,'int');


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


% LOAD Trial 
% CREATE Trial directory in specified location
% PRINT intent of trial processing
Trial.load('nq');

unitSubset = [unitsInt];    
create_directory(fullfile(figDir,Trial.filebase));
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
%pch = fet_HB_pitchB(Trial,sampleRate);

% CREATE lowpass filtered xyz object
% COMPUTE basis vector aligned to the head
fxyz = filter(copy(xyz),'ButFilter',3,20,'low');    
hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);

% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',thetaRef);
phz = lfp.phase([5,13]);    
phz.data = unwrap(phz.data);
phz.resample(xyz);    
%phz.data = mod(phz.data+pi,2*pi)-pi;
phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; 
phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;
lfp.resample(xyz);


% SELECT high quality units with place fields
% COMPUTE placefields' statistics
% LOAD placefields for theta state
% LOAD spikes 

if isempty(unitSubset), return;end;
pft = pfs_2d_theta(Trial,                                          ...
                   unitSubset,                                     ...
                   false,                                          ...
                   'pfsArgsOverride',struct('numIter',1,           ...
                                            'halfsample',false),   ...
                   'overwrite',overwrite);    
%pfstats = compute_pfstats_bs(Trial,unitSubset);
pfstats = [];
spk = Trial.spk.copy();
spk.create(Trial,              ...  
           sampleRate,         ...
           'theta',            ...
           unitSubset,         ...
           '');    

% LOAD placefields and subsampled estimate
pfs = pfs_2d_states(Trial,                                         ...
                    unitSubset,                                    ...
                    [],                                            ...
                    states(1:end-1),                               ...
                    'pfsArgsOverride',struct('numIter',1,          ...
                                             'halfsample',false),  ...
                    'overwrite',overwrite);
pfs{end+1} = pft;

% COMPUTE behavioral state occupancy maps
stsocc = cell(size(states));
for sts = 1:numStates,    
    [stsocc(sts),stsoccBins] = xyocc(Trial,stc{states{sts}});
end
maxocc = max(cellfun(@(r) max(r(:)),stsocc));

% COMPUTE drz and ddz
% $$$ drz = compute_drz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
% $$$ ddz = compute_ddz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
% $$$ hrz = compute_hrz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
% $$$ ghz = compute_ghz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);
% $$$ gdz = compute_gdz(Trial,unitSubset,pft,pfstats,[],marker,[],xyz);    

% COMPUTE unit auto correlogram
[accg,tbins] = autoccg(Trial);



for u = 1:numel(unitSubset),  tic        
    %figure(hfig);
    clf();
    hfig.PaperPositionMode = 'auto';        
    sp = gobjects([1,0]);
    unit = unitSubset(u);
    
    maxPfsRate = max(cell2mat(cf(@(p,u) maxRate(p,u,false,'prctile99',0.5),...
                                 [pfs],repmat({unit},[1,numel(pfs)]))));
    pfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,true,'mean',0.5)),pfs,...
                              repmat({unit},[1,numStates])));
% $$$     dfsMaxRates = cell2mat(cf(@(p,u) max(p.maxRate(u,false,'mean',0.5)),dfs,...
% $$$                               repmat({unit},[1,numel(dfs)])));
    
    if all(pfsMaxRates==0),
        continue;
    end
    
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

    FigInfo = uicontrol('Parent',hfig,                                                           ...
                        'Style','text',                                                          ...
                        'String',{['Unit: ',num2str(unit),                                       ...
                                     ' El:',num2str(Trial.spk.map(unit,2)),                      ...
                                     ' eclu:',num2str(Trial.spk.map(unit,3))],                   ...
                                   [Trial.filebase],                                             ...
                                   ['stcMode: ',Trial.stc.mode],                                 ...
                                   ['eDist:   ',num2str(Trial.nq.eDist(unit))],                  ...
                                   ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],          ...
                                   ['SNR:     ',num2str(Trial.nq.SNR(unit))],                    ...
                                   ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],                 ...
                                   ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]            ...
                                  },...
                        'Units','centimeters',...
                        'Position',[2,18.5,6,6]);

    
    % row1col1 - auto ccg
    sp(end+1) = subplot2(ny,numStates+2,3:4,1);
    bar(tbins,accg(:,unit));axis tight;


    % row1col1 - pitch hl loc distrib
% $$$     sp(end+1) = subplot2(ny,numStates+2,5:6,1);hold('on');
% $$$     ind = [stc{'lloc'}];
% $$$     hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
% $$$     hax.FaceColor = 'c';
% $$$     hax.EdgeColor = 'c';
% $$$     hax.FaceAlpha = 0.5;
% $$$     hax.EdgeAlpha = 0.5;
% $$$     ind = [stc{'hloc'}];
% $$$     hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
% $$$     hax.FaceColor = 'r';
% $$$     hax.EdgeColor = 'r';
% $$$     hax.FaceAlpha = 0.5;
% $$$     hax.EdgeAlpha = 0.5;
% $$$ 
% $$$     % row1col1 - pitch hl pause distrib
% $$$     sp(end+1) = subplot2(ny,numStates+2,7:8,1);hold('on');
% $$$     ind = [stc{'lpause'}];
% $$$     hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
% $$$     hax.FaceColor = 'c';
% $$$     hax.EdgeColor = 'c';
% $$$     hax.FaceAlpha = 0.5;
% $$$     hax.EdgeAlpha = 0.5;
% $$$     ind = [stc{'hpause'}];
% $$$     hax = bar(linspace(-pi/2,pi/2,50),histc(pch(ind,1),linspace(-pi/2,pi/2,50)),'histc');
% $$$     hax.FaceColor = 'r';
% $$$     hax.EdgeColor = 'r';
% $$$     hax.FaceAlpha = 0.5;
% $$$     hax.EdgeAlpha = 0.5;
    
    
    
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
        if numel(res) >10,
            res(res>xyz.size(1))=[];
        else
            continue;
        end
        
        % PLACEFIELDS MTAAknnpfs
        sp(end+1) = subplot2(ny,numStates+2,[1,2],s+1);
        plot(pfs{s},unit,1,'none',mpfsRate,true,[],false,interpPar,@jet);
        set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
        title(sprintf('Max Rate: %3.2f',pfsMaxRates(s)));
        
        
        % TRANSITION Positions
        sp(end+1) = subplot2(ny,numStates+2,[3,4],s+1);        

        hold('on');
        imagesc(stsoccBins{1},stsoccBins{2},stsocc{s}');
        axis('xy');
        set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});                
        caxis([0,maxocc])
        title(statesCcg{s});                            
        
        % Theta phase preference
        sp(end+1) = subplot2(ny,numStates+2,[5,6],s+1);
        rose(phz(res),36);
        %plot(pfs{s},unit,'mean','none',mpfsRate,true,0.5,false,interpPar,@jet);
        %set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
        %title(sprintf('Max Rate: %3.2f',pfsMaxRatesMean(s)));
        
        % TRANITION triggered histogram onset
        % TRANITION triggered histogram offset
% $$$             sp(end+1) = subplot2(ny,numStates+2,[15],s+1); hold('on');            
% $$$             plot(bhvccg{sts}.tbin,RectFilter(sq(bhvccg{s}.ccg(:,unit==bhvccg{s}.cluMap(:,1),:,1))))
% $$$             sp(end+1) = subplot2(ny,numStates+2,[16],s+1); hold('on');            
% $$$              %plot(bhvccg{s}.tbin,RectFilter(sq(bhvccg{s}.ccg(:,unit==bhvccg{s}.cluMap(:,1),:))));
% $$$             ylim([0,mccgRate]);
        
    end

    
    
% $$$     for s = 1:numel(dfs),
% $$$         sp(end+1) = subplot2(ny,numStates+2,[s*2:s*2+1]-1,numStates+2);
% $$$         dfs{s}.plot(unit,'mean',true,mpfdRate,false,0.85,false,[],@jet);
% $$$         hax = colorbar(sp(end));
% $$$         hax.Position(1) = hax.Position(1) + 0.05;
% $$$         title(dfst{s});
% $$$     end
    


    pause(0.01);
    
    FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit,'%04.f')];
    %print(hfig,'-depsc2',fullfile(figDir,Trial.filebase,[FigName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,Trial.filebase,[FigName,'.png']));
    toc
    delete(sp);        

end%for unitSubset

%---------------------------------------------------------------------------------------------------
