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


FigDir = create_directory('/storage/gravio/figures/placefields'); 


% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
             'theta-groom-sit'};

% PROCESS session data
% $$$ cf(@(t)  pfs_2d_theta(t,'overwrite',true),  Trials);
% $$$ cf(@(t)  compute_pfstats_bs(t,'overwrite',true),  Trials);
% $$$ cf(@(t)  MjgER2016_drzfields(t,[],true), Trials);

% FOR each Trial -------------------------------------------------------------

hfig = figure();
hfig.Position = [273, 54, 1269, 681];
hfig.PaperPositionMode = 'auto';
ny = 12;


for t = 1:numel(Trials),
    clf();

% LOAD Trial 
% CREATE Trial directory in specified location
% PRINT intent of trial processing
    Trial = Trials{t};    
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
    tbp_phase = lfp.phase([6,12]);

% SELECT high quality units with place fields
% COMPUTE placefields' statistics
% LOAD placefields for theta state
% LOAD spikes 
    units = select_placefields(Trial);
    if isempty(units), continue;end;
    pft = pfs_2d_theta(Trial,'overwrite',true);    
    pfstats = compute_pfstats_bs(Trial,units);
    spk = Trial.spk.copy();
    spk.create(Trial,xyz.sampleRate,'theta',units,'deburst');    
    
% LOAD placefields and subsampled estimate
    for sts = 1:numel(states),        
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = units;
        defargs.states = states{sts};
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end
    for sts = 1:numel(states),
        defargs = get_default_args('MjgER2016','MTAApfs','struct');
        defargs.units = units;
        defargs.states = states{sts};
        defargs.numIter = 1000;
        defargs.halfsample = true;
        defargs.overwrite = false;
        defargs = struct2varargin(defargs);
        pfs{sts} = MTAApfs(Trial,defargs{:});
    end    
    

% LOAD DRZ fields
    dfs = cell([1,3]);
    [dfs{:}] = MjgER2016_drzfields(Trial,units,false,pfstats);
    dfst = {'pitch','height','rhm'};
    
% COMPUTE drz and ddz
    drz = compute_drz(Trial,units,pft,pfstats);
    ddz = compute_ddz(Trial,units,pft,pfstats);
    
% COMPUTE place fields and subsampled estimate
    for sts = 1:numel(states),
        [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5);
        sper{sts}{1}(sper{sts}{1}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
        sper{sts}{2}(sper{sts}{2}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
    end
    
    [accg,tbins] = autoccg(Trial);
    
    mCom = pfstats.peakPatchCOM(:,:,ismember(pfstats.cluMap,units),:);
    mCom(mCom==0) = nan;
    mCom = sq(mean(pfstats.peakPatchCOM,2,'omitnan'));
    

        
    for u = 1:numel(units),  tic        
        
        %clf();
        sp = [];
        unit = units(u);

        
        mpfsRate = max(cell2mat(cf(@(p,u) max(p.maxRate(u)),...
                                   pfkbs,repmat({unit},[1,numel(bhvccg)]))));
        if mpfsRate<=0,mpfsRate=1;end % default to 1 if 0
        mccgRate = max(cell2mat(cf(@(c,u) max(max(c.ccg(:,u,:))),...
                                   bhvccg,repmat({unit},[1,numel(bhvccg)]))));
        if mccgRate<=0,mccgRate=1;end % default to 1 if 0

        uResw = spkw.res(spkw.clu==unit);
        uSpkw = spkw.spk(spkw.clu==unit,:,:);
        [~,sInd] = SelectPeriods(uResw,[stc{states{end},1}],'d',1,0);
        if numel(sInd)<=1, mspkt = zeros([size(uSpkw,2),size(uSpkw,3)]);
        else,              mspkt = bsxfun(@plus,sq(mean(uSpkw(sInd,:,:)))',fliplr(linspace(1000,8000,size( uSpkw,2))));
        end

        % row1col1 - auto ccg
        sp(end+1) = subplot2(ny,numel(states)+2,1:2,1);
        bar(tbins,accg(:,unit));axis tight;
        title(['Trial: ',Trial.filebase,' Unit: ',num2str(unit)]);

        % row1col1 - pitch hl loc distrib
        sp(end+1) = subplot2(ny,numel(states)+2,3:4,1);hold('on');
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
        sp(end+1) = subplot2(ny,numel(states)+2,5:6,1);hold('on');
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
        sp(end+1) = subplot2(ny,numel(states)+2,[9:12],1);
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

        
        
        for s = 1:numel(states),
% GET current state
% SELECT spikes within current state 
% SKIP if spike count is less than 50
            state = [stc{states{s},spk.sampleRate}];
            res = spk(unit);
            res = res(WithinRanges(res,state.data));
            if numel(res) >50,
                res(res>xyz.size(1))=[];
                drzspk = drz(res,unit==units);
                %ddzspk = sqrt(abs(ddz(res,unit==units))./pi).*sign(ddz(res,unit==units));
                ddzspk = ddz(res,unit==units);
                phzspk = tbp_phase(res,spk.map(spk.map(:,1)==units(u),2));
                gind = ~isnan(drzspk)&~isnan(phzspk);                
            else
                res = [];
                drzspk=[];
                ddzspk=[];
                phzspk=[];
                gind=[];
            end
            

% PLACEFIELDS MTAAknnpfs
            sp(end+1) = subplot2(ny,numel(states)+2,[1,2],s+1);hold('on');
            plot(pfkbs{s},unit,'mean',[],mpfsRate);
            title(statesCcg{s});            
            set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});            
            plot(mCom(s,u,2),mCom(s,u,1),'*m');

% PLACEFIELDS MTAApfs
            sp(end+1) = subplot2(ny,numel(states)+2,[3,4],s+1);hold('on');
            plot(pfs{s},unit,'mean',[],mpfsRate);
            set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});            
            plot(mCom(s,u,2),mCom(s,u,1),'*m');
            
% TRANSITION Positions
            if unit==units(1),
                subplot2(ny,numel(states)+2,[5,6],s+1);
                hold('on');
                plot(xyz(ceil(sper{s}{1}./1250.*xyz.sampleRate),7,1),...
                     xyz(ceil(sper{s}{1}./1250.*xyz.sampleRate),7,2),'*g');
                plot(xyz(ceil(sper{s}{2}./1250.*xyz.sampleRate),7,1),...
                     xyz(ceil(sper{s}{2}./1250.*xyz.sampleRate),7,2),'*r');
                xlim([-500,500]);ylim([-500,500]);
                set(gca,'YTickLabel',{});set(gca,'XTickLabel',{});
            end                
% TRAJECTORIES Positions
% $$$             if unit==units(1),
% $$$                 subplot2(ny,numel(states)+2,[3,4],s+1);
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
            sp(end+1) = subplot2(ny,numel(states)+2,[7,8],s+1); hold('on');
            if sum(gind)>10,
                plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
                plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
                xlim([-1,1]),
                ylim([-180,540])
            end

            sp(end+1) = subplot2(ny,numel(states)+2,[9,10],s+1); hold('on');
            if sum(gind)>10,
                plot(ddzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
                plot(ddzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
                xlim([-300,300]),
                %xlim([-20,20]),                
                %xlim([-90000,90000]),
                ylim([-180,540])
            end
            
            
% TRANITION triggered histogram onset
% $$$             sp(end+1) = subplot2(ny,numel(states)+2,[5,6],s+1);
% $$$             plot(bhvccg{s},unit,1);axis('tight');        title('onset');
% $$$             ylim([0,mccgRate]);
% $$$             
% TRANITION triggered histogram offset
% $$$             sp(end+1) = subplot2(ny,numel(states)+2,[7,8],s+1);
% $$$             plot(bhvccg{s},unit,2);axis('tight');        title('offset');
% $$$             ylim([0,mccgRate]);

% WAVEFORM of unit
% $$$             sp(end+1) = subplot2(ny,numel(states)+2,[9:12],s+1);
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
            sp(end+1) = subplot2(ny,numel(states)+2,[s*2:s*2+1]-1,numel(states)+2);
            dfs{s}.plot(unit,'maxRate',mpfsRate,'isCircular',false);
            hax = colorbar();
            hax.Position(1) = hax.Position(1) + 0.05;
            title(dfst{s});
        end
        


        pause(0.01);
        
        FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit)];
        %print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(FigDir,Trial.filebase,[FigName,'.png']));
toc
        delete(sp);        

    end%for units
end%for Trials



% TESTING CODE


% $$$     if testRequested,
% $$$     htfig = figure();        
% $$$         for u = 1:numel(units),
% $$$             subplot(1,numel(pfkbs)+1,1);
% $$$             plot(pft,units(u));
% $$$             for sts = 1:numel(pfkbs),
% $$$                 subplot(1,numel(pfkbs)+1,sts+1);
% $$$                 plot(pfkbs{sts},units(u),'mean');
% $$$             end
% $$$             waitforbuttonpress();
% $$$         end
% $$$     end
% $$$ 
