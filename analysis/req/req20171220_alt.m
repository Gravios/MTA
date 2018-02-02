

sessionList = get_session_list('MjgER2016');

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
units = cf(@(t)  select_placefields(t),  Trials);
sessionList(cellfun(@isempty,units)) = [];
Trials(cellfun(@isempty,units)) = [];
units(cellfun(@isempty,units)) = [];


states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesLabels = {'loc','lloc','hloc','rear',     ...
          'pause','lpause','hpause',            ...
          'theta-groom-sit'};


numUnits = sum(cell2mat(cf(@(u) numel(u)  ,units)));
numStates = numel(states);


P         = cell([numel(sessionList),numStates]);
phzStats  = cell([numel(sessionList),numStates]);
Rmax      = cell([numel(sessionList),numStates]);
drzHCnt   = cell([numel(sessionList),numStates]);

PP        = cell([numel(sessionList),numStates]);
phzStatsP = cell([numel(sessionList),numStates]);
RmaxP     = cell([numel(sessionList),numStates]);
drzHCntP  = cell([numel(sessionList),numStates]);

gPN       = cell([numel(sessionList),numStates]);
phzStatsN = cell([numel(sessionList),numStates]);
RmaxN     = cell([numel(sessionList),numStates]);
drzHCntN  = cell([numel(sessionList),numStates]);
         
for t = 17:numel(sessionList),
    Trial = MTATrial.validate(sessionList(t));
    stc = Trial.stc.copy();
    xyz = preproc_xyz(Trial);
    pft = pfs_2d_theta(Trial);
    drz = compute_drz(Trial, units{t}, pft);
    ddz = compute_ddz(Trial, units{t}, pft);
    %lfp = Trial.load('lfp',[Trial.ephys.electrode(:).CA1pyrThetaRef]);
    lfp = Trial.load('lfp',sessionList(t).thetaRef);
    lfp.resample(xyz);
    phz = lfp.phase([6,12]);

    drzp = drz;    drzp(drzp<0)=nan;
    ddzp = ddz;    ddzp(ddzp<0)=nan;
    drzn = drz;    drzn(drzn>0)=nan;
    ddzn = ddz;    ddzn(ddzn>0)=nan;
    
    for s = 1:numStates,
        spk = Trial.spk.copy();
        spk.create(Trial,xyz.sampleRate,states{s},units{t},'deburst');
        [P{t,s},phzStats{t,s},Rmax{t,s},drzHCnt{t,s}] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units{t});
        [PP{t,s},phzStatsP{t,s},RmaxP{t,s},drzHCntP{t,s}] = MjgER2016_phasePrecession(Trial,drzp,ddzp,phz,spk,units{t});
        [PN{t,s},phzStatsN{t,s},RmaxN{t,s},drzHCntN{t,s}] = MjgER2016_phasePrecession(Trial,drzn,ddzn,phz,spk,units{t});
    end
end


pa = reshape(permute(cat(1,P{:}),[1,4,2,3]),numUnits,numStates,size(P{1},2),size(P{1},3));
ra = reshape(permute(cat(1,Rmax{:}),[1,3,2]),numUnits,numStates,size(Rmax{1},2));
da = reshape(permute(cat(1,drzHCnt{:}),[1,3,2]),numUnits,numStates,size(drzHCnt{1},2));

%da = bsxfun(@rdivide,da,sum(da,3));
dax = sum(da>10,3)>8;

pan = reshape(permute(cat(1,PN{:}),[1,4,2,3]),numUnits,numStates,size(PN{1},2),size(PN{1},3));
ran = reshape(permute(cat(1,RmaxN{:}),[1,3,2]),numUnits,numStates,size(RmaxN{1},2));
dan = reshape(permute(cat(1,drzHCntN{:}),[1,3,2]),numUnits,numStates,size(drzHCntN{1},2));
danx = sum(dan>10,3)>4;

pap = reshape(permute(cat(1,PP{:}),[1,4,2,3]),numUnits,numStates,size(PP{1},2),size(PP{1},3));
rap = reshape(permute(cat(1,RmaxP{:}),[1,3,2]),numUnits,numStates,size(RmaxP{1},2));
dap = reshape(permute(cat(1,drzHCntP{:}),[1,3,2]),numUnits,numStates,size(drzHCntP{1},2));
dapx = sum(dap>10,3)>4;


% Plot JPDF of linear-circular-resultant
FigDir = create_directory('/storage/gravio/figures/placefields_phase_precession');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hax = gobjects([numStates,3]);
for s = 1:numStates,
    ind = ':';
    ind = dax(:,s);
    hax(s,1) = subplot2(numStates,3,s,1);    
    g = histogram2(ra(ind,s,1),pa(ind,s,1,1),...
                   linspace(0,1,25),...
                   linspace(-1,1,25),...
                   'DisplayStyle','tile'); 
    g.EdgeAlpha = 0;
    if s==round(numStates/2), ylabel('linear-circular r'); end    
    hax(s,1).XTickMode = 'manual';    
    hax(s,1).XTick = [0,0.25,0.5,0.75,1];    
    hax(s,1).YTickMode = 'manual';    
    hax(s,1).YTick = [-3,-2,-1,0,1,2,3];
    if s~=numStates,  hax(s,1).XTickLabel = {}; end
    if s==1, title({'DRZ [-1,1]'}); end
    ylabel(statesLabels{s});
    
    hax(s,2) = subplot2(numStates,3,s,2);
    g = histogram2(ran(ind,s,1),pan(ind,s,1,1),...
                   linspace(0,1,25),...
                   linspace(-1,1,25),...
                   'DisplayStyle','tile');
    g.EdgeAlpha = 0;
    if s==numStates, xlabel('mean resultant length'); end
    hax(s,2).XTickMode = 'manual';    
    hax(s,2).XTick = [0,0.25,0.5,0.75,1];        
    if s~=numStates, hax(s,2).XTickLabel = {};  end    
    hax(s,2).YTickMode = 'manual';    
    hax(s,2).YTick = [-3,-2,-1,0,1,2,3];    
    hax(s,2).YTickLabel = {};    
    if s==1, title({'DRZ [-1,0]'}); end    

    hax(s,3) = subplot2(numStates,3,s,3);
    g = histogram2(rap(ind,s,1),pap(ind,s,1,1),...
                   linspace(0,1,25),...
                   linspace(-1,1,25),...
                   'DisplayStyle','tile');
    g.EdgeAlpha = 0;
    hax(s,3).XTickMode = 'manual';    
    hax(s,3).XTick = [0,0.25,0.5,0.75,1];            
    if s~=numStates, 
        hax(s,3).XTickLabel = {};
    end
    hax(s,3).YTickMode = 'manual';    
    hax(s,3).YTick = [-3,-2,-1,0,1,2,3];    
    hax(s,3).YTickLabel = {};    
    if s==1, title({'DRZ [0,1]'}); end        

end
hfig.Position = [0.5,0.5,12,40];
suptitle('Phase Precession Slope');
af(@(h)  set(h,'Units','centimeters'),  hax);
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax);
FigName = ['phasePrecession_resLength_X_linearCircularSlope'];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));



figure();
for s = 1:numStates,
    subplot2(numStates,3,s,1);
    g = histogram2(ra(:,s,1),pa(:,s,1,2),linspace(0,1,25),linspace(-pi,pi,25),'DisplayStyle','tile'); g.EdgeAlpha = 0;
    subplot2(numStates,3,s,2);
    g = histogram2(ran(:,s,1),pan(:,s,1,2),linspace(0,1,25),linspace(-pi,pi,25),'DisplayStyle','tile');g.EdgeAlpha = 0;
    subplot2(numStates,3,s,3);
    g = histogram2(rap(:,s,1),pap(:,s,1,2),linspace(0,1,25),linspace(-pi,pi,25),'DisplayStyle','tile');g.EdgeAlpha = 0;
end


figure();
plot(pan(:,s,1,1),pa(:,s,1,1),'.')


% Fit
figure();
plot(fitRange,R);   xlim([-pi,pi]); hold('on');

figure();
t = 1
s = 2;
u = 3;
for u = 1:numel(units),
    clf();
    unit = units(u);



    subplot(131);
    plot(pft,unit,'mean');
    hold('on');
    plot(mrp(u,1),mrp(u,2),'*m');
    title(num2str(unit));


    subplot(132);
    hold('on');
    state = [stc{states{s},spk.sampleRate}];
    res = spk(unit);
    res = res(WithinRanges(res,state.data));
    if numel(res) >10,
        res(res>xyz.size(1))=[];
        drzspk = drz(res,u);
        phzspk = phz(res,spk.map(spk.map(:,1)==units(u),2));
        gind = ~isnan(drzspk)&~isnan(phzspk);                
    else
        res = [];
        drzspk=[];
        ddzspk=[];
        phzspk=[];
        gind=[];
    end
    plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'b.');
    plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'b.');
    xlim([-1,1]),
    ylim([-180,540])
    plot([-1,1],circ_rad2ang(2*pi*P{t,s}(u,1,1)*[-1,1]+P{t,s}(u,1,2)),'-m','LineWidth',2)
    plot([-1,1],circ_rad2ang(2*pi*P{t,s}(u,1,1)*[-1,1]+P{t,s}(u,1,2))+360,'-m','LineWidth',2)

    plot([0,1],circ_rad2ang(2*pi*PP{t,s}(u,1,1)*[0,1]+PP{t,s}(u,1,2)),'-r','LineWidth',2)
    plot([0,1],circ_rad2ang(2*pi*PP{t,s}(u,1,1)*[0,1]+PP{t,s}(u,1,2))+360,'-r','LineWidth',2)

    plot([-1,0],circ_rad2ang(2*pi*PN{t,s}(u,1,1)*[-1,0]+PN{t,s}(u,1,2)),'-g','LineWidth',2)
    plot([-1,0],circ_rad2ang(2*pi*PN{t,s}(u,1,1)*[-1,0]+PN{t,s}(u,1,2))+360,'-g','LineWidth',2)

    waitforbuttonpress();
end



