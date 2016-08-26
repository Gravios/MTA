MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;
trialList = 'Ed10VR_teleport';
OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';

T = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');
T(3).offsets = [15,-90];

if overwriteSession,
    Session = MTASession(T(1).sessionName,  ...
                         T(1).mazeName,     ...
                         true,              ...
                         T(1).TTLValue,     ...
                         'vicon',           ...
                         'nlx',             ...
                         T(1).xyzSampleRate ...
    );

    xyz = Session.load('xyz');
    xyz.data(:,:,1) = xyz.data(:,:,1)+T(1).xOffSet;
    xyz.data(:,:,2) = xyz.data(:,:,2)-T(1).yOffSet;
    xyz.save;
end

if overwriteTrials
    QuickTrialSetup(T,'overwrite',true); 
end


Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);


%% Basic Threshold BHV and LFP Labeling

if overwriteStc
    if isempty(Trial.stc.gsi('v')),
        xyz = Trial.load('xyz');
        xyz.filter('ButFilter',3,2.4);
        fvxy = xyz.vel(1,[1,2]);
        fvxy.data(fvxy.data<1e-3)=1e-3;
        fvxy.data = log10(fvxy.data);
        vper = ThreshCross(fvxy.data,0.5,round(.25*xyz.sampleRate));
        Trial.stc.addState(Trial.spath,...
                           Trial.filebase,...
                           vper,...
                           xyz.sampleRate,...
                           Trial.sync.copy,...
                           Trial.sync.data(1),...
                           'velthresh','v');
    end

    if isempty(Trial.stc.gsi('h')),
        xyz = Trial.load('xyz');
        xyz.filter('ButFilter',3,2.4);
        fvxy = xyz.vel(6,[1,2]);
        fvxy.data(fvxy.data<1e-3)=1e-3;
        fvxy.data = log10(fvxy.data);
        vper = ThreshCross(fvxy.data,0.5,round(.25*xyz.sampleRate));
        Trial.stc.addState(Trial.spath,...
                           Trial.filebase,...
                           vper,...
                           xyz.sampleRate,...
                           Trial.sync.copy,...
                           Trial.sync.data(1),...
                           'velHthresh','h');
    end

    if isempty(Trial.stc.gsi('r')),
        rper = rear(Trial,'com',45);
        Trial.stc.addState(Trial.spath,...
                           Trial.filebase,...
                           rper,...
                           xyz.sampleRate,...
                           Trial.sync.copy,...
                           Trial.sync.data(1),...
                           'rear','r');
    end


    if isempty(Trial.stc.gsi('n')),
        Trial.stc.states{end+1} = Trial.stc{'v'}-(Trial.stc{'r',120}+[-.5,.5]);
        Trial.stc.states{end}.key = 'n'; 
        Trial.stc.states{end}.label = 'NRvel';    
        Trial.stc.states{end}.updateFilename([Trial.filebase,'.sst.',...
                            Trial.stc.states{end}.label,'.',...
                            Trial.stc.states{end}.key,'.mat']);
    end

    if isempty(Trial.stc.gsi('t')),
        Trial = labelTheta(Trial,[],32);
    end
    events = LoadEvents(fullfile(Trial.spath, [Trial.name '.all.evt']));

    
    Trial.stc.save(1);
end

% Calculate and plot
Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

overwrite = false;
units = [];

% Generate unit auto correlogram
[accg,tbin] = autoccg(Trial,units,'theta');




%% Gererate unit rate maps (Place Fields)
binDims = [40,40];
smoothingWeights = [1.2,1.2];
%numIter = 100;
%bootstrap = 0.9;
numIter = 1;
bootstrap = 0;
pfs = {};
for t = 1:nt-1
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial); 
   for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite, ...
                           'binDims',binDims,...
                           'SmoothingWeights',smoothingWeights,...
                           'numIter',numIter,...
                           'bootstrap',bootstrap);
    end
end


t = nt;
Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
Trial.stc = Stc.copy;
Trial.stc.load(Trial); 
for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,[states{i},'-shifted'],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'numIter',numIter,...
                       'bootstrap',bootstrap);

end
t = nt+1;
for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,['shifted&',states{i}],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'numIter',numIter,...
                       'bootstrap',bootstrap);
end






%% Select Units with Firing rate greater than 3Hz in all.theta
t = 1;
mRate = [];
for t = 1:nt;
    for i = 1:nsts,

        mRate(t,i,:) = pfs{t,i}.maxRate;
        
    end
end
units = find(sq(mRate(1,1,:))>3);




%% Gather Place field statistics
pfstats = {};
pfshuff = {};
% Test this version should be able to run multiple units at once
for u = 1:numel(units)
    for t = 1:nt
        Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
        Trial.stc = Stc.copy;
        Trial.stc.load(Trial);
        for i = 1:nsts,        
            [pfstats{t,i,u},pfshuff{t,i,u}] = PlaceFieldStats(Trial,pfs{t,i},units(u));
        end
    end
end
for u = 1:numel(units)
    for i = 1:nsts,        
        [pfstats{t+1,i,u},pfshuff{t+1,i,u}] = PlaceFieldStats(Trial,pfs{t+1,i},units(u));
    end
end
for t = 1:nt+1
    for i = 1:nsts,        
        pcom = cellfun(@(x) sq(x.patchCOM(1,1,find(max(x.patchPFR)==x.patchPFR),:)),pfstats(t,i,:),'UniformOutput',false);
        pind = ~cellfun(@isempty,pcom);
        peakPatchCOM(t,i,pind,:) = sq(cell2mat(pcom(pind)))';
        peakPatchRate(t,i,:) = sq(cellfun(@(x) max(x.patchPFR),pfstats(1,1,:)));
        parea = cellfun(@(x) sq(x.patchArea(1,1,find(max(x.patchPFR)==x.patchPFR),:)),pfstats(t,i,:),'UniformOutput',false);
        pind = ~cellfun(@isempty,parea);
        peakPatchArea(t,i,pind) = sq(cell2mat(parea(pind)))';
    end
end





%% Plot Unit spatial rate maps
spOpts.width  = 4;
spOpts.height = 2;
spOpts.ny = numel(T)+1;
spOpts.nx = numel(states);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';n
figOpts.headerPadding = 2;
figOpts.footerPadding = 2;
figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                     (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];


mkdir(fullfile(OwnDir,'Ed10-20140820-shift_teleport_pfs_bcx'));
sp = [];
autoincr = false;
%autoincr = true;

figHnum = 666999;
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',figOpts.position)
set(hfig,'PaperPositionMode','auto');

unit = units(1);
while unit~=-1,
    clf
    for t = 2:nt+1,
        for i = 1:nsts,
            sp(t,i) = axes('Units',spOpts.units,...
                           'Position',[(spOpts.width +round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                                spOpts.width,...
                                spOpts.height]...
                           );
            hold('on')
            pf = pfs{t,i};

            ratemap = pf.plot(unit,'isCircular',false);
            ratemap(isnan(ratemap)) = -1;
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
            text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
                 sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
            colormap([0,0,0;parula]);
            caxis([-1,sq(mRate(1,i,unit)).*1.5]);
            plot(peakPatchCOM(t,i,unit==units,1),...
                 peakPatchCOM(t,i,unit==units,2),'*k');
            xlim([-600,600]),ylim([-350,350])                    
            title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);
        end
    end

    t = t+1;
    i = 1;
    sp(t,i) = axes('Units',spOpts.units,...
                   'Position',[(spOpts.width +round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                        (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                        spOpts.width,...
                        spOpts.height]...
                   );
    bar(tbin,accg(:,unit));
    xlim([min(tbin),max(tbin)]);
    title([' AutoCCG: Unit ',num2str(unit)]);
% $$$ 
% $$$     print(gcf,'-depsc2',fullfile(OwnDir,'Ed10-20140820-shift_teleport_pfs',...
% $$$                                  ['pfs_',num2str(unit),'.eps']));
% $$$     print(gcf,'-dpng',  fullfile(OwnDir,'Ed10-20140820-shift_teleport_pfs',...
% $$$                                  ['pfs_',num2str(unit),'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end





%% Collect shift information
for u = 1:numel(units)
    for i = 1:nsts,            
        for t = 1:nt
            for j = 1:nt
            pfShift(t,j,i,u,:) = peakPatchCOM(t,i,u,:)-peakPatchCOM(j,i,u,:);
            end
        end        
    end
end
t = t+1;
for u = 1:numel(units)
    for i = 1:nsts,            
        for j = 1:nt
            pfShift(t,j,i,u,:) = peakPatchCOM(t,i,u,:)-peakPatchCOM(j,i,u,:);
        end        
    end
end
t = t-1;
j = j+1;
for u = 1:numel(units)
    for i = 1:nsts,            
        for t = 1:nt
            pfShift(t,j,i,u,:) = peakPatchCOM(t,i,u,:)-peakPatchCOM(j,i,u,:);
        end        
    end
end
t = t+1;
for u = 1:numel(units)
    for i = 1:nsts,            
        pfShift(t,j,i,u,:) = peakPatchCOM(t,i,u,:)-peakPatchCOM(j,i,u,:);
    end
end

t = 1;
i = 1;

spCohere = pfs{t,i}.spatialCoherence(units);


figure,hold on,
ind =spCohere>0.99;
subplot(131)
plot(sq(pfShift(2,3,1,ind))/10,sq(pfShift(2,3,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
xlim([-100,100]),ylim([-100,100])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');
subplot(132)
plot(sq(pfShift(4,5,1,ind))/10,sq(pfShift(4,5,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
xlim([-100,100]),ylim([-100,100])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');
subplot(133)
plot(sq(pfShift(6,7,1,ind))/10,sq(pfShift(6,7,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
xlim([-100,100]),ylim([-100,100])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');


figure
boxplot([sq(pfShift(2,3,1,ind))/10;...
         sq(pfShift(4,5,1,ind))/10;...
         sq(pfShift(6,7,1,ind))/10],...
        [ones([sum(ind),1]);...
         ones([sum(ind),1]).*2;...
         ones([sum(ind),1]).*3]);

figure,hold on,
plot(spCohere,sq(pfShift(2,3,1,:))/10,'.')
plot(spCohere,sq(pfShift(4,5,1,:))/10,'.g')
plot(spCohere,sq(pfShift(6,7,1,:))/10,'.m')

Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');
figure,
subplot(311);hist(sq(pfShift(2,3,1,log10(sq(peakPatchArea(2,1,:)))<4.9))/10,50),Lines(20,[],'r');
title('Fst. control -> shift')
subplot(312);hist(sq(pfShift(4,5,1,log10(sq(peakPatchArea(2,1,:)))<4.9))/10,50),Lines(20,[],'r');
title('Snd. control -> shift')
subplot(313);hist(sq(pfShift(6,7,1,log10(sq(peakPatchArea(2,1,:)))<4.9))/10,50),Lines(20,[],'r');
title('Alt. control -> shift')
ForAllSubplots('xlim([-120,120])')




%% time X xyz X ufr

Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);

xyz = Trial.load('xyz');
xyz.resample(10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'theta',units,0.5,true);

tper = Stc{'t'};
tper.cast('TimeSeries',xyz);
ind = logical(tper.data);

sper = Stc{'x',1};

t =[1:xyz.size(1)]./10;

set(0,'defaultAxesFontSize',10,...
      'defaultTextFontSize',10)

figHnum = 666994;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[11 0 32 25])
set(hfig,'PaperPositionMode','auto');

try,mkdir(fullfile(OwnDir,'/Ed10-20140820-shift_teleport_ufr_pfs_scatter_bcx'));end

autoincr = true;
unit = units(2);
while unit~=-1,
    clf
    subplot2(3,6,[1,2],[1:6]);
    hold on
    scatter(t(~ind),xyz(~ind,6,1),20,ones([sum(~ind),3]).*0.75,'filled');
    scatter(t(ind),xyz(ind,6,1),20,ufr(ind,unit==units),'filled');
    Lines(round((Trial.sync.data(:)-Trial.sync.data(1)))+1,[],'m');
    Lines(sper.data(:,1),[],'g');
    Lines(sper.data(:,2),[],'b');
    ylim([-600,600])    
    title(['Unit(',num2str(unit),') Firing Rate overlayed on Ocuupancy of X axis'])
    i = 1;
    for tr = 1:nt,
        pf = pfs{tr+1,i};
        subplot2(3,6,3,tr);
        hold('on')
        ratemap = pf.plot(unit,'isCircular',false);
        minv = min(ratemap(:));
        maxv = sq(mRate(1,i,unit));
        ratemap(isnan(ratemap)) = minv-((maxv-minv)/20);
        imagesc(pf.adata.bins{2},pf.adata.bins{1},ratemap);    
        text(pf.adata.bins{2}(1)+30,pf.adata.bins{1}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        ca = caxis;
        caxis([ca(1),sq(mRate(1,i,unit)).*1.5]);
        title([pf.session.trialName ':' pf.parameters.states]);
        plot(peakPatchCOM(tr+1,i,unit==units,2),...
             peakPatchCOM(tr+1,i,unit==units,1),'*k');
        ylim([-600,600]),xlim([-350,350])        
        %if tr~=1, set(gca,'YTickLabel',{}),end
    end
    print(gcf,'-depsc2',fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_scatter',...
                                 ['ufr_pfs_',num2str(unit),'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_scatter',...
                                 ['ufr_pfs_',num2str(unit),'.png']));
    unit = figure_controls(hfig,unit,units,autoincr);
end


xoc(:,f,c) = accumarray(vbs(nniz(vbs)&nniz(vys)),vys(nniz(vbs)&nniz(vys),f,c,c),[vbins,1],@nanmean);
uoc(:,f,c) = accumarray(vbs(nniz(vbs)&nniz(vys)),vys(nniz(vbs)&nniz(vys),f,c,c),[vbins,1],@nanmean);



xp = {}
for t = 1:nt
    Trial = MTATrial(T(t+1).sessionName,T(t+1).mazeName,T(t+1).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial);

    txyz = Trial.load('xyz');
    
    tufr = Trial.ufr.copy;
    tufr = tufr.create(Trial,txyz,'theta',units,1.5,true);    

    for u = 1:size(ufr,2)    
        uper = ThreshCross(tufr(:,u),3,2);
        xp{t,u} = [];
        for up = uper'
            xp{t,u}(end+1) = mean(txyz(up(1):up(2),1,1));
        end
    end
end 


mkdir(fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_xtrj_mpos_bcx'))
figHnum = 66;
hfig = figure(figHnum);clf


autoincr = false;
unit =units(1);
i = 1;
while unit~=-1,
    for t = 1:5,

        subplot2(5,2,t,1);            
        hist(xp{t,unit==units},linspace(-600,600,100));
        xlim([-600,600])

        subplot2(5,2,t,2);   
        pf = pfs{t+1,i};

        hold on
        ratemap = pf.plot(unit,'isCircular',false);
        minv = min(ratemap(:));
        maxv = sq(mRate(1,i,unit));
        ratemap(isnan(ratemap)) = minv-((maxv-minv)/20);
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
        text(pf.adata.bins{1}(1)+950,pf.adata.bins{2}(end)-70,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        ca = caxis;a
        caxis([ca(1),sq(mRate(1,i,unit)).*1.5]);
        plot(peakPatchCOM(t+i,i,unit==units,1),...
             peakPatchCOM(t+i,i,unit==units,2),'*k');
        xlim([-600,600]),ylim([-350,350])                    
        title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);

    end
% $$$     print(gcf,'-depsc2',fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_xtrj_mpos',...
% $$$                                  ['ufr_pfs_',num2str(unit),'.eps']));
% $$$     print(gcf,'-dpng',  fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_xtrj_mpos',...
% $$$                                  ['ufr_pfs_',num2str(unit),'.png']));
    unit = figure_controls(hfig,unit,units,autoincr);
end




%% final figure ori - landscape
sunits = [5,18,107];

spOpts.width  = 2;
spOpts.height = 3;
spOpts.ny = numel(sunits)+2;
spOpts.nx = numel(T);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 8;
figOpts.position = [1,1,(spOpts.height+round(spOpts.padding/2))*spOpts.ny,...
                     (spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2)+figOpts.headerPadding+figOpts.footerPadding];



mkdir(fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_xySdist'))
figHnum = 399329239;;
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',figOpts.position)
set(hfig,'PaperPositionMode','auto');

hfig = figure(399329239);
unit =units(1);
i = 1;
clf

for unit = sunits
    for t = 1:nt,
        pf = pfs{t+1,i};
        uind = find(unit==sunits);
        sp(t,i) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width+round(spOpts.padding/2))*(t-1)+round(spOpts.padding/2),...
                            (spOpts.height +round(spOpts.padding/2))*(uind-1)+round(spOpts.padding/2)+figOpts.footerPadding,...
                            spOpts.width,...
                            spOpts.height]...
                       );
        
        hold('on')
        ratemap = pf.plot(unit,'isCircular',false);
        minv = min(ratemap(:));
        maxv = sq(mRate(1,i,unit));
        ratemap(isnan(ratemap)) = minv-((maxv-minv)/20);
        imagesc(pf.adata.bins{2},pf.adata.bins{1},ratemap);    
        text(pf.adata.bins{2}(1)+30,pf.adata.bins{1}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        ca = caxis;
        caxis([ca(1),sq(mRate(1,i,unit)).*1.5]);
        title([pf.session.trialName ':' pf.parameters.states]);
        plot(peakPatchCOM(t+1,i,unit==units,2),...
             peakPatchCOM(t+1,i,unit==units,1),'*k');
        Lines([],-200.*mod(t+1,2),'k');
        ylim([-600,600]),xlim([-350,350])        
        set(gca,'YTickLabel',{})
        set(gca,'XTickLabel',{})
    end
end



spCohere = pfs{1,1}.spatialCoherence(units);
ind =spCohere>0.99;

for i = 0:2:8;
axes('Units',spOpts.units,...
     'Position',[(i/2)*3+2,3.5,2.5,2.5]);
hold on
plot(sq(pfShift(i/2+2,i/2+3,1,ind))/10,sq(pfShift(i/2+2,i/2+3,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');
if i>0,set(gca,'YTickLabel',{}),end
end


% $$$     print(gcf,'-depsc2',fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_xySdist',...
% $$$                                  ['ufr_pfs_',num2str(unit),'.eps']));
% $$$     print(gcf,'-dpng',  fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_xySdist',...
% $$$                                  ['ufr_pfs_',num2str(unit),'.png']));



%% Portrait

sunits = [1,5,18];

spOpts.width  = 4;
spOpts.height = 2;
spOpts.ny = numel(sunits);
spOpts.nx = numel(T);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 8;
figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                     (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];


figSet = 1;
FigDir = 'Ed10-20140820-shift_teleport_ufr_pfs_xySdist_hvel';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 399329239;;
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',figOpts.position)
set(hfig,'PaperPositionMode','auto');

hfig = figure(399329239);
unit =units(1);
i = 3;
clf

for unit = sunits
    for t = 1:nt,
        pf = pfs{t+1,i};
        uind = find(unit==sunits);
        % Create axes
        sp(t,i) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width +round(spOpts.padding/2))*(uind-1)+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1+2)+round(spOpts.padding/2),...
                            spOpts.width,...
                            spOpts.height]...
                       );
        
        hold('on')

        % Correct color of nans and plot place field
        ratemap = pf.plot(unit,'isCircular',false);
        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
        text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        caxis([-1,sq(mRate(1,i,unit)).*1.5]);
        
        if uind==1,ylabel([pf.session.trialName ':' pf.parameters.states]);end

        % Plot cross over 
        plot(peakPatchCOM(t+1,i,unit==units,1),...
             peakPatchCOM(t+1,i,unit==units,2),'*k');
        Lines(-200.*mod(t+1,2),[],'k');
        xlim([-600,600]),ylim([-350,350])        
        set(gca,'YTickLabel',{})
        set(gca,'XTickLabel',{})
        
        % set colormap 
        
    end
end



spCohere = pfs{1,1}.spatialCoherence(units);
%ind =spCohere>0.99;
ind = sq(all(mRate(:,3,units)>3))&spCohere>0.99;

% Inter trial x shift dist
for i = 0:2:8;
    axes('Units',spOpts.units,...
         'Position',[18,(4-i/2)*3+2.5,2.5,2.5]);
    hold on
    xlim([-40,40]),ylim([-40,40])
    Lines([],0,[.75,.75,.75]);
    Lines(0,[],[.75,.75,.75]);
    Lines(20,[],'r');
    Lines(-20,[],'r');
% $$$ Lines([],20,'r');
% $$$ Lines([],-20,'r');
    plot(-sq(pfShift(i/2+2,i/2+3,3,ind,1))/10,-sq(pfShift(i/2+2,i/2+3,3,ind,2))/10,'.')
    hold on,
    error_ellipse(cov([-sq(pfShift(i/2+2,i/2+3,3,ind,1))/10,-sq(pfShift(i/2+2,i/2+3,3,ind,2))/10]))

    if i~=8,set(gca,'XTickLabel',{}),end
end


% Fist to current x shift dist
for i = 0:2:8;
    axes('Units',spOpts.units,...
         'Position',[22,(4-i/2)*3+2.5,2.5,2.5]);
    hold on
    xlim([-40,40]),ylim([-40,40])
    Lines([],0,[.75,.75,.75]);
    Lines(0,[],[.75,.75,.75]);
    Lines(20,[],'r');
    Lines(-20,[],'r');
% $$$ Lines([],20,'r');
% $$$ Lines([],-20,'r');
    plot(-sq(pfShift(2,i/2+3,3,ind,1))/10,-sq(pfShift(2,i/2+3,3,ind,2))/10,'.')
    error_ellipse(cov([-sq(pfShift(2,i/2+3,3,ind,1))/10,-sq(pfShift(2,i/2+3,3,ind,2))/10]))
    if i~=8,set(gca,'XTickLabel',{}),end
end


print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['ufr_pfs_',num2str(figSet),'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['xShiftD_',num2str(figSet),'.png']));




units = units(ind);

%% Gererate unit rate maps (Place Fields)

MTAstartup('vr_exp');
triallist = 'Ed10VR_teleport';
OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';
T = SessionList(triallist,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');
T(3).offsets = [15,-90];


Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);



units =[1,5,7,9,16,18,22,28,29,99,101,104,107,110,122,134,158,168,184,185]';

Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

    
binDims = [20,20];
numIter = 1001;
nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 2;
sampleRate = 30;
pfk = {};
overwrite = false;
i = 3;

for t = 2:nt-1
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial); 
    xyz = Trial.load('xyz');
    xyz.resample(sampleRate);

    pfk{t-1} = MTAAknnpfs_bs(Trial,units,states{i},overwrite, ...
                             'binDims',binDims,...
                             'nNearestNeighbors',nNearestNeighbors,...
                             'ufrShufBlockSize',ufrShufBlockSize,...
                             'distThreshold',distThreshold,...
                             'pos',xyz,...
                             'numIter',numIter);
end


t = nt;
Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
Trial.stc = Stc.copy;
Trial.stc.load(Trial); 
xyz = Trial.load('xyz');
xyz.resample(sampleRate);


pfk{t-1} = MTAAknnpfs_bs(Trial,units,[states{i},'-shifted'],overwrite, ...
                       'binDims',binDims,...
                       'nNearestNeighbors',nNearestNeighbors,...
                       'ufrShufBlockSize',ufrShufBlockSize,...
                       'distThreshold',distThreshold,...
                       'pos',xyz,...
                       'numIter',numIter);

t = nt+1;
pfk{t-1} = MTAAknnpfs_bs(Trial,units,['shifted&',states{i}],overwrite, ...
                       'binDims',binDims,...
                       'nNearestNeighbors',nNearestNeighbors,...
                       'ufrShufBlockSize',ufrShufBlockSize,...
                       'distThreshold',distThreshold,...                       
                       'pos',xyz,...
                       'numIter',numIter);



%% Place Field Statistics 


pfkstats = {};
pfkshuff = {};
% Test this version should be able to run multiple units at once
for t = 2:nt
    for u = 1:numel(units)        
        Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
        Trial.stc = Stc.copy;
        Trial.stc.load(Trial);
        [pfkstats{t-1,u},pfkshuff{t-1,u}] = PlaceFieldStats(Trial,pfk{t-1},units(u));
    end
end


for u = 1:numel(units)
    [pfkstats{t,u},pfkshuff{t,u}] = PlaceFieldStats(Trial,pfk{t},units(u));
end



for t = 1:nt
    for k = 1:numIter,
        % Retrieve the patch center of mass from patch with the highest firing rate
        pcom = ...
        %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
            cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                    pfkshuff(t,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);

        pind = ~cellfun(@isempty,pcom);
        peakPatchCOM(t,k,pind,:) = ...
            cell2mat(cellfun(@(x) x(:,1),...
                             pcom(pind),...
                             'uniformoutput',false))';

        
        % Retrieve the max patch Firing rate
        peakPatchRate(t,k,:) = ...
            sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                       pfkshuff(t,:),...
                       repmat({k},[1,numel(units)])));

        
        % Retrieve the patch area from patch with highest firing rate
        parea = ...
            cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                    pfkshuff(t,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);
        pind = ~cellfun(@isempty,parea);
        peakPatchArea(t,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                      parea(pind),...
                                                      'uniformoutput',false))');
    end
end


% Units with double field requires special treatment 
% units =[1,5,7,9,16,18,22,28,29,99,101,104,107,110,122,134,158,168,184,185]';
%  
dunits = [22,29,99,104,134];
dims = 1;
cens = 2;
for d = dunits,
    fprintf('double placefield unit: %i',d);
    for t = 1:nt

        pcomx = reshape(pfkshuff{t,d==units}.patchCOM(1,:,:,2),[],1);
        pcomy = reshape(pfkshuff{t,d==units}.patchCOM(1,:,:,1),[],1);
        pcind = nniz(pcomx)&nniz(pcomy);
        pcomx = pcomx(pcind);
        pcomy = pcomy(pcind);

        figure,
        plot(pcomx,pcomy,'.')
        xlim([-600,600,]);        
        ylim([-350,350,]);            
        [cids] = ClusterPP(gcf);    
        delete(gcf)
        
        pcomx = pcomx(cids==1);
        pcomy = pcomy(cids==1);
        
        rsind = randperm(numel(pcomx));
        if numel(rsind)>numIter
            rsind = rsind(1:numIter);
        end
        
        peakPatchCOM(t,:,d==units,2) = zeros;
        peakPatchCOM(t,1:numel(rsind),d==units,2) = pcomx(rsind);
        peakPatchCOM(t,:,d==units,1) = zeros;
        peakPatchCOM(t,1:numel(rsind),d==units,1) = pcomy(rsind);
        
        delete(gcf)
        
    end
end



eds = linspace([-600,600,500]);

display = true;
if display,
    hfig = figure(394929939);
    hold on
    set(0,'defaultAxesFontSize',8,...
          'defaultTextFontSize',8)
    set(hfig,'Units','centimeters')
    set(hfig,'Position',[15,0,5,25]);
    set(hfig,'PaperPositionMode','auto');
    
    FigDir = 'Ed10-20140820-shift_teleport_pfk_hvel_bs_patchSelected';
    mkdir(fullfile(OwnDir,FigDir))
    
    fpind = 1:3:7;
    for u = 1:numel(units);
        clf
        for i = 1:2:5,

            k = floor(i/2)+1;
            j = ceil(i/2)+1;
            
            subplot(9,1,fpind(k));hold on,        
            imagesc(pfk{i}.adata.bins{1},...
                    pfk{i}.adata.bins{2}, ...
                    reshape(pfk{i}.data.rateMap(:,pfk{i}.data.clu==units(u),1), ...
                            fliplr(pfk{i}.adata.binSizes')));
            axis xy
            %hold on,plot(sq(peakPatchCOM(i,:,u,2)),sq(peakPatchCOM(i,:,u,1)),'*m');
            
            xlim([-600,600]);
            ylim([-350,350]);
            if i ==1,
                title(['BS pfk patch pos, unit: ',num2str(units(u))]);
            end
            

            % patchCOM distribution along x axis
            subplot(9,1,i+k); hold on,
            hs = bar(eds,histc(sq(peakPatchCOM(i,nniz(sq(peakPatchCOM(i,:,u,2))'),u,2)),eds),'histc');
            hs.FaceAlpha = 0.5;
            hs.FaceColor = 'c';
            hs.EdgeColor = 'c';
            hs.EdgeAlpha = 0.5;
            
            ho = bar(eds,histc(sq(peakPatchCOM(i+1,nniz(sq(peakPatchCOM(i+1,:,u,2))'),u,2)),eds),'histc');
            ho.FaceAlpha = 0.5;
            ho.FaceColor = 'r';
            ho.EdgeColor = 'r';
            ho.EdgeAlpha = 0.5;
            xlim([-600,600]);
            
            subplot(9,1,i+j); hold on,    
            imagesc(pfk{i+1}.adata.bins{1},...
                    pfk{i+1}.adata.bins{2}, ...
                    reshape(pfk{i+1}.data.rateMap(:,pfk{i+1}.data.clu==units(u),1), ...
                            fliplr(pfk{i+1}.adata.binSizes')));
            axis xy
            %hold on,plot(sq(peakPatchCOM(i+1,:,u,2)),sq(peakPatchCOM(i+1,:,u,1)),'*m');        
            xlim([-600,600]);
            ylim([-350,350]);

        end
        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfs_bs_patch',num2str(units(u)),'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfs_bs_patch',num2str(units(u)),'.png']));
        
    end
end




figure,imagesc(pfk{1}.adata.bins{1},pfk{1}.adata.bins{2}, ...
               reshape(pfk{1}.data.rateMap(:,1,1),fliplr(cellfun(@numel,pfk{1}.adata.bins)))),axis xy
hold on,plot(pfkshuff{1,1}.patchCOM(1,1,1,2),pfkshuff{1,1}.patchCOM(1,1,1,1),'*w')
hold on,plot(sq(peakPatchCOM(1,:,1,2)),sq(peakPatchCOM(1,:,1,1)),'*m');

figure,hold on,
hs = bar(eds,histc(sq(peakPatchCOM(1,:,1,2)),eds),'histc');
hs.FaceAlpha = 0.5;
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
hs.EdgeAlpha = 0.5;
xlim([-600,600,]);



% Calculate the dprime for each distribit
dprime = @(x,y) (mean(x(nniz(x')))-mean(y(nniz(y'))))/sqrt(0.5*(var(x(nniz(x')))+var(y(nniz(y')))));
dprx = nan([nt-1,numel(units)]);
dpry = nan([nt-1,numel(units)]);

for i = 1:5,
    for u = 1:numel(units);    
        dprx(i,u) = dprime(peakPatchCOM(i,:,u,2),peakPatchCOM(i+1,:,u,2));
        dpry(i,u) = dprime(peakPatchCOM(i,:,u,1),peakPatchCOM(i+1,:,u,1));
    end    
end



FigDir = 'Ed10-20140820-shift_teleport_dprime_xyshift_patchSelected';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 84827378;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[2,0,14,6])
set(hfig,'PaperPositionMode','auto');
for i = 1:3,
    axes('Units','centimeters',...
         'Position',[2+(i-1)*(2.5+2),2,2.5,2.5]);
    plot(dprx(i,:),dpry(i,:),'.')
    xlim([-20,20]);
    ylim([-20,20]);   
    %Lines(nanmean(dprx(i,:)),[],'k');
    %Lines([],nanmean(dpry(i,:)),'k');
    xlabel('x-shift (d-prime)');
    ylabel('y-shift (d-prime)');
    grid on
    title({T(1).sessionName,[T(i+1).trialName,' vs ',T(i+2).trialName]});    
end
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_dprime.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_dprime.png']));




%% Plot dprx vs xshift

FigDir = 'Ed10-20140820-shift_teleport_dprx_xshift_patchSelected';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 84827377;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[2,0,17,6])
set(hfig,'PaperPositionMode','auto');
for i = 1:3,
    axes('Units','centimeters',...
         'Position',[2+(i-1)*(2.5+2),2,2.5,2.5]);
    plot(-sq(peakPatchCOM(i,1,:,2)-peakPatchCOM(i+1,1,:,2))/10,dprx(i,:),'.');
    xlim([-40,40]);
    ylim([-10,10]);
    Lines([],nanmean(dprx(i,:)),'r');
    Lines(nanmean(-sq(peakPatchCOM(i,1,:,2)-peakPatchCOM(i+1,1,:,2))/10),[],'r');
    Lines([],0,'k');
    Lines(0,[],'k');    
    xlabel('pfk x Shift (cm)');
    ylabel('dprime between conditions');
    grid on
    title({T(1).sessionName,[T(i+1).trialName,' vs ',T(i+2).trialName]});    
end
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_dprx_x.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_dprx_x.png']));




for i = 1:4,
    for u = 1:numel(units),
        mpcom(i,u,2) = mean(peakPatchCOM(i,nniz(sq(peakPatchCOM(i,:,u,2))'),u,2));
        mpcom(i,u,1) = mean(peakPatchCOM(i,nniz(sq(peakPatchCOM(i,:,u,1))'),u,1));    
    end
end

    
FigDir = 'Ed10-20140820-shift_teleport_dprx_mean_xshift_patchSelected';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 84827377;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[2,0,17,6])
set(hfig,'PaperPositionMode','auto');
for i = 1:3,
    axes('Units','centimeters',...
         'Position',[2+(i-1)*(2.5+2),2,2.5,2.5]);
    plot(-sq(mpcom(i,:,2)-mpcom(i+1,:,2))/10,dprx(i,:),'.');
    xlim([-40,40]);
    ylim([-10,10]);
    Lines([],nanmean(dprx(i,:)),'r');
    Lines(nanmean(-sq(mpcom(i,:,2)-mpcom(i+1,:,2))/10),[],'r');
    Lines([],0,'k');
    Lines(0,[],'k');    
    xlabel('pfk x Shift (cm)');
    ylabel('dprime between conditions');
    grid on
    title({T(1).sessionName,[T(i+1).trialName,' vs ',T(i+2).trialName]});    
end
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_dprx_x.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_dprx_x.png']));








%% Portrait with pfknn



sunits = [1,5,18];
t = 1;
mRate = [];
for t = 1:nt;
    mRate(t,:) = pfk{t}.maxRate(sunits);
end

spOpts.width  = 4;
spOpts.height = 2;
spOpts.ny = numel(sunits);
spOpts.nx = numel(T);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 8;
figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                     (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];


figSet = 1;
FigDir = 'Ed10-20140820-shift_teleport_ufr_pfk_xySdist_hvel';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 399329239;;
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',figOpts.position)
set(hfig,'PaperPositionMode','auto');

hfig = figure(399329239);
unit =units(1);
i = 3;
clf

for unit = sunits
    for t = 1:4,
        pf = pfk{t};
        uind = find(unit==sunits);
        % Create axes
        sp(t,i) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width +round(spOpts.padding/2))*(uind-1)+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1+2)+round(spOpts.padding/2),...
                            spOpts.width,...
                            spOpts.height]...
                       );
        
        hold('on')

        % Correct color of nans and plot place field
        ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,1),fliplr(pf.adata.binSizes'));
        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
        text(pf.adata.bins{1}(end)-200,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        %caxis([-1,sq(mRate(1,unit==sunits)).*1.5]);
        caxis([-1,8]);        
        if uind==1,ylabel([pf.session.trialName ':' pf.parameters.states]);end

        % Plot cross over 
        plot(mpcom(t,unit==units,2),...
             mpcom(t,unit==units,1),'*k');
        Lines(-200.*mod(t+1,2),[],'k');
        xlim([-600,600]),ylim([-350,350])        
        set(gca,'YTickLabel',{})
        set(gca,'XTickLabel',{})
        
        % set colormap 
        
    end
end




% Inter trial x shift dist
for i = 0:2:4;
    axes('Units',spOpts.units,...
         'Position',[18,(4-i/2)*3+2.5,2.5,2.5]);
    hold on
    xlim([-40,40]),ylim([-40,40])
    Lines([],0,[.75,.75,.75]);
    Lines(0,[],[.75,.75,.75]);
    Lines(20,[],'r');
    Lines(-20,[],'r');
    plot(-(mpcom(i/2+1,:,2)-mpcom(i/2+2,:,2))/10,-(mpcom(i/2+1,:,1)-mpcom(i/2+2,:,1))/10,'.')

    if i~=4,set(gca,'XTickLabel',{}),end
end



print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_xyshift',num2str(figSet),'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_xyshift',num2str(figSet),'.png']));




% Intra-trial variance

figure,
subplot(211);bar(-200:5:200,hist(sq(peakPatchCOM(1,2:2:1001,2,2)-peakPatchCOM(1,3:2:1001,2,2)),-200:5:200),'histc')
subplot(212);bar(-200:5:200,hist(sq(peakPatchCOM(2,2:2:1001,2,2)-peakPatchCOM(2,3:2:1001,2,2)),-200:5:200),'histc')


figure,hist(sq(peakPatchCOM(1,2:2:1001,2,1)-peakPatchCOM(2,3:2:1001,2,1)),100)

figure,hist(sq(peakPatchCOM(1,2:2:1001,4,1)),100)



unit = 1;
t = 1;

X = sq(peakPatchCOM(t  ,2:2:1001, unit,2));
Y = sq(peakPatchCOM(t+1,3:2:1001, unit,2));
nxind = nniz(X(:));
nyind = nniz(Y(:));

figure,
subplot(411);
pf = pfk{t};
ratemap = reshape(pf.data.rateMap(:,unit,1),fliplr(pf.adata.binSizes'));
imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
axis xy
subplot(412);hist(X(nxind),100),xlim([-600,600])
subplot(413);hist(Y(nyind),100),xlim([-600,600])
subplot(414);
pf = pfk{t+1};
ratemap = reshape(pf.data.rateMap(:,unit,1),fliplr(pf.adata.binSizes'));
imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
axis xy
 
[H,P,CI] = ttest2(X(nxind),Y(nyind),'tail','both','vartype','unequal')
[P,H,STATS] = ranksum(X(nxind),Y(nyind),'tail','both')
[H,P,CI] = ttest2(X(nxind),mean(Y(nyind)),'tail','both')

[H,P,CI] = vartest2(X(nxind),Y(nyind),'tail','both')
[H,P,CI] = ztest(X(nxind),nanmean(Y(nyind)),nanstd(Y(nyind)),'tail','both')



%% Compute permutations

Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);

offsets =  [15,-15;...
            15,-90;...
            15,-15;...
            15,-15];

% add trials as states
trialLabels = {'Coreg1PEVR1','Shift1PEVR2','Coreg2PEVR3','Shift2PEVR4'};

for t = 1:4,
aper = resample(Trial.sync.copy,txyz.sampleRate);
aper.data = aper.data-aper.data(1)+1;
aper = aper+offsets(t,:);
Stc.addState(Trial.spath,...
                   Trial.filebase,...
                   aper.data,...
                   txyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'gper','');
end

states = cellfun(@horzcat,trialLabels,repmat({'&velHthresh'},size(trialLabels)));
states = {'Coreg1PEVR1&velHthresh',...
          'Shift1PEVR2&velHthresh',...
          'Coreg2PEVR3&velHthresh',...
          'Shift2PEVR4&velHthresh'};

units =[1,5,7,9,16,18,22,28,29,99,101,104,107,110,122,134,158,168,184,185]';

Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

    
binDims = [20,20];
numIter = 1001;
nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 2;
sampleRate = 30;
pfk = {};
overwrite = false;

xyz = Trial.load('xyz');
xyz.resample(sampleRate);

for s1 = 1:numel(states)
    for s2 = s1+1:numel(states)

        pfk{s1,s2} = MTAAknnpfs_perm(Trial,units,states([s1,s2]),overwrite, ...
                                     'binDims',binDims,...
                                     'nNearestNeighbors',nNearestNeighbors,...
                                     'ufrShufBlockSize',ufrShufBlockSize,...
                                     'distThreshold',distThreshold,...
                                     'pos',xyz,...
                                     'numIter',numIter);
    end
end

        



