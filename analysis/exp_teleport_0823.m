MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;
trialList = 'Ed10VR_20160823';
OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';

T = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');


if overwriteSession,
    Session = MTASession(T(1).sessionName,  ...
                         T(1).mazeName,     ...
                         true,              ...
                         T(1).TTLValue,     ...
                         'vicon',           ...
                         'nlx',             ...
                         T(1).xyzSampleRate ...
    );
    %Update spk object
    Session.spk.create(Session);

    pXY(Session)    
    xyz = Session.load('xyz');
    xyz.data(:,:,1) = xyz.data(:,:,1)+T(1).xOffSet;
    xyz.data(:,:,2) = xyz.data(:,:,2)+T(1).yOffSet;
    xyz.save;
    pXY(Session)        
    Session.save;
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

% $$$     if isempty(Trial.stc.gsi('r')),
% $$$         rper = rear(Trial,'com',45);
% $$$         Trial.stc.addState(Trial.spath,...
% $$$                            Trial.filebase,...
% $$$                            rper,...
% $$$                            xyz.sampleRate,...
% $$$                            Trial.sync.copy,...
% $$$                            Trial.sync.data(1),...
% $$$                            'rear','r');
% $$$     end
% $$$ 
% $$$ 
% $$$     if isempty(Trial.stc.gsi('n')),
% $$$         Trial.stc.states{end+1} = Trial.stc{'v'}-(Trial.stc{'r',120}+[-.5,.5]);
% $$$         Trial.stc.states{end}.key = 'n'; 
% $$$         Trial.stc.states{end}.label = 'NRvel';    
% $$$         Trial.stc.states{end}.updateFilename([Trial.filebase,'.sst.',...
% $$$                             Trial.stc.states{end}.label,'.',...
% $$$                             Trial.stc.states{end}.key,'.mat']);
% $$$     end

    if isempty(Trial.stc.gsi('t')),
        Trial = labelTheta(Trial,[],32);
    end

    
    Trial.stc.save(1);
end

% Calculate and plot
Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

overwrite = false;
%units = [];
units = [13,17,18,20,21,28,32,42,48,56,57,58,59,62,64,68,77,80,87,107,112,116,120,142,147,160,164,167,172,174];

% Generate unit auto correlogram
[accg,tbin] = autoccg(Trial,units,'theta');




%% Gererate unit rate maps (Place Fields)
binDims = [40,40];
smoothingWeights = [1.2,1.2];
pfs = {};
for t = 1:nt
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial); 
   for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite, ...
                           'binDims',binDims,'SmoothingWeights',smoothingWeights);
    end
end



%% Select Units with Firing rate greater than 3Hz in all.theta
t = 1;
mRate = [];
for i = 1:nsts,
    if t==1, 
        mRate(t,i,:) = pfs{t,i}.maxRate;
    end
    
end
%units = find(sq(mRate(1,1,:))>3);




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
for t = 1:nt,
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
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 2;
figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                     (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];


mkdir(fullfile(OwnDir,[Trial.name '-shift_teleport_pfs']));
sp = [];
%autoincr = false;
autoincr = true;

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
    for t = 1:nt,
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
            ratemap = ratemap;            
            ratemap(isnan(ratemap)) = -1;
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    

            text(pf.adata.bins{1}(1)+950,pf.adata.bins{2}(end)-70,...
                 sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
            plot(peakPatchCOM(t,i,unit==units,1),...
                 peakPatchCOM(t,i,unit==units,2),'*k');
            xlim([-600,600]),ylim([-350,350])                    
            title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);
        end
    end
    ForAllSubplots('colormap([0,0,0;parula])');
    ForAllSubplots(['caxis([-1,',num2str(sq(mRate(1,i,unit))),'.*1.5])']);
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

    print(gcf,'-depsc2',fullfile(OwnDir,[Trial.name '-shift_teleport_pfs'],...
                                 ['pfs_',num2str(unit),'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,[Trial.name '-shift_teleport_pfs'],...
                                 ['pfs_',num2str(unit),'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end



%% Old stuff starts here



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



t = 2;
i = 3;
spCohere = pfs{t,i}.spatialCoherence(units);

figure,hold on,
ind =spCohere>0.99;
subplot(121)
plot(sq(pfShift(2,3,1,ind))/10,sq(pfShift(2,3,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');
subplot(122)
plot(sq(pfShift(3,4,1,ind))/10,sq(pfShift(3,4,2,ind))/10,'.')
xlim([-40,40]),ylim([-40,40])
Lines([],0,'k');
Lines([],20,'r');
Lines([],-20,'m');


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






% final fig 


sunits = [13,18,42];

spOpts.width  = 4;
spOpts.height = 2;
spOpts.ny = numel(sunits)+2;
spOpts.nx = numel(T);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 8;
figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                     (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];


FigDir = 'Ed10-20140823-shift_teleport_ufr_pfs_xySdist';
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
    for t = 1:nt-1,
        pf = pfs{t+1,i};
        uind = find(unit==sunits);
        sp(t,i) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width +round(spOpts.padding/2))*(uind-1)+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                            spOpts.width,...
                            spOpts.height]...
                       );
        
        hold('on')
        ratemap = pf.plot(unit,'isCircular',false);
        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
        text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        caxis([-1,sq(mRate(1,i,unit)).*1.5]);
        
        if uind==1,ylabel([pf.session.trialName ':' pf.parameters.states]);end
        plot(peakPatchCOM(t+1,i,unit==units,1),...
             peakPatchCOM(t+1,i,unit==units,2),'*k');
        Lines(-200.*mod(t+1,2),[],'k');
        xlim([-600,600]),ylim([-350,350])        
        set(gca,'YTickLabel',{})
        set(gca,'XTickLabel',{})
    end
end



spCohere = pfs{1,1}.spatialCoherence(units);
%ind =spCohere>0.99;
ind = sq(mRate(1,3,units)>3)&spCohere>0.99;

for i = 0:2:2;
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
plot(-sq(pfShift(i/2+2,i/2+3,2,ind))/10,-sq(pfShift(i/2+2,i/2+3,1,ind))/10,'.')

if i~=8,set(gca,'XTickLabel',{}),end
end



print(gcf,'-depsc2',fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_xySdist',...
                             ['ufr_pfs_',num2str(unit),'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,'Ed10-20140820-shift_teleport_ufr_pfs_xySdist',...
                             ['ufr_pfs_',num2str(unit),'.png']));


figure,plot(sq(peakPatchCOM(2,3,ind,1)),-sq(pfShift(2,3,3,ind,1))/10,'.')
 grid on
 
hold on,plot(sq(peakPatchCOM(3,3,ind,1)),-sq(pfShift(3,4,3,ind,1))/10,'.')

figure,plot(sq(peakPatchCOM(3,3,ind,1)),-sq(pfShift(3,4,3,ind,1))/10,'.')
grid on









%% xy decode

addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/fwdrebayesiandecodingtools/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/tempScripts/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/aux/

dunits = units(spCohere>0.99);


ratemap=[];
for u = 1:numel(dunits),
    ratemap(:,:,u) = pfs{2,1}.plot(dunits(u),'isCircular',false);
end


Trial = MTATrial.validate(T(1));
Trial.stc = Stc.copy;
Trial.stc.load(Trial); 
xyz = Trial.load('xyz');
xyz.resample(10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'theta',dunits,1.5,true);

E = decode_bayesian_poisson(ratemap,ufr.data');

xbins = pfs{2,1}.adata.bins{1};
ybins = pfs{2,1}.adata.bins{2};

pos = nan([size(E,3),2]);
for tind = 1:size(E,3),
    try,
        xyi = LocalMinimaN(-E(:,:,tind),0,100);        
        pos(tind,:) = [xbins(xyi(1)),ybins(xyi(2))];
    end
end    

tper = Trial.stc{'t',xyz.sampleRate};
tper.cast('TimeSeries');
tper.data = logical(tper.data);

txyz = xyz.copy;
txyz.data(~tper.data,:,:)=nan;

mpos = xyz.copy;
mpos.data = pos;
mpos.data(~tper.data,:)=nan;

figure,plot(txyz(:,5,1))
hold on,
plot(mpos(:,1))
Lines(round((Trial.sync(:)-Trial.sync(1))*10)+1,[],'r')




figure,plot(r(Trial.stc{'t'}))


xp = sq(xyz(:,5,[1,2]))-pos;
r = xyz.copy;
[th,r.data] = cart2pol(xp(:,1),xp(:,2));



figure,hist(r(Trial.stc{'t'}),100)



figure,
plot(r)
hold on
plot(sum(ufr.data,2))


figure,plot(sq(xyz(:,5,[1])))
hold on,
[~,ttt] = min(-mean(E,2));
%ttt(ttt==0)=1;
plot(xbins(sq(ttt)));



xd = sq(xyz(:,5,[1]))-pos(:,1);

figure,
plot(xd.*double(sum(ufr.data,2)>10))
hold on
plot(sum(ufr.data,2))


tind = 800:1300;
figure,

plot(xyz(tind,5,1),xyz(tind,5,2),'.')
hold on
plot(pos(tind,1).*sum(ufr.data(tind,:),2)~=0),pos(tind,2).*sum(ufr.data(tind,:),2)==0),'.')


Trial = MTATrial.validate(T(3));
pxyz = Trial.load('xyz');
figure,plot(pxyz(:,5,1))


%% PFK bootstrap patch COM

MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;
trialList = 'Ed10VR_20160823';
OwnDir = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/';

T = SessionList(trialList,...
                '/storage/gravio/data/processed/xyz/Ed10/',...
                '/storage/eduardo/data/processed/nlx/Ed10/');

Trial = MTATrial.validate(T(1));
Trial.load('stc',T(1).stcMode);

%units =[1,5,7,9,16,18,22,28,29,94,99,101,104,107,110,122,134,158,168,184,185]';
units = [13,17,18,20,21,28,32,42,48,56,57,58,59,62,64,68,77,80,87,107,112,116,120,142,147,160,164,167,172,174];

Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

    
binDims = [20,20];
numIter = 1000;
nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 1;
sampleRate = 30;
pfk = {};
overwrite = false;
i = 3;

for t = 2:nt
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


for t = 1:nt-1
    for k = 1:numIter,
        % Retrieve the patch center of mass from patch with the highest firing rate
        pcom = ...
            cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
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





% $$$ figure,imagesc(pfk{1}.adata.bins{1},pfk{1}.adata.bins{2}, ...
% $$$                reshape(pfk{1}.data.rateMap(:,1,1),fliplr(cellfun(@numel,pfk{1}.adata.bins)))),axis xy
% $$$ hold on,plot(pfkshuff{1,1}.patchCOM(1,1,1,2),pfkshuff{1,1}.patchCOM(1,1,1,1),'*w')
% $$$ hold on,plot(sq(peakPatchCOM(1,:,1,2)),sq(peakPatchCOM(1,:,1,1)),'*m');
% $$$ 
% $$$ figure,hold on,
% $$$ hs = bar(eds,histc(sq(peakPatchCOM(1,:,1,2)),eds),'histc');
% $$$ hs.FaceAlpha = 0.5;
% $$$ hs.FaceColor = 'c';
% $$$ hs.EdgeColor = 'c';
% $$$ hs.EdgeAlpha = 0.5;
% $$$ xlim([-600,600,]);




eds = linspace([-600,600,500]);

display = true;
if display,
    hfig = figure(394929939);
    hold on
    set(0,'defaultAxesFontSize',8,...
          'defaultTextFontSize',8)
    set(hfig,'Units','centimeters')
    set(hfig,'Position',[15,0,20,16]);
    set(hfig,'PaperPositionMode','auto');
    
    FigDir = 'Ed10-20140823-shift_teleport_pfk_pfs_hvel_bs_patch';
    mkdir(fullfile(OwnDir,FigDir))
    
    k = 1:2:6;
    for u = 1:numel(units);
        clf
        for i = 1:3

            % Setup subplot
            sp(i) = subplot2(6,3,[k(i),k(i)+1],1);
            hold('on');

            
            pf = pfk{i};
            % Plot ratemap
            ratemap = reshape(pf.data.rateMap(:,pf.data.clu==units(u),1), ...
                            fliplr(pf.adata.binSizes'));
            ratemap(isnan(ratemap)) = -1;
            imagesc(pf.adata.bins{1},...
                    pf.adata.bins{2}, ...
                    ratemap);
            colormap([0,0,0;parula]);
            caxis([-1,sq(mRate(1,3,u)).*1.5]);
            axis xy

            % Plot placefield center
            plot(sq(nanmean(peakPatchCOM(i,nniz(peakPatchCOM(i,:,u,2)')&nniz(peakPatchCOM(i,:,u,1)'),u,2))),...
                 sq(nanmean(peakPatchCOM(i,nniz(peakPatchCOM(i,:,u,1)')&nniz(peakPatchCOM(i,:,u,2)'),u,1))),'*m');
            text(pf.adata.bins{1}(end)-250,pf.adata.bins{2}(end)-50,...
                 sprintf('%2.1f',max(pf.data.rateMap(:,pf.data.clu==units(u),1))),'Color','w','FontWeight','bold','FontSize',10)

            
            xlim([-600,600]);
            ylim([-350,350]);
            if i ==1,
                title(['BS pfk patch pos, unit: ',num2str(units(u))]);
            end
            
        end

        for i = 1:3

            % Plot place fields
            sp(i+3) = subplot2(6,3,[k(i),k(i)+1],2);
            pf = pfs{i+1,3};
            ratemap = pf.plot(units(u),'isCircular',false);
            ratemap(isnan(ratemap)) = -1;
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    

            text(pf.adata.bins{1}(1)+950,pf.adata.bins{2}(end)-70,...
                 sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
            title([pf.session.trialName ':' pf.parameters.states,': ',num2str(units(u))]);
            colormap([0,0,0;parula]);
            caxis([-1,sq(mRate(1,3,u)).*1.5]);
            
            axis xy
            hold on,
            plot(sq(nanmean(peakPatchCOM(i,nniz(peakPatchCOM(i,:,u,2)')&nniz(peakPatchCOM(i,:,u,1)'),u,2))),...
                 sq(nanmean(peakPatchCOM(i,nniz(peakPatchCOM(i,:,u,1)')&nniz(peakPatchCOM(i,:,u,2)'),u,1))),'*m');

            xlim([-600,600]);
            ylim([-350,350]);
            if i ==1,
                title(['BS pfs unit: ',num2str(units(u))]);
            end
            
        end
        

        
        for i = 1:2
            % patchCOM distribution along x axis
            sp(i+6) = subplot2(6,3,[k(i),k(i)+1]+1,3);
            hold on,
            hs = bar(eds,histc(sq(peakPatchCOM(i,:,u,2)),eds),'histc');
            hs.FaceAlpha = 0.5;
            hs.FaceColor = 'c';
            hs.EdgeColor = 'c';
            hs.EdgeAlpha = 0.5;
            
            ho = bar(eds,histc(sq(peakPatchCOM(i+1,:,u,2)),eds),'histc');
            ho.FaceAlpha = 0.5;
            ho.FaceColor = 'r';
            ho.EdgeColor = 'r';
            ho.EdgeAlpha = 0.5;
            xlim([-600,600]);
        end        

        for i = 1:8,
            sp(i).Units = 'centimeters';
            sp(i).Position([3,4]) = [4,2];
        end        
        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfs_bs_patch',num2str(units(u)),'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfs_bs_patch',num2str(units(u)),'.png']));
        %waitforbuttonpress        
    end
end

dprime = @(x,y) (mean(x(nniz(x')))-mean(y(nniz(y'))))/sqrt(0.5*(var(x(nniz(x')))+var(y(nniz(y')))));
dprx = nan([nt-2,numel(units)]);
dpry = nan([nt-2,numel(units)]);

for i = 1:2,
    for u = 1:numel(units);    
        dprx(i,u) = dprime(peakPatchCOM(i,:,u,2),peakPatchCOM(i+1,:,u,2));
        dpry(i,u) = dprime(peakPatchCOM(i,:,u,1),peakPatchCOM(i+1,:,u,1));
    end    
end




FigDir = 'Ed10-20140823-shift_teleport_dprime_xyshift';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 84827377;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[2,0,16,6])
set(hfig,'PaperPositionMode','auto');
for i = 1:2,
    axes('Units','centimeters',...
         'Position',[2+(i-1)*(2.5+2),2,2.5,2.5]);
    plot(dprx(i,:),dpry(i,:),'.');
    xlabel('x dprime')
    ylabel('y dprime')    
    xlim([-10,10]),ylim([-10,10]) ;   
    Lines(nanmean(dprx(i,:)),[],'k');
    Lines([],nanmean(dpry(i,:)),'k');
    grid on
    title({[T(1).sessionName                        ],...
           [T(i+1).trialName,' vs ',T(i+2).trialName],...
           ['Placefield bootstrapped center dprime' ]});    
end
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_dprime.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_dprime.png']));




%% Plot dprx vs xshift

FigDir = 'Ed10-20140823-shift_teleport_dprx_xshift';
mkdir(fullfile(OwnDir,FigDir))
figHnum = 84827377;
hfig = figure(figHnum);clf
set(hfig,'units','centimeters')
set(hfig,'Position',[2,0,18,6])
set(hfig,'PaperPositionMode','auto');
for i = 1:2,
    axes('Units','centimeters',...
         'Position',[2+(i-1)*(2.5+2),2,2.5,2.5]);
    plot(-sq(peakPatchCOM(i,1,:,2)-peakPatchCOM(i+1,1,:,2))/10,dprx(i,:),'.')
    xlim([-10,10]),ylim([-10,10])    
    %Lines(nanmean(-sq(peakPatchCOM(i,1,:,2)-peakPatchCOM(i+1,1,:,2))/10),[],'k')
    %Lines([],nanmean(dprx(i,:)),'k')
    xlabel('pfk Shift')
    grid on
    title({T(1).sessionName,[T(i+1).trialName,' vs ',T(i+2).trialName]});    
end
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['pfk_dprx_x.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,['pfk_dprx_x.png']));







%% Portrait with pfknn


for i = 1:3,
    for u = 1:numel(units),
        mpcom(i,u,2) = mean(peakPatchCOM(i,nniz(sq(peakPatchCOM(i,:,u,2))'),u,2));
        mpcom(i,u,1) = mean(peakPatchCOM(i,nniz(sq(peakPatchCOM(i,:,u,1))'),u,1));    
    end
end


sunits = [20,28,160];
t = 1;
mRate = [];
for t = 1:nt-1;
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
FigDir = 'Ed10-20140823-shift_teleport_ufr_pfk_xySdist_hvel';
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
    for t = 1:3,
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
        text(pf.adata.bins{1}(end)-250,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        %caxis([-1,sq(mRate(1,unit==sunits)).*1.5]);
        caxis([-1,15]);        
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
for i = 0:2:2;
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




