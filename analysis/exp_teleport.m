MTAstartup('vr_exp');
overwriteSession = false;
overwriteTrials  = false;
overwriteStc     = false;
trialList = 'Ed10VR_teleport';

T = SessionList(trialList);

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
    QuickTrialSetup(T); 
end


Trial = MTATrial.validate(T(1)).;
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

    
    Stc.save(1);
end

% Calculate and plot
Stc = Trial.stc.copy;
nt = numel(T);
states = {'theta','velthresh','velHthresh'};
nsts = size(states,2);

display = true;
overwrite = false;
units = 1:185;

% Generate unit auto correlogram
[accg,tbin] = autoccg(Trial,units,'theta');


% Gererate unit rate maps 
%binDims = [60,60];
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



for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,['shifted&',states{i}],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'type','xy');
end

t = t+1;
for i = 1:nsts,
    pfs{t,i} = MTAApfs(Trial,units,[states{i},'-shifted'],overwrite, ...
                       'binDims',binDims,...
                       'SmoothingWeights',smoothingWeights,...
                       'type','xy');
end



t = 1;
mRate = [];
for i = 1:nsts,
    if t==1, 
        mRate(t,i,:) = pfs{t,i}.maxRate;
    end
    
end

units = find(sq(mRate(1,1,:))>3);

if display,    
    spOpts.width  = 4;
    spOpts.height = 2;
    spOpts.ny = numel(T)+1+1;
    spOpts.nx = numel(states);
    spOpts.padding = 2;
    spOpts.units = 'centimeters';

    figOpts.units = 'centimeters';
    figOpts.headerPadding = 4;
    figOpts.footerPadding = 4;
    figOpts.position = [1,1,(spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding];

    
    sp = [];
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
        for t = 1:nt+1,
            for i = 1:nsts,
                sp(t,i) = axes('Units',spOpts.units,...
                               'Position',[(spOpts.width +round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                           (spOpts.height+round(spOpts.padding/2))*(spOpts.ny-t+1)+round(spOpts.padding/2),...
                                            spOpts.width,...
                                            spOpts.height]...
                );

                pfs{t,i}.plot(unit,[],true,[0,max(mRate(1,i,unit)).*1.5],false);
                title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
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

        print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'figures/vr_exp/Ed10-20140820-shift_teleport',...
                                     ['pfs_',num2str(unit),'.eps']));
        print(gcf,'-dpng',fullfile(getenv('PROJECT'),'figures/vr_exp/Ed10-20140820-shift_teleport',...
                                     ['pfs_',num2str(unit),'.png']));
% $$$ % $$$         
% $$$         reportfig('/storage/gravio/figures/',...
% $$$                   hfig,...
% $$$                   ['exp_teleport-' Trial.name],...
% $$$                   'vr_exp',...
% $$$                   false,...
% $$$                   ['Unit: ' num2str(unit)],...  Tag
% $$$                   '',...                Comment
% $$$                   100,...                Resolution
% $$$                   false,...             SaveFig
% $$$                   'png',...             Format
% $$$                   8,...                 Width
% $$$                   12,...                 Height
% $$$                   unit...               Id
% $$$         );
        unit = figure_controls(hfig,unit,units,autoincr);
    
    end

end



unts = [4,25,27,234,235];
hfig = figure(666989);
unit = units(1);

i = 4;

for t=1:5,
for unit = unts,
    sp = subplot2(nt,numel(unts),t,find(ismember(unts,unit)));cla
    pfs{t,i}.plot(unit,[],true);
    title([pfs{t,i}.session.trialName ':' pfs{t,i}.parameters.states,': ',num2str(unit)]);
end
end


reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
          ['exp_teleport_selected_units-' Trial.name] ,false,Trial.name);
    



pfx = {};
for t = 1:nt
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial);
    for i = 1:nsts,
        pfx{t,i} = MTAApfs(Trial,units,states{i},overwrite, ...
                           'binDims',[60],...
                           'SmoothingWeights',[0.8],...
                           'type','x');
    end
end

for i = 1:nsts,
    pfx{t,i} = MTAApfs(Trial,units,[states{i},'&shifted'],overwrite, ...
                       'binDims',[60],...
                       'SmoothingWeights',[0.8],...
                       'type','x');
end

t = t+1
for i = 1:nsts,
    pfx{t,i} = MTAApfs(Trial,units,[states{i},'-shifted'],overwrite, ...
                       'binDims',[60],...
                       'SmoothingWeights',[0.8],...
                       'type','x');
end





t = 1;
i = 1;
u = 3;
unit = 5;
pfstats = {};
pfshuff = {};
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

% Test this version should be able to run multiple units at once
for t = 1:nt
    Trial = MTATrial(T(t).sessionName,T(t).mazeName,T(t).trialName);    
    Trial.stc = Stc.copy;
    Trial.stc.load(Trial);
    for i = 1:nsts,        
        [pfstats{t,i},pfshuff{t,i}] = PlaceFieldStats(Trial,pfs{t,i},units);
    end
end



%pfstats{t,i}.patchCom

pfstats{2,i}-pfstats{3,i}
pfstats{4,i}-pfstats{5,i}



% 1D version of place fields
unit = units(1);
autoincr = false;
figHnum = 666995;
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',figOpts.position)
set(hfig,'PaperPositionMode','auto');
while unit~=-1,
    clf;
    hold('on');
    i = 3;
    for t = 2:nt-1,
        plot(pfx{t,i}.adata.bins{1},pfx{t,i}.data.rateMap(:,pfx{t,i}.data.clu==unit));
    end
    legend({T(2:end-1).trialName});
    unit = figure_controls(hfig,unit,units,autoincr);
end


% time X xyz X ufr
xyz = Trial.load('xyz');
xyz.resample(10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'theta',units,0.2);
unit = units(3);
c = zeros([xyz.size(1),3]);
c(ufr(:,unit)>=5,1) = 1;
figure, clf
pscatter(1:xyz.size(1),xyz(:,6,1),10,c);
Lines(round((Trial.sync.data(:)-Trial.sync.data(1)).*10)+1,[],'m');


