


Slist = get_session_list('ER06');

i = 4;
QuickSessionSetup(Slist(i)); 

s = MTASession.validate(Slist(i));
pXY(s);
pZ(s);
NeuronQuality(s);
PlotSessionErrors(s);

% $$$ xyz.data(:,:,1) = xyz.data(:,:,1)+Slist(2).xOffSet;
% $$$ xyz.data(:,:,2) = xyz.data(:,:,2)+Slist(2).yOffSet;
% $$$ xyz.save;

Tlist = get_session_list('exp_flight');
QuickTrialSetup(s,'includeSyncInd',[Tlist(:).includeSyncInd],'overwrite',true); % Write the session trial 
QuickTrialSetup(Tlist,'overwrite',true);

% Load trial
Trial = MTATrial.validate('ER06-20130614.cof.all');
% Label trial behavior with neural network model
%Trial = labelBhv_NN(Trial,'NN0317');
% $$$ if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end
Stc = Trial.load('stc','default');



%% Place field vars
% SELECT pyramidal units
units = select_units(Trial,18,'pyr');
overwrite = false;
% Initialize place field variable
pfs = {};


%% Pre flight exploration 
% Load pre flight exploration trial
Trial = MTATrial.validate('ER06-20130614.cof.gnd');
Trial.maze.boundaries(end) = 450;

% Label periods of hippocampal theta activity 

Trial.stc = Stc;
Trial.load('stc');
% $$$ Trial.load('stc','NN0317_PP');
%Trial.load('stc','default');
% $$$ if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end
%Trial = labelTheta(Trial,[],73)







pfs{1} = MTAApfs(Trial,units,'theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);


%% Flight p
% Load Flight Trial
Trial = MTATrial.validate('ER06-20130614.cof.fly');
Trial.maze.boundaries(end) = 450;
%if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end
%Trial = labelFlight(Trial,[],1,'set');
Trial.stc = Stc.copy;
Trial.load('stc');



if isempty(Trial.stc.gsi('f')),Trial = labelFlight(Trial);end
fstc = Trial.load('stc','flight');

pZ(Trial)
xyz = Trial.load('xyz');
hold on,plot(xyz(:,1,3));
Lines(fstc{'f'}(:),[],'r');



pfs{2} = MTAApfs(Trial,units,'theta-flight',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);

pfs{3} = MTAApfs(Trial,units,'flight&theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);





%% setup figure paths
OwnDir = '/storage/gravio/nextcloud/';
FigDir = ['MjgER2016/figures/figure3/flight_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));




%% Plot placefields across height
slices = 1:2:min(cellfun(@(x) x.adata.binSizes(end),pfs));

spOpts.width  = 2;
spOpts.height = 2;
spOpts.ny = 3;
spOpts.nx = numel(slices);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 0;
figOpts.position = [1,1,(spOpts.height+round(spOpts.padding/2))*spOpts.ny+round(spOpts.padding/2),...
                     (spOpts.width+round(spOpts.padding/2)) *spOpts.nx+figOpts.headerPadding+figOpts.footerPadding];

width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;
    
    
mRates = [];
for p = 1:3
    mRates(:,p) = pfs{p}.maxRate;
end
mRates = max(mRates,[],2);



hfig = figure(393929);
autoincr = true;
unit = units(1);
set(hfig,'Units','centimeters');
set(hfig,'Position',[1,1,40,10]);
set(hfig,'PaperPositionMode','auto');
i = 1;

while unit~=-1,
    clf
    for i  = 1:3
        pf = pfs{i};
        ratemap = pf.plot(unit,'mazeMaskFlag',false);
        ratemap(isnan(ratemap)) = -1;
        for s = 1:numel(slices)
            sp(i,s) = axes('Units',spOpts.units,...
                           'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2)-4,...
                                       (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                        spOpts.width,...
                                        spOpts.height]...
                           );
            hold('on')
            
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap(:,:,slices(s)).*mask');    
            axis xy
            colormap([0,0,0;parula]);
            %caxis([-1,max(ratemap(:).*reshape(repmat(mask,[1,1,size(ratemap,3)]),[],1))]);
            caxis([-1,mRates(unit==units).*0.5]);
            title(num2str(round(pf.adata.bins{3}(slices(s))))) 
            text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
                 sprintf('%2.1f',max(max(ratemap(:,:,slices(s))))),'Color','w','FontWeight','bold','FontSize',10)
        if s==1, ylabel({pf.session.trialName,pf.parameters.states,['unit: ',num2str(unit)]});end
        end
    end

    FigName = ['pfs_unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end








if display,    
    set(0,'defaultAxesFontSize',8,...
          'defaultTextFontSize',8)

    autoincr = false;
    hfig = figure(38384);
    hfig.Units = 'centimeters';
    hfig.Position = [1,1,40,24];    
    hfig.PaperPositionMode = 'auto';

    while unit~=-1,
        sp = [];
        for i=1:numel(pfs),
            mrxz(i) =  max(max(max(pfs{i}.plot(unit,'xz',1,'mazeMaskFlag',false))));
        end

        Nmrxz = max(mrxz);
        for i=1:numel(pfs),
            subplot2(numel(pfs),2,i,1);
            pfs{i}.plot(unit,'xz',1,[0,mrxz],'mazeMaskFlag',false); 
            title([pfs{i}.session.trialName,':',pfs{i}.parameters.states,': ',num2str(unit)]);
            subplot2(numel(pfs),2,i,2);
            pfs{i}.plot(unit,'xy',1,[0,mrxz],'mazeMaskFlag',false); 
        end
        saveas(gcf,fullfile('/gpfs01/sirota/home/gravio/',...
                            'figures','SFN2014',...
                            ['pfsFly2_' Trial.filebase '-' num2str(unit) '.png']),'png');

        unit = figure_controls(hfig,unit,units,autoincr);
    end

end







pfs{s} = MTAApfs(Trial,units,'theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);




%% projecting the max rates along the x y and z axis 

%% setup figure paths
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016';
FigDir = ['pfs_flight_xz_yz_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));



%% Plot placefields across height


spOpts.width  = 2;
spOpts.height = 2;
spOpts.ny = 3;
spOpts.nx = 3;
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 0;
figOpts.position = [1,1,(spOpts.height+round(spOpts.padding/2))*spOpts.ny+round(spOpts.padding/2),...
                     (spOpts.width+round(spOpts.padding/2)) *spOpts.nx+figOpts.headerPadding+figOpts.footerPadding];


% rectangular mask
width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;
mask = repmat(mask,[1,1,pfs{1}.adata.binSizes(3)]);    

    

unit = units(1);
hfig = figure(393929);
autoincr = true;
unit = units(1);
set(hfig,'Units','centimeters');
set(hfig,'Position',[1,1,20,10]);
set(hfig,'PaperPositionMode','auto');
i = 1;

mdim = 'xyz';
neye = 1:3;

mRates = [];
for p = 1:3
    mRates(:,p) = pfs{p}.maxRate;
end
mRates = max(mRates,[],2);

unit = units(1);
i = 1;
while unit~=-1,
    clf
    for i  = 1:3
        pf = pfs{i};
        ratemap = pf.plot(unit,'mazeMaskFlag',false);
        ratemap = ratemap.*mask;
        ratemap(isnan(ratemap)) = -1; 
        for s = 1:3
            sp(i,s) = axes('Units',spOpts.units,...
                           'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2),...
                                       (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                        spOpts.width,...
                                        spOpts.height]...
                           );
            hold('on')
            dd = neye~=s;
            %mrm = sq(nanmean(ratemap,s));
            mrm = sq(nanmax(ratemap,[],s));
            imagesc(pf.adata.bins{find(dd,1,'first')},pf.adata.bins{find(dd,1,'last')},mrm');    
            axis xy
            colormap([0,0,0;parula]);
            caxis([-1,mRates(unit==units)]);
            title(mdim(dd))
            text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-900,...
                 sprintf('%2.1f',max(mrm(:))),'Color','w','FontWeight','bold','FontSize',10)
        if s==1, ylabel({pf.session.trialName,pf.parameters.states,['unit: ',num2str(unit)]});end
        end
    end

    FigName = ['pfs_unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end




% Flight of Ed10-20140822

MTAstartup('vr_exp');

QuickSessionSetup('Ed10-20140822.cof');
QuickTrialSetup('Ed10-20140822.cof');

s = MTASession.validate('Ed10-20140822.cof.all');

NeuronQuality(s)
figure();pZ(s,gcf,'t');
figure();pXY(s);

Trial = MTATrial.validate('Ed10-20140822.cof.all');

figure();pZ(Trial,gcf,'t');
figure();pXY(s);
pft = pfs_2d_theta(Trial,[],[],false); % overwrite

figure();
for u = pft.data.clu,
    clf();
    plot(pft,u);
    title(['unit: ',num2str(u),'  el:',pft.data.el(u)]);
    waitforbuttonpress();    
end

%compute_pfstats_bs(Trial,'Ed10-20140822.cof',)


Trial = MTATrial.validate('Ed10-20140822.cof.all');


Trial = MTATrial.validate('Ed10-20140822.rov.all');


pfkbs = {};

Trial = MTATrial.validate('Ed10-20140822.rov.stf');
fprintf('processing state: %s...\n',states{sts});
defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
defargs.overwrite = true;
defargs.states = 'velthresh&theta';
defargs.nNodes = 12;
defargs = struct2varargin(defargs);        
pfkbs{1} = MTAAknnpfs_bs(Trial,defargs{:});      

Trial = MTATrial.validate('Ed10-20140822.rov.vobj');
Trial = MTATrial.validate('Ed10-20140822.rov.blk');
Trial = MTATrial.validate('Ed10-20140822.rov.fly');


