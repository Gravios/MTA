
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)


%% Setup jg05-20120311
% modified entry within SessionList for jg05-20120311

% setup session
QuickSessionSetup(SessionList('jg05-20120311'));

% load session
s = MTASession('jg05-20120311');

% correct marker swaps
ERCOR_fillgaps_RidgidBody(s,1234000);

% plot timeseries for z axis
pZ(s);
% plot position within xy plane
pXY(s);





% create trial wrapper for full session
trialName = 'jg05-20120311.cof.all';
%QuickTrialSetup(trialName)
Trial = MTATrial.validate(trialName);

% plot timeseries for z axis
pZ(Trial);
% plot position within xy plane
pXY(Trial);

Trial = labelTheta(Trial,[],68,1);
Trial = labelTheta(Trial);

Stc = Trial.stc.copy;


xyz = Trial.load('xyz');
spk = Trial.spk.copy;
spk.create(Trial,xyz.sampleRate);
figure,plot(spk(10),'.')



% 3.5
trialName = 'jg05-20120311.cof.spk';
%QuickTrialSetup(trialName)
Trial = MTATrial.validate(trialName);
Trial = labelTheta(Trial);


% plot timeseries for z axis
pZ(Trial);
% plot position within xy plane
pXY(Trial);


%Trial.stc = Stc.copy;
%Trial.stc.load(Trial); 

overwrite = 0;
units = [];

% compute 3d place fields for the theta state
pf = MTAApfs(Trial,units,'theta',...
             overwrite,...
             'numIter',1,...
             'binDims',[20,20,20],...
             'type','xyz',...
             'SmoothingWeights',[2.2,2.2,2.2]);



% setup figure paths
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016';
FigDir = ['pfs_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));


slices = 1:2:17;

spOpts.width  = 2;
spOpts.height = 2;
spOpts.ny = 3;
spOpts.nx = numel(slices);
spOpts.padding = 2;
spOpts.units = 'centimeters';
figOpts.units = 'centimeters';
figOpts.headerPadding = 2;
figOpts.footerPadding = 0;
figOpts.position = [1,1,(spOpts.height+round(spOpts.padding/2))*spOpts.ny+figOpts.headerPadding+figOpts.footerPadding,...
                     (spOpts.width+round(spOpts.padding/2)) *spOpts.nx+round(spOpts.padding/2)];

width = pf.adata.binSizes(1);
height = pf.adata.binSizes(2);
radius = round(pf.adata.binSizes(1)/2)-find(pf.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;
    

units = Trial.spk.map(:,1);

hfig = figure(393929);    
autoincr = true;
unit = units(1);
set(hfig,'Units','centimeters');
set(hfig,'Position',[1,1,40,6]);
set(hfig,'PaperPositionMode','auto');
i = 1;


while unit~=-1,
    clf
    ratemap = pf.plot(unit,'isCircular',false);
    ratemap(isnan(ratemap)) = -1;
    for s = 1:numel(slices)
        sp(i,s) = axes('Units',spOpts.units,...
                       'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2),...
                            (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                            spOpts.width,...
                            spOpts.height]...
                       );
        hold('on')
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap(:,:,slices(s)).*mask');    
        axis xy
        colormap([0,0,0;parula]);
        caxis([-1,max(ratemap(:).*reshape(repmat(mask,[1,1,size(ratemap,3)]),[],1))]);
        title(num2str(round(pf.adata.bins{3}(slices(s)))))
        text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(max(ratemap(:,:,slices(s))))),'Color','w','FontWeight','bold','FontSize',10)
    end

    FigName = ['pfs_theta_unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

    unit = figure_controls(hfig,unit,units,autoincr);    
end

