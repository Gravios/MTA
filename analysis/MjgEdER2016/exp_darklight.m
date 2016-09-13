


% Set overwrite permisions
overWriteSessions = false;
overWriteTrials   = false;
overWriteStc      = false;

% Set Default figure parameters
set(0,'defaultAxesFontSize',8,...
      'defaultTextFontSize',8)


% Retrive Session and Trial lists
sessionList = 'Ed10_cof';
trialList   = 'Ed10_darklight';


% Set Figure directory 
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
mkdir(OwnDir);


% Compile Sessions and Trials from data 
Slist = SessionList(sessionList);
if overWriteSessions, 
    QuickSessionSetup(Slist,[],[],false); 
    QuickTrialSetup(Slist);
end
Tlist = SessionList(trialList);
if overWriteTrials,   
    QuickTrialSetup(Tlist,'overwrite',true); 
end


% Dis
s = MTASession('Ed10-20140817');
xyz = s.load('xyz');
pXY(s)
pZ(s)
PlotSessionErrors(s)
s.spk.create(s)


Trial = MTATrial.validate('Ed10-20140817.cof.all');
Trial = labelBhvBasic(Trial);
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end

% $$$ fxyz = xyz.copy;
% $$$ fxyz.filter('ButFilter',3,2.4,'low');
% $$$ vxy = fxyz.vel([1,7],[1,2]);
% $$$ vxy.data(vxy.data<1e-3) = 1e-3;
% $$$ vxy.data = log10(vxy.data);
% $$$ 
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ 
% $$$ ind = vxy(:,1)>0.5;
% $$$ figure
% $$$ subplot(121),hist(ang(ind,'head_back','head_front',2),100)
% $$$ subplot(122),hist(xyz(ind,'head_back',3),100)


Trial = MTATrial.validate('Ed10-20140817.cof.dark');
Trial.load('stc','auto_wbhr');



units = 1:73;

pfs = {};
pfs{1} = MTAApfs(Trial,units,'velHthresh&theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);





Trial = MTATrial.validate('Ed10-20140817.cof.light');
Trial = labelBhvBasic(Trial);
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end


pfs{2} = MTAApfs(Trial,units,'velHthresh&theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);



%% Plot placefields across height
FigDir = 'exp_darklight_pfs_3d';
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

width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;
    
figHnum = 393929;
autoincr = false;
unit = units(1);
hfig = figure(figHnum);clf
set(hfig,'units',figOpts.units)
set(hfig,'Position',[1,1,40,10]);%figOpts.position)
set(hfig,'PaperPositionMode','auto');


while unit~=-1,
    for i  = 1:2
        pf = pfs{i};
        ratemap = pf.plot(unit,'isCircular',false);
        ratemap(isnan(ratemap)) = -1;
        for s = 1:numel(slices)
            if unit~=units(1),
                axes(sp(i,s));
                cla;
                hold('on');
            else
                sp(i,s) = axes('Units',spOpts.units,...
                               'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2),...
                                    (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                    spOpts.width,...
                                    spOpts.height]...
                               );
                hold('on')
            end
            
            imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap(:,:,slices(s)).*mask');    
            axis xy
            colormap([0,0,0;parula]);
            caxis([-1,max(ratemap(:).*reshape(repmat(mask,[1,1,size(ratemap,3)]),[],1))]);
            title(num2str(round(pf.adata.bins{3}(slices(s)))))
            text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
                 sprintf('%2.1f',max(max(ratemap(:,:,slices(s))))),'Color','w','FontWeight','bold','FontSize',10)
            if i==1&s==1,
                ylabel(['Unit: ' num2str(unit)]);
            end
            
        end
    end

    
    
    FigName = ['pfs_3d_unit',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
    print(hfig,'-dpng',fullfile(OwnDir,FigDir,[FigName,'.png']));

    fprintf('%s: fig unit %i\n',FigName,unit)n
    unit = figure_controls(hfig,unit,units,autoincr);
    %fprintf('complete\n')
end





