function pf = pfs_3d_theta(Trial,varargin)
[numIter,reportFig] = DefaultArgs(varargin,{1,1},1);

Trial= MTATrial.validate(Trial);    
units = Trial.spk.map(:,1);

% load theta periods
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end


%% Setup figure paths
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = ['pfs_3d_theta_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));


%% compute 3d place fields for the theta state
pf = MTAApfs(Trial,[],'theta',...
             false,...
             'numIter',1,...
             'binDims',[20,20,20],...
             'type','xyz',...
             'SmoothingWeights',[2.2,2.2,2.2]);



if reportFig

    %% setup figure
    slices = 1:2:pf.adata.binSizes(end);%min(cellfun(@(x) x.binSizes(end),pfs));

    
% $$$     if strcmp(Trial.maze.shape,'rectangle'),
% $$$         spOpts.width  = 4;
% $$$     else
        spOpts.width  = 2;
% $$$     end
    
    if strcmp(Trial.maze.shape,'rectangle'),
        spOpts.height  = 1;
    else
        spOpts.height  = 2;
    end
    %spOpts.height = 2;
    spOpts.ny = 3;
    spOpts.nx = numel(slices);
    spOpts.padding = 2;
    spOpts.units = 'centimeters';


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
    if strcmp(Trial.maze.shape,'rectangle'),
        mask(:) = 1;
    end    



    %% Plot place fields sliced along z axis

    hfig = figure(393929);    
    hfig.Units    = 'centimeters';
    hfig.Position = [1,1,40,6];
    hfig.PaperPositionMode = 'auto';

    autoincr = true;
    unit = units(1);
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
            axis tight            
            text(pf.adata.bins{1}(round(numel(pf.adata.bins{1}).*0.7)),...
                 pf.adata.bins{2}(round(numel(pf.adata.bins{2}).*0.9)),...
                 sprintf('%2.1f',max(max(ratemap(:,:,slices(s))))),...
                 'Color','w',...
                 'FontWeight','bold',...
                 'FontSize',10)
            if s==1, ylabel({Trial.filebase,['unit: ',num2str(unit)]});end
        end


        FigName = ['pfs_3d_theta_unit-',num2str(unit)];
        %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

        unit = figure_controls(hfig,unit,units,autoincr);    
    end

end