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
pf = MTAApfs(Trial,                               ... MTATrial
             [],                                  ... units
             'theta',                             ... state label
             false,                               ... overwrite
             'numIter',1,                         ... number of iterations
             'binDims',[100,100,30],              ... feature dimensions
             'type','xyz',                        ... string to identify save file
             'SmoothingWeights',[0.8,0.8,0.8]);   ... smoothing weights (bins per std)



if reportFig

    interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                   linspace(-500,500,50),...
                                   linspace( -60,360,60)}},...
                          'nanMaskThreshold', 0.1,...
                          'methodNanMap',     'cubic',...
                          'methodRateMap',    'cubic');

    %% setup figure
    slices = 1:4:numel(interpParPfs.bins{end});%min(cellfun(@(x) x.binSizes(end),pfs));

    
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



    %% Plot place fields sliced along z axis

    hfig = figure(393929);    
    hfig.Units    = 'centimeters';
    hfig.Position = [1,1,40,6];
    hfig.PaperPositionMode = 'auto';
    interpParPfs = struct('bins',{{linspace(-500,500,50),...
                                   linspace(-500,500,50),...
                                   linspace( -60,360,60)}},...
                          'nanMaskThreshold', 0.1,...
                          'methodNanMap',     'cubic',...
                          'methodRateMap',    'cubic');

    
    autoincr = true;
    unit = units(1);
    i = 1;
    while unit~=-1,
        clf
        ratemap = pf.plot(unit,'mean','mazeMaskFlag',false,'interpPar',interpParPfs);
        ratemap(isnan(ratemap)) = -1;
        for s = 1:numel(slices)
            sp(i,s) = axes('Units',spOpts.units,...
                           'Position',[(spOpts.width+round(spOpts.padding/2))*(s+1)+round(spOpts.padding/2),...
                                (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                spOpts.width,...
                                spOpts.height]...
                           );
            hold('on')
            imagescnan({interpParPfs.bins{1},interpParPfs.bins{2},ratemap(:,:,slices(s))'});    
            axis xy
            colormap([0,0,0;parula]);
            caxis([-1,max(ratemap(:))]);
            title(num2str(round(interpParPfs.bins{3}(slices(s)))))
            axis tight            
            text(interpParPfs.bins{1}(round(numel(interpParPfs.bins{1}).*0.7)),...
                 interpParPfs.bins{2}(round(numel(interpParPfs.bins{2}).*0.9)),...
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