function pfs = pfs_3d_states(Trial,varargin)
[stcMode,states,numIter,reportFig] = DefaultArgs(varargin,{'',{'theta','walk','rear','turn','pause','groom','sit'},1,1},1);


Trial= MTATrial.validate(Trial);    
units = Trial.spk.map(:,1);
nsts = numel(states);

% load theta periods
% if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end


%% Setup figure paths
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = ['pfs_3d_states_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));


%% compute 3d place fields for the theta state
for s = 1:numel(states),
    pfs{s} = MTAApfs(Trial,[],states{s},...
                     false,...
                     'numIter',1,...
                     'binDims',[20,20,20],...
                     'type','xyz',...
                     'SmoothingWeights',[2.2,2.2,2.2]);
end


if reportFig

    %% setup figure
    slices = 1:2:pfs{1}.adata.binSizes(end);
       
    spOpts.width  = 2;
    spOpts.height = 2;
    spOpts.ny = nsts;
    spOpts.nx = numel(slices);
    spOpts.padding = 1;
    spOpts.units = 'centimeters';


    figOpts.headerPadding = 2;
    figOpts.footerPadding = 0;
    figOpts.position = [1,1,(spOpts.width +spOpts.padding*2)*spOpts.nx,...
                            (spOpts.height+spOpts.padding*2)*spOpts.ny+...
                                    figOpts.headerPadding+...
                                    figOpts.footerPadding];


    width = pfs{1}.adata.binSizes(1);
    height = pfs{1}.adata.binSizes(2);
    radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
    mask(mask==0)=nan;
    



    %% Plot place fields sliced along z axis

    hfig = figure(393929);    
    hfig.Units    = 'centimeters';
    hfig.Position = figOpts.position;
    hfig.PaperPositionMode = 'auto';

    autoincr = true;
    unit = units(1);
    i = 1;
    while unit~=-1,
        clf
        for i = 1:nsts
            pf = pfs{i};
            ratemap = pf.plot(unit,'isCircular',false);
            ratemap(isnan(ratemap)) = -1;
            
            for s = 1:numel(slices)
                sp(i,s) = axes('Units',spOpts.units,...
                               'Position',[(spOpts.width+round(spOpts.padding/2))*(s)+round(spOpts.padding/2),...
                                    (spOpts.height+round(spOpts.padding/2))*(i)+round(spOpts.padding/2),...
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
                if s==1&&i==1, ylabel({Trial.filebase,states{i},['unit: ',num2str(unit)]});end
                if s==1&&i~=1, ylabel({states{i},['unit: ',num2str(unit)]});end            
            end
        end
        

        FigName = ['pfs_3d_states_unit-',num2str(unit)];
        %print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

        unit = figure_controls(hfig,unit,units,autoincr);    
    end

end