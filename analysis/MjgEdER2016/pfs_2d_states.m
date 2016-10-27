function pfs = pfs_2d_states(Trial,varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('stcMode',       'NN0317R',                                                     ...
                 'states',        {{'loc','lloc','hloc','rear','pause','lpause','hpause'}},      ...
                 'reportFig',     true,                                                          ...
                 'tag',           '',                                                            ...
                 'overwrite',     false                                                          ...
);%-------------------------------------------------------------------------------------------------

[stcMode,states,reportFig,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');

% modify states
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states = cat(2,{'theta-sit-groom'},states);




% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
%
%    stcMode 
%    states 
if isempty(tag),
    tag = DataHash(struct('stcMode',stcMode,'states',{states}));
end
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

Trial= MTATrial.validate(Trial);    

% load labeled behavior
try,
    Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
catch err
    disp(err)
    Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
end


% Reduce clu list based on theta pfs max rate
pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;
units = select_units(Trial,18);
units = units(mrt(pft.data.clu(units))>1);

nsts = numel(states);


%% Setup figure paths
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = ['pfs_2d_states_',tag,'_',Trial.filebase];
mkdir(fullfile(OwnDir,FigDir));



%% compute 3d place fields for the theta state
pfs = {};
for s = 1:nsts
    defargs = get_default_args_MjgEdER2016('MTAApfs','struct');
    defargs.units = units;
    defargs.states = states{s};
    defargs.numIter = 1;
    defargs = struct2varargin(defargs);        
    pfs{s} = MTAApfs(Trial,defargs{:});      
end


if reportFig

    %% setup figure
    spOpts.width  = 2;
    spOpts.height = 2;
    spOpts.ny = 1;
    spOpts.nx = nsts;
    spOpts.padding = 2;
    spOpts.units = 'centimeters';
    figOpts.units = 'centimeters';
    figOpts.headerPadding = 2;
    figOpts.footerPadding = 0;
    figOpts.position = [1,...
                        1,...
                        (spOpts.width+round(spOpts.padding/2)+1)*spOpts.nx,...
                        (spOpts.height+round(spOpts.padding/2))*spOpts.ny+...
                          figOpts.headerPadding+figOpts.footerPadding];

    width = pfs{1}.adata.binSizes(1);
    height = pfs{1}.adata.binSizes(2);
    radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
    mask(mask==0)=nan;
    

    maxPfsRate = [];
    for s = 1:nsts,
        maxPfsRate(s,:) = pfs{s}.maxRate(units);
    end

    %% Plot place fields sliced along z axis
    hfig = figure(393929);    
    hfig.Units    = figOpts.units;
    hfig.Position = figOpts.position;
    hfig.PaperPositionMode = 'auto';

    autoincr = true;
    unit = units(1);
    i = 1;
    while unit~=-1,
        clf
        for s = 1:nsts
            pf = pfs{s};
            ratemap = pf.plot(unit,'isCircular',false);
            ratemap(isnan(ratemap)) = -1;
            
            sp(i,s) = axes('Units',spOpts.units,...
                           'Position',[2+(spOpts.width+round(spOpts.padding/2))*(s)+round(spOpts.padding/2),...
                                (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                spOpts.width,...
                                spOpts.height]...
                           );
            hold('on')
            imagesc(pf.adata.bins{1},pf.adata.bins{2},(ratemap.*mask)');    
            axis xy
            colormap([0,0,0;parula]);
            caxis([-1,max(maxPfsRate(:,unit==units))]);
            title([states{s} ' ' num2str(maxPfsRate(s,unit==units))]);
% $$$             text(pf.adata.bins{1}(end)-300,pf.adata.bins{2}(end)-50,...
% $$$                  sprintf('%2.1f',max(max(ratemap))),'Color','w','FontWeight','bold','FontSize',8)
        end
        
        FigName = ['pfs_2d_states_unit-',num2str(unit)];
        FigInfo = uicontrol('Parent',hfig,...
                            'Style','text',...
                            'String',{FigName,Trial.filebase},...
                            'Units','centimeters',...
                            'Position',[.5,(spOpts.height+round(spOpts.padding/2))*spOpts.ny-3,4,3]);

        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

        unit = figure_controls(hfig,unit,units,autoincr);    
    end

end