function pfs = pfs_2d_states(Trial,varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',          [],                                                           ...
                 'stcMode',        'msnn_ppsvd_raux',                                            ...
                 'states',         {{'loc&theta','lloc&theta','hloc&theta','rear&theta',         ...
                                     'pause&theta','lpause&theta','hpause&theta'}},              ...
                 'reportFig',      false,                                                        ...
                 'tag',            '',                                                           ...
                 'overwrite',      false,                                                        ...
                 'pfsArgsOverride',[]                                                            ...
);
[units,stcMode,states,reportFig,tag,overwrite, pfsArgsOverride] =                                ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
%
%    stcMode 
%    states 
if isempty(tag),
    tag = DataHash(struct('stcMode',stcMode,'pfsArgsOverride',pfsArgsOverride));
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

Trial= MTATrial.validate(Trial);    

% LOAD labeled behavior
if ~strcmp(Trial.stc.mode,stcMode),
    try,
        Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
    catch err
        disp(err)
        Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
    end
end


% REDUCE clu list based on theta pfs max ratea
if isempty(units)
    units = select_placefields(Trial,18);
elseif ischar(units)&strcmp(units,'all')
    units = Trial.spk.map(:,1);
end

nsts = numel(states);

%% Setup figure paths
OwnDir = '/storage/gravio/figures/placefields/';
FigDir = ['pfs_2d_states_',Trial.filebase,'_',tag];
mkdir(fullfile(OwnDir,FigDir));


%% compute 3d place fields for the theta state
pfs = {};
for s = 1:nsts
    if ischar(states{s}),
        disp(['MTA:analysis:placefields:pfs_2d_states:processing:',Trial.filebase,':',states{s}]);
    else
        disp(['MTA:analysis:placefields:pfs_2d_states:processing:',Trial.filebase,':',states{s}.label]);
    end
    
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units = units;
    pargs.states = states{s};
    
    if ~isempty(pfsArgsOverride) && isstruct(pfsArgsOverride),
        for f = fieldnames(pfsArgsOverride)'
            pargs.(f{1}) = pfsArgsOverride.(f{1});
        end
    end
    
    pargs.overwrite = overwrite;    
    %pargs.tag = DataHash({states{s},tag});
    pfsArgs = struct2varargin(pargs);        
    pfs{s} = MTAApfs(Trial,pfsArgs{:});      
end



% FIGURE 
if reportFig
    
% COMPUTE unit auto correlograms 
    [accg,tbins] = autoccg(MTASession.validate(Trial.filebase));
    
% LOAD unit quality information
    Trial.load('nq');


% SETUP figure options
    spOpts.width  = 2;
    spOpts.height = 2;
    spOpts.ny = 1;
    spOpts.nx = nsts+1;
    spOpts.padding = 1;
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
    maxPfsSI = [];
    for s = 1:nsts,
        for u = 1:numel(units),
            maxPfsRate(s,u) = pfs{s}.maxRate(units(u));
            maxPfsSI(s,u) = pfs{s}.data.si(pfs{s}.data.clu==units(u));
        end        
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
            ratemap = pf.plot(unit,'mazeMaskFlag',false);
            ratemap(isnan(ratemap)) = -1;
            
            sp(i,s) = axes('Units',spOpts.units,...
                           'Position',[3+(spOpts.width+round(spOpts.padding/2))*(s)+round(spOpts.padding/2),...
                                (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                spOpts.width,...
                                spOpts.height]...
                           );
            hold('on')
            imagesc(pf.adata.bins{1},pf.adata.bins{2},(ratemap.*mask)');    
            axis xy
            colormap([0,0,0;parula]);
            caxis([-1,max(maxPfsRate(:,unit==units))]);
            xlabel(states{s})
            title({['max rate: ',num2str(round(maxPfsRate(s,unit==units),2))],...
                   [' si: ',num2str(round(maxPfsSI(s,unit==units),2))]});
            if s ~= 1,
                set(gca,'XTickLabels',{});
                set(gca,'YTickLabels',{});
            else
                axes('Units',spOpts.units,...
                     'Position',[3+(spOpts.width+round(spOpts.padding/2))*(nsts+1)+round(spOpts.padding/2),...
                                 (spOpts.height+round(spOpts.padding/2))*(i-1)+round(spOpts.padding/2),...
                                 spOpts.width,...
                                 spOpts.height]...
                     );
                 bar(tbins,accg(:,unit));axis tight;            
            end
            
        end
        
        FigName = ['pfs_2d_states_unit-',num2str(unit)];
        FigInfo = uicontrol('Parent',hfig,...
                            'Style','text',...
                            'String',{FigName,Trial.filebase,...
                                      ['stcMode: ',Trial.stc.mode],...
                                      ['eDist:   ',num2str(Trial.nq.eDist(unit))],...
                                      ['Refrac:  ',num2str(log10(Trial.nq.Refrac(unit)))],...
                                      ['SNR:     ',num2str(Trial.nq.SNR(unit))],...
                                      ['AmpSym:  ',num2str(Trial.nq.AmpSym(unit))],...
                                      ['SpkWidthR:  ',num2str(Trial.nq.SpkWidthR(unit))]...
                                     },...
                            'Units','centimeters',...
                            'Position',[.2,(spOpts.height+round(spOpts.padding/2))*spOpts.ny-3,6,4]);

        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

        unit = figure_controls(hfig,unit,units,autoincr);    
        pause(0.2);
    end

end
