function pf = pfs_2d_theta_traj(Trial,varargin)
%function pf = pfs_2d_theta(Trial,varargin)
%
%  varargin:
%    xyzTemporalShift
%    figureFormats
%    overwrite
%


% DEFARGS --------------------------------------------------------------------------------------
defargs = struct('xyzTemporalShift',                0,                                       ...
                 'marker',                          'hcom',                                  ...                 
                 'figureFormats',                   {{}},                                    ...
                 'overwrite',                       0                                        ...
);
[xyzTemporalShift,marker,figureFormats,overwrite] = DefaultArgs(varargin,defargs,'--struct');
% END DEFARGS -----------------------------------------------------------------------------------



% DEFVARS ----------------------------------------------------------------------
Trial= MTATrial.validate(Trial);    
units = Trial.spk.map(:,1);
if isempty(Trial.stc.gsi('t')),Trial = labelTheta(Trial);end
% END DEFVARS ----------------------------------------------------------------------



% MAIN -------------------------------------------------------------------------    
%% compute 2d place fields for the theta state

% RESTRICT computation to theta periods without groom or sit
try, 
    state = 'theta-sit-groom';
    Trial.stc{state};
catch
    state = 'theta';
end

% LOAD xyz
xyz = preproc_xyz(Trial,'trb');
xyz.data = circshift(sq(xyz(:,marker,[1,2])),-round(xyzTemporalShift*xyz.sampleRate));

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.tag       = ['theta2dtrajshift',num2str(xyzTemporalShift),'m-',marker];
defargs.numIter   = 1;
defargs.units     = units;
defargs.states    = state;
defargs.overwrite = overwrite;
defargs.xyzp      = xyz;
defargs.trackingMarker = marker;
pfsArgs           = struct2varargin(defargs);
pf = MTAApfs(Trial,pfsArgs{:});



if ~isempty(figureFormats),
% CREATE figure
    FigDir = create_directory(fullfile('/storage/gravio/figures/placefields/pfs_2d_theta_traj/',Trial.filebase));
    
    hfig = figure(393929);    
    hfig.Units    = 'centimeters';
    hfig.Position = [1,1,14,8];
    hfig.PaperPositionMode = 'auto';
    
    autoincr = true;
    unit = units(1);
    i = 1;
    hax = gobjects([1,1]);
    while unit~=-1,
% PLOT placefield
        hax(1) = subplot(1,1,1);
        plot(pf,unit),
        text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',pf.maxRate(unit)),'Color','w','FontWeight','bold','FontSize',10)

% PRINT figure
        FigName = ['pfs_2d_theta_traj_unit-',num2str(unit)];                
        for format = figureFormats,
            switch format{1},
              case 'eps',
                print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
              case 'png',
                print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));
            end
        end
        
% ITERATE to next unit
        unit = figure_controls(hfig,unit,units,autoincr);    
    end

end

% END MAIN -------------------------------------------------------------------------    