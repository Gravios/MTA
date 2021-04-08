function [pfs] = req20201117(Trial,unitsInt,tag,displayFlag,waitFlag,overwrite)
% req20201117(Trial)
%  Tags: interneuron ratemaps
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Project: MjgER2016: BhvPlaceCode
%  Description: Interneurons are behaviorally modulated. Plot 4d ratemap in xyhb space
%  Protocol: 
%    1. select interneurons
%    2. compute ratemaps
%    3. select fully occupied spatial bins
%    4. compute erpPCA

%Trial = MTATrial.validate('jg05-20120312.cof.all');
Trial = MTATrial.validate(Trial);
%Trial = MTATrial.validate('jg05-20120310.cof.all');
baseFigPath = '/storage/gravio/figures/analysis';
figDir = create_directory(...
    fullfile(baseFigPath,['ratemaps_xyhb_interneurons_2020/',Trial.filebase]));

%unitsInt = select_units(Trial,'int');


sampleRate = 16;
tper =[Trial.stc{'loc+rear&theta-groom-sit'}];
%tper =[Trial.stc{'theta-groom-sit'}];

xyzp = [];
if overwrite,
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);
    fet = fet_HB_pitchB(Trial,sampleRate);
    xyzp = xyz.copy();
    xyzp.data = [xyz(:,'hcom',1),xyz(:,'hcom',2),fet.data];
end

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.tag              = tag;
pargs.states           = tper;
pargs.units            = unitsInt;
pargs.numIter          = 1;
pargs.halfsample       = false;
pargs.overwrite        = overwrite;
pargs.boundaryLimits   = [-600,600;-600,600;-2,0.8;-0.8,2];
pargs.binDims          = [200,200,0.1,0.1];
pargs.SmoothingWeights = [0.5,0.5,1.8,1.8];
pargs.xyzp             = xyzp;
pargs.autoSaveFlag     = true;
%pargs.compute_pfs     = @PlotPF_Int;
pargs.compute_pfs     = @PlotPF;
pfsArgs = struct2varargin(pargs);
pfs = MTAApfs(Trial,pfsArgs{:});



if displayFlag,
    
    pft = pfs_2d_theta(Trial,unitsInt);
    [accg,tbins] = autoccg(Trial);

    nbins = 6;
    % SETUP figure
    hfig = figure();
    hfig.Units = 'Centimeters';
    hfig.PaperPositionMode = 'auto';
    hfig.Position = [0,0,17.5,18.5];
    % GENERATE axes
    hax = tight_subplot(nbins+1,nbins,0,0.1);

    for u = 1:numel(unitsInt);

        % PLOT theta place field
        axes(hax(1));
        plot(pft,unitsInt(u),1);
        clear_axes_ticks(gca);
        Lines(-500:200:500,[],'m');
        Lines([],-500:200:500,'m');    
        % PLOT ccg
        axes(hax(2));
        bar(tbins,accg(:,unitsInt(u)));axis tight;
        clear_axes_ticks(gca);

        % GET ratemap and get min max of
        ratemap = plot(pfs,unitsInt(u),'mean','',[],false);    
        maxrate = prctile(ratemap(:),95);
        minrate = prctile(ratemap(:),5);
        % REPORT min and max rate of the colormap used in the mosaic
        axes(hax(3));
        text(0.05,0.4,{'Min Rate:',num2str(minrate)});
        axes(hax(4));
        text(0.05,0.4,{'Max Rate:',num2str(maxrate)});
        % PLOT 4D fielnd    
        for x = 1:nbins,
            for y = 1:nbins,
                axes(hax((x)*nbins+y));
                imagescnan({pfs.adata.bins{3},pfs.adata.bins{4},sq(ratemap(y,nbins+1-x,:,:))'},... 
                           [minrate,maxrate],...
                           'colorMap',@jet);
                axis('xy');
                clear_axes_ticks(gca);
% $$$             xlim(pfs.adata.bins{3}([3,end-2]));
% $$$             ylim(pfs.adata.bins{4}([2,end-2]));       
            end
        end

        drawnow();
        % SAVE png and eps of figure
        figName = ['pfs','_',Trial.filebase,'_unit-',num2str(unitsInt(u))];
        %print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
        print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));

        if waitFlag,
            waitforbuttonpress();    
        end

        % CLEAR all axes
        ForAllSubplots('cla();');    
    end
    close(hfig)

end