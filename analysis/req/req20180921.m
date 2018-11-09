function req20180921(Trial)
% req20180921
%  Tags: interneuron ratemaps
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Project: MjgER2016:placefields
%  Description: Interneurons are behaviorally modulated. Plot 4d ratemap in xyhb space
%  Protocol: 
%    1. select interneurons
%    2  compute ratemaps


Trial = MTATrial.validate(Trial);
%Trial = MTATrial.validate('jg05-20120310.cof.all');

waitFlag = false;
figDir = create_directory(['/storage/gravio/figures/analysis/ratemaps_xyhb_interneurons/',...
                    Trial.filebase]);

unitsInt = select_units(Trial,'int');

sampleRate = 16;
xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
fet = fet_HB_pitchB(Trial,sampleRate);
tper =[Trial.stc{'theta-groom-sit'}];
xyzp = xyz.copy();
xyzp.data = [xyz(:,'hcom',1),xyz(:,'hcom',2),fet.data];

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.tag              = 'interneurons_xyhb';
pargs.states           = tper;
pargs.units            = unitsInt;
pargs.numIter          = 1;
pargs.halfsample       = false;
pargs.overwrite        = true;
pargs.boundaryLimits   = [-500,500;-500,500;-2.25,1;-0.5,2];
pargs.binDims          = [100,100,0.2,0.2];
pargs.SmoothingWeights = [0.8,0.8,1.1,1.1];
pargs.xyzp             = xyzp;
pargs.autoSaveFlag     = true;
pfsArgs = struct2varargin(pargs);
pfs = MTAApfs(Trial,pfsArgs{:});

pft = pfs_2d_theta(Trial,unitsInt);
[accg,tbins] = autoccg(Trial);

% SETUP figure
hfig = figure();
hfig.Units = 'Centimeters';
hfig.PaperPositionMode = 'auto';
hfig.Position = [0,0,17.5,18.5];
% GENERATE axes
hax = tight_subplot(10,10,0,0.1);
for u = 1:numel(unitsInt);
% PLOT theta place field
    axes(hax(1));
    plot(pft,unitsInt(u),1);
    clear_axes_ticks(gca);
% PLOT ccg
    axes(hax(2));
    bar(tbins,accg(:,unitsInt(u)));axis tight;
    clear_axes_ticks(gca);

% GET ratemap and get min max of
    ratemap = plot(pfs,unitsInt(u),'mean','',[],false);    
    maxrate = prctile(ratemap(:),98);
    minrate = prctile(ratemap(:),2);
% REPORT min and max rate of the colormap used in the mosaic
    axes(hax(3));
    text(0.05,0.4,{'Min Rate:',num2str(minrate)});
    axes(hax(4));
    text(0.05,0.4,{'Max Rate:',num2str(maxrate)});
% PLOT 4D field    
    for x = 1:9,
        for y = 2:10,
            axes(hax((x)*10+y));
            imagescnan({pfs.adata.bins{3},pfs.adata.bins{4},sq(ratemap(y,11-x,:,:))'},...
                       [minrate,maxrate],...
                       'colorMap',@jet);
            axis('xy');
            clear_axes_ticks(gca);
            xlim(pfs.adata.bins{3}([3,end-2]));
            ylim(pfs.adata.bins{4}([2,end-2]));       
        end
    end
    drawnow();
% SAVE png and eps of figure
    figName = ['pfs','_',Trial.filebase,'_unit-',num2str(unitsInt(u))];
    print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
    
    if waitFlag,
        waitforbuttonpress();    
    end
    
% CLEAR all axes
    ForAllSubplots('cla();');    
end
close(hfig)
