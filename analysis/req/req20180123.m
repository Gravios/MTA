% p(s|phi,z)


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';



FigDir = create_directory('/storage/gravio/figures/placefields_nonSpatialFeatures'); 

% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
             'theta-groom-sit'};

numStates = numel(states);

cf(@(t) t.load('nq'), Trials);

Trial = Trials{15}; %15,16,17,18

units = select_placefields(Trial);

xyz = preproc_xyz(Trial,'trb');
pft = pfs_2d_theta(Trial,'overwrite',false);    


ang = create(MTADang,Trial,xyz);

% $$$ 
% $$$ ind = [Trial.stc{'x+p+r'}];
% $$$ figure();
% $$$ subplot(121);plot(ang(ind,'hcom','nose',2),ang(ind,'spine_middle','spine_upper',2),'.');
% $$$ subplot(122);plot(ang(ind,'hcom','nose',2),ang(ind,'spine_middle','hcom',2),'.');


pch = fet_HB_pitch(Trial);
map_to_reference_session(pch,Trial,pitchReferenceTrial);    


drz = compute_drz(Trial,units,pft);%,pfstats);


tper = [Trial.stc{'theta-groom-sit'}];
tper.resample(xyz);

overwrite = true;

drzState = {};
for u = 1:numel(units);
    dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...
                     xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
    drzState{u} = dper&tper;    
end


% CE DRZxPITCHxHEIGHT
% $$$ pargs = get_default_args('MjgER2016','MTAApfs','struct');
% $$$ pargs.units      = units;
% $$$ pargs.numIter = 1;
% $$$ pargs.halfsample = false;
% $$$ pargs.tag            = 'DRZxPITCHxHEIGHT';
% $$$ pargs.boundaryLimits = [-pi/2-0.3,pi/2+0.3;-10,350];
% $$$ pargs.binDims        = [0.1,10];
% $$$ if overwrite,
% $$$ pargs.overwrite  = true;    
% $$$ pargs.xyzp = MTADxyz('data',[pch(:,3),xyz(:,'nose',3)],'sampleRate',xyz.sampleRate);
% $$$ for u = 1:numel(units);
% $$$     pargs.units  = units(u);        
% $$$     pargs.states = drzState{u};
% $$$     pfsArgs = struct2varargin(pargs);
% $$$     pfs_dp = MTAApfs(Trial,pfsArgs{:});
% $$$ end
% $$$ end
% $$$ pargs.states    = 'tdrz';
% $$$ pargs.units     = units;
% $$$ pargs.overwrite = false;
% $$$ pfsArgs = struct2varargin(pargs);
% $$$ pfs_dp = MTAApfs(Trial,pfsArgs{:});
% $$$ 

% CE DRZxPITCHxHEIGHT
% $$$ pargs = get_default_args('MjgER2016','MTAApfs','struct');
% $$$ pargs.units      = units;
% $$$ pargs.numIter = 1001;
% $$$ pargs.halfsample = true;
% $$$ pargs.tag            = 'DRZxHPITCHxBPITCH';
% $$$ pargs.boundaryLimits = [-pi/2,pi/2;-pi/2,pi/2];
% $$$ pargs.binDims        = [0.1,0.1];
% $$$ pargs.SmoothingWeights = [1.5,1.5];
% $$$ if overwrite,
% $$$ pargs.overwrite  = true;    
% $$$ pargs.xyzp = MTADxyz('data',[ang(:,'hcom','nose',2),ang(:,'spine_middle','spine_upper',2)],'sampleRate',xyz.sampleRate);
% $$$ for u = 1:numel(units);
% $$$     pargs.units  = units(u);        
% $$$     pargs.states = drzState{u};
% $$$     pfsArgs = struct2varargin(pargs);
% $$$     pfs_dp = MTAApfs(Trial,pfsArgs{:});
% $$$ end
% $$$ end
% $$$ pargs.states    = 'tdrz';
% $$$ pargs.units     = units;
% $$$ pargs.overwrite = false;
% $$$ pfsArgs = struct2varargin(pargs);
% $$$ pfs_bh = MTAApfs(Trial,pfsArgs{:});
% $$$ 


% CE DRZxPITCHxHEIGHT
pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units      = units;
pargs.numIter = 1001;
pargs.halfsample = true;
pargs.tag            = 'DRZxHBPITCHxBPITCH';
pargs.boundaryLimits = [-pi/2,pi/2;-pi/2,pi/2];
pargs.binDims        = [0.1,0.1];
pargs.SmoothingWeights = [1.5,1.5];
if overwrite,
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',[circ_dist(pch(:,3),pch(:,1)),pch(:,1)],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units  = units(u);        
        pargs.states = drzState{u};
        pfsArgs = struct2varargin(pargs);
        pfs_hbh = MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.states    = 'tdrz';
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfs_hbh = MTAApfs(Trial,pfsArgs{:});



dsxyz = xyz.copy();
dsxyz.resample(12);
dsang = create(MTADang,Trial,dsxyz);

hfig = figure(666002);
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,35,16];
hfig.PaperPositionMode = 'auto';
ny = 12;
hax = gobjects([1,4]);
for u = 1:numel(units), 
    clf();    
% PLOT placefield rate map
    hax(1) = subplot(221);  hold('on');  plot(pft,units(u));  colorbar();
    plot(xyz(drzState{u},'nose',1),xyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta Place Field, unit:',num2str(units(u))]);
    
% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(2) = subplot(222);  
    hold('on');  
    plot(pfs_bh,units(u),'isCircular',false);
    colorbar();
    
    plot(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2),'.m','MarkerSize',1),
    xlabel('head pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('head height (mm)');  ylim([-pi/2,pi/2]);
    title('RateMap');

% PLOT placefield rate map
    hax(3) = subplot(223);  hold('on');  plot(pft,units(u),'snr');  colorbar();
    plot(xyz(drzState{u},'nose',1),xyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta SNR Field, unit:',num2str(units(u))]);
    
% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(4) = subplot(224);  
    hold('on');  
    plot(pfs_bh,units(u),'snr',false,5,'isCircular',false);%1.5*pfs_bh.maxRate(units(u),false)
    colorbar();
    plot(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2),'.m','MarkerSize',1),
    xlabel('head pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('head height (mm)');  ylim([-pi/2,pi/2]);
    title('SNR Map');
    
% SAVE figure
    drawnow();
    FigName = ['rateMap_HPITCHxBPITCH','_',Trial.filebase,'_unit-',num2str(units(u))];
    print(hfig,'-depsc2',fullfile(FigDir,[FigName,'.eps']));        
    print(hfig,'-dpng',  fullfile(FigDir,[FigName,'.png']));
end



bhvccg = {}; sper = {};
for sts = 1:numStates,
    [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5,units,1);
    sper{sts}{1}(sper{sts}{1}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
    sper{sts}{2}(sper{sts}{2}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
end

bhvccgp = {}; sperp = {};
for sts = 1:numStates,
    [bhvccgp{sts},sperp{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5,units,4);
    sperp{sts}{1}(sperp{sts}{1}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
    sperp{sts}{2}(sperp{sts}{2}>ceil(size(xyz,1)./xyz.sampleRate.*1250)-1) = [];
end

figure();
plot(bhvccg{4},109,1)


eds = linspace(-pi/2,pi/2,31);

figure,
for u = 1:numel(units),
    clf();
mrate = reshape(mean(sq(pfs_hbh.data.rateMap(:,u,:)),2,'omitnan'),31,31);
subplot(221); hold('on');
    imagescnan({eds,eds,mrate'});
    axis('xy');
    plot(circ_dist(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2)),...
         dsang(drzState{u},'spine_middle','spine_upper',2),...
         '.m','MarkerSize',1),
    xlim([eds([1,end])]);
    ylim([eds([1,end])]);

scount = reshape(sum(~isnan(sq(pfs_hbh.data.rateMap(:,u,:))),2),31,31);
subplot(222); 
    imagescnan({eds,eds,scount'});
    axis('xy');    

mask = scount>990;
mrate(~mask(:)) = nan;
subplot(224); 
    imagescnan({eds,eds,mrate'});
    axis('xy');    
    hold('on');
    plot(circ_dist(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2)),...
         dsang(drzState{u},'spine_middle','spine_upper',2),...
         '.m','MarkerSize',1),
    
subplot(223); 
    hold('on');
    plot(pft,units(u),'mean',false,[],true,true);
    plot(xyz(drzState{u},'nose',1),xyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    colorbar();
    title(['Theta Place Field, unit:',num2str(units(u))]);
waitforbuttonpress();
end
    
%



hfig = figure(666002);
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,35,25];
hfig.PaperPositionMode = 'auto';
ny = 12;
hax = gobjects([1,4]);
for u = 1:numel(units), 
    clf();    
% PLOT placefield rate map
    hax(1) = subplot(221);  hold('on');  plot(pft,units(u));  colorbar();
    plot(xyz(drzState{u},'nose',1),xyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta Place Field, unit:',num2str(units(u))]);
    
% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(2) = subplot(222);  
    hold('on');  
    plot(pfs_hbh,units(u),'isCircular',false,'sigThresh',0.99);
    colorbar();    
    plot(circ_dist(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2)),...
         dsang(drzState{u},'spine_middle','spine_upper',2),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('body pitch (rad)');  ylim([-pi/2,pi/2]);
    title('RateMap');

% PLOT placefield rate map
    hax(3) = subplot(223);  
    hold('on');  
    plot(pft,units(u),'snr');  colorbar();
    plot(xyz(drzState{u},'nose',1),xyz(drzState{u},'nose',2),'.m','MarkerSize',1),
    xlabel('mm');  xlim([-500,500]);
    ylabel('mm');  ylim([-500,500]);
    title(['Theta SNR Field, unit:',num2str(units(u))]);
    
% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
    hax(4) = subplot(224);  
    hold('on');  
    plot(pfs_hbh,units(u),'snr',false,5,'isCircular',false,'sigThresh',0.99);
    colorbar();
    plot(circ_dist(dsang(drzState{u},'hcom','nose',2),dsang(drzState{u},'spine_middle','spine_upper',2)),...
         dsang(drzState{u},'spine_middle','spine_upper',2),'.m','MarkerSize',1),
    xlabel('head-body pitch (rad)');  xlim([-pi/2,pi/2]);
    ylabel('body pitch (rad)');  ylim([-pi/2,pi/2]);
    title('SNR Map');
    
% SAVE figure
    drawnow();
    FigName = ['rateMap_BHPITCHxBPITCH','_',Trial.filebase,'_unit-',num2str(units(u))];
    print(hfig,'-depsc2',fullfile(FigDir,[FigName,'.eps']));        
    print(hfig,'-dpng',  fullfile(FigDir,[FigName,'.png']));
end

