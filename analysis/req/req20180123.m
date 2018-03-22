


% SET Analysis Parameters
version = '3';
overwrite = true;
display = true;
sessionListName = 'MjgER2016';
figDir = create_directory('/storage/gravio/figures/analysis/placefields_nonSpatialFeatures'); 
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta',            ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
             'theta-groom-sit'};

% LOAD session list
% LOAD Trials
% SELECT units for analysis
% LOAD theta placefields
sessionList = get_session_list(sessionListName);

Trials  = af(@(S)  MTATrial.validate(S),   sessionList);
          cf(@(T)  T.load('nq'),           Trials);
units   = cf(@(T)  select_placefields(T),  Trials);

pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);



% TRACK parameter counts 
numStates = numel(states);
numTrials = numel(Trials);


%% pfd_BPITCHxHPITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analDir = ['BPITCHxHPITCH-v',version];
create_directory(fullfile(figDir,analDir));

pfd = {};
pfindex = 1;
for tind = 1:numTrials,
    Trial = Trials{tind}; 

% LOAD vars
    xyz = preproc_xyz(Trial,'trb');
    pch = fet_HB_pitchB(Trial);
    drz = compute_drz(Trial,units{tind},pft{tind});%,pfstats);
    tper =[Trial.stc{'theta-groom-sit'}];
    tper.resample(xyz);

% COMPUTE HPITCH x BPITCH | DRZ [-0.5,0.5]
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units   = units{tind};
    pargs.numIter = 1001;
    pargs.halfsample = true;
    pargs.tag            = ['DRZxHBPITCHxBPITCH_v',version];
    pargs.boundaryLimits = [-2,2;-2,2];
    pargs.binDims        = [0.1,0.1];
    pargs.SmoothingWeights = [2,2];
    if overwrite,  
        pargs.overwrite = true;    
        pargs.xyzp = MTADxyz('data',pch.data,'sampleRate',xyz.sampleRate);
        drzState = {};
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            MTAApfs(Trial,pfsArgs{:});
        end
    end
    pargs.states    = 'tdrz';
    pargs.units     = units{tind};
    pargs.overwrite = false;
    pfsArgs = struct2varargin(pargs);
    pfd{tind,pfindex} = MTAApfs(Trial,pfsArgs{:});

% VISUALIZE HPITCH x BPITCH | DRZ [-0.5,0.5]
    if display,
        dspch = pch.copy();
        dspch.resample(5);
        dsxyz = xyz.copy();
        dsxyz.resample(5);

        hfig = figure(666002);
        hfig.Units = 'centimeters';
        %hfig.Position = [0.5,0.5,35,25];
        hfig.Position = [0.5,0.5,16,12];
        hfig.PaperPositionMode = 'auto';
        ny = 12;
        hax = gobjects([1,4]);
        for u = 1:numel(units{tind}), 
            clf();    
            maxPfsRate = max([pft{tind}.maxRate(units{tind}(u)),pfd{tind,pfindex}.maxRate(units{tind}(u),'mazeMaskFlag',false)]);

% PLOT placefield rate map
            hax(1) = subplot(221);  hold('on');  plot(pft{tind},units{tind}(u),'mean',true,maxPfsRate,false,0.99);
            plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
            xlabel('mm');  xlim([-500,500]);
            ylabel('mm');  ylim([-500,500]);
            title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
            hax(2) = subplot(222);  
            hold('on');  
            plot(pfd{tind,pfindex},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
            colorbar();    
            plot(dspch(drzState{u},2),...
                 dspch(drzState{u},1),'.m','MarkerSize',1),
            xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
            ylabel('body pitch (rad)');  ylim([-pi/2,2]);
            title('RateMap');

% PLOT placefield rate map
            hax(3) = subplot(223);  
            hold('on');  
            plot(pft{tind},units{tind}(u),'snr',true,[],false,0.99);
            plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
            xlabel('mm');  xlim([-500,500]);
            ylabel('mm');  ylim([-500,500]);
            title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

% PLOT Rate map PITCH x HEIGHT | DRZ[-0.5,0.5]
            hax(4) = subplot(224);  
            hold('on');  
            plot(pfd{tind,pfindex},units{tind}(u),'snr',true,5,false,0.85,false);
            plot(dspch(drzState{u},2),...
                 dspch(drzState{u},1),'.m','MarkerSize',1),
            xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
            ylabel('body pitch (rad)');  ylim([-2,pi/2]);
            title('SNR Map');

% FORMAT figure
            af(@(h) set(h,'Units','centimeters'),            hax);    
            af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax);
            af(@(h) set(h.Title,'Units','pixels'),           hax);
            af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);

% SAVE figure
            drawnow();
            figName = ['rateMap_BHPITCHxBPITCH_',version,'_',Trial.filebase,'_unit-',num2str(units{tind}(u))];
            print(hfig,'-depsc2',fullfile(figDir,analDir,[figName,'.eps']));        
            print(hfig,'-dpng',  fullfile(figDir,analDir,[figName,'.png']));
        end%for u
    end%if display
end%for tind



% COMPUTE erpPCA_HPITCHxBPITCH eigenvectors ------------------------------------------------------------------
rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,1));
rmaps = cat(2, rmaps{:});
si    = cf(@(p) p.data.si(:,:,1), pfd(:,1));
si = cat(2, si{:});
[~,sind] = sort(si);
[smnd,sind] = sort(sum(isnan(rmaps),1));
nind = ~isnan(rmaps(:,sind(1)));
bins = pfd{1,1}.adata.bins;


zrmaps = rmaps(:,sind(smnd>600));
zdims = size(zrmaps);
zrmaps(isnan(zrmaps)) = 0;
validDimsInds = sum(zrmaps==0,2)<zdims(2)/3;
zrmaps = zrmaps(validDimsInds,:);

[LU,LR,FSr,VT] = erpPCA(zrmaps',5);
V = LR;

reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';

hfig = figure(666003);clf();
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,20,6];
hfig.PaperPositionMode = 'auto';

% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] = bhv_contours(sessionListName,                              ... sessionListName
                                                  'fet_HB_pitchB',                                  ... featureSet
                                                  [1,2],                                            ... featureInd
                                                  {{linspace(-2,2,50),linspace(-2,2,50)}},          ... featureBin
                                                  'Ed05-20140529.ont.all',                          ... referenceTrial
                                                  {{'lloc+lpause&theta','hloc+hpause&theta',        ... states
                                                    'rear&theta'}},                                 ...
                                                  'wcr'                                             ... stateColors
);


nV = 5;
hax = gobjects([1,5]);
for i = 1:nV,
    hax(i) = subplot(1,nV,i);
    fpc = nan([zdims(1),1]);
    fpc(validDimsInds) = V(:,i);
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc,pfd))},[-0.5,3.5],'linear',false,[0,0,0],1,1);       % PRINT eigenvectors
    colorbar();                                           
    axis('xy');
    axis('tight');
    hold('on');
    for s = 1:numel(states),                              % OVERLAY state Contours
        copyobj(H{s},gca);
    end
    caxis([-0.5,3.5])
end
af(@(h) set(h,'Units','centimeters'),            hax);    
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
ForAllSubplots('daspect([1,1,1])');

figName = ['erpPCA_HPITCHxBPITCH_',version,'_',sessionListName];
print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));

% ERPPCA HPITCHxBPITCH END ---------------------------------------------------------------------------





%% erpPCA_HPITCHxBSPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analDir = ['BPITCHxBSPEED-v',version];
create_directory(fullfile(figDir,analDir));
display = false;
overwrite = false;
pfindex = 2;

for tind = 1:numTrials,
    Trial = Trials{tind}; %15,16,17,18

    if display || overwrite,
        xyz = preproc_xyz(Trial,'trb');
        xyz.filter('RectFilter');
        pft = pfs_2d_theta(Trial,'overwrite',false);
        pch = fet_HB_pitchB(Trial);
        vxy = xyz.vel(['spine_lower'],[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);
        fet = pch.copy('empty');
        fet.label = 'fet_VP';       
        fet.data = [pch(:,2),vxy(:)];
        drz = compute_drz(Trial,units{tind},pft{tind});%,pfstats);
        tper =[Trial.stc{'theta-groom-sit-rear'}];
        tper.resample(xyz);
    end

    % CE HPITCHxVEL
    pargs = get_default_args('MjgER2016','MTAApfs','struct');
    pargs.units   = units{tind};
    pargs.numIter = 101;
    pargs.halfsample = true;
    pargs.tag            = ['DRZxHBPITCHxVEL_v',version];
    pargs.boundaryLimits = [-2,2;-2,2];
    pargs.binDims        = [0.1,0.1];
    pargs.SmoothingWeights = [2,2];
    if overwrite,  
        pargs.overwrite = true;    
        pargs.xyzp = MTADxyz('data',fet.data,'sampleRate',xyz.sampleRate);
        drzState = {};
        for u = 1:numel(units{tind});
            dper = MTADepoch([],[],ThreshCross(-0.5<drz(:,u)&drz(:,u)<0.5,0.5,1),...% SELECT periods where drz
                             xyz.sampleRate,xyz.sync.copy(),xyz.origin,'TimePeriods','sts',[],'tdrz','d');
            drzState{u} = dper&tper;
            pargs.units  = units{tind}(u);
            pargs.states = drzState{u};
            pfsArgs = struct2varargin(pargs);
            MTAApfs(Trial,pfsArgs{:});
        end
    end
    pargs.states    = 'tdrz';
    pargs.units     = units{tind};
    pargs.overwrite = false;
    pfsArgs = struct2varargin(pargs);
    pfd{tind,pfindex} = MTAApfs(Trial,pfsArgs{:});


    if display,
        dspch = pch.copy();
        dspch.resample(5);
        dsxyz = xyz.copy();
        dsxyz.resample(5);        
        dsvxy = dsxyz.vel(['spine_lower'],[1,2]);
        dsvxy.data(dsvxy.data<1e-3) = 1e-3;
        dsvxy.data = log10(dsvxy.data);
        dsfet = dspch.copy('empty');
        dsfet.label = 'fet_VP';       
        dsfet.data = [dspch(:,2),dsvxy(:)];
        
        

        hfig = figure(666003);
        hfig.Units = 'centimeters';
        %hfig.Position = [0.5,0.5,35,25];
        hfig.Position = [0.5,0.5,16,12];
        hfig.PaperPositionMode = 'auto';
        ny = 12;
        hax = gobjects([1,4]);
        for u = 1:numel(units{tind}), 
            clf();    
            maxPfsRate = max([pft{tind}.maxRate(units{tind}(u)),pfd{tind,pfindex}.maxRate(units{tind}(u),'mazeMaskFlag',false)]);
            
% PLOT theta placefield rate map            
% PLOT theta placefield SNR map            
            hax(1) = subplot(221);  hold('on');  plot(pft{tind},units{tind}(u),'mean',true,maxPfsRate,false,0.99);
            plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
            xlabel('mm');  xlim([-500,500]);
            ylabel('mm');  ylim([-500,500]);
            title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

            hax(3) = subplot(223);  
            hold('on');  
            plot(pft{tind},units{tind}(u),'snr',true,[],false,0.99);
            plot(dsxyz(drzState{u},'nose',1),dsxyz(drzState{u},'nose',2),'.m','MarkerSize',1),
            xlabel('mm');  xlim([-500,500]);
            ylabel('mm');  ylim([-500,500]);
            title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

            
% PLOT rate map HPITCH x BSPEED | DRZ[-0.5,0.5]
% PLOT SNR  map HPITCH x BSPEED | DRZ[-0.5,0.5]
            hax(2) = subplot(222);  
            hold('on');  
            plot(pfd{tind,pfindex},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
            colorbar();    
            plot(dsfet(drzState{u},1),...
                 dsfet(drzState{u},2),'.m','MarkerSize',1),
            xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
            ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
            title('RateMap');

            hax(4) = subplot(224);  
            hold('on');  
            plot(pfd{tind,pfindex},units{tind}(u),'snr',true,5,false,0.85,false);
            plot(dsfet(drzState{u},1),...
                 dsfet(drzState{u},2),'.m','MarkerSize',1),
            xlabel('head-body pitch (rad)');  xlim([-2,pi/2]);
            ylabel('BL speed (log10(cm/s))');  ylim([-2,2]);
            title('SNR Map');
            
% SET figure parameters
            af(@(h) set(h,'Units','centimeters'),            hax);    
            af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
            af(@(h) set(h.Title,'Units','pixels'),           hax);
            af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);
            % SAVE figure
            drawnow();
            figName = ['rateMap_BHPITCHxVEL_',version,'_',Trial.filebase,'_unit-',num2str(units{tind}(u))];
            print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
            print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
        end%for u

    end%if display

end%for tind

bins = pfd{1,pfindex}.adata.bins;

% CREATE SELCETED VECTOR SUBSPACE
rmaps = cf(@(p) p.data.rateMap(:,:,1), pfd(:,pfindex));
si    = cf(@(p) p.data.si(:,:,1),      pfd(:,pfindex));
rmaps = cat(2, rmaps{:});
si    = cat(2, si{:}   );
[~,sind] = sort(si);
[smnd,sind] = sort(sum(isnan(rmaps),1));
nind = ~isnan(rmaps(:,sind(1)));
zrmaps = rmaps(:,sind(smnd>800));
zdims = size(zrmaps);                               % NOTE the dimensions of the original vector space
zrmaps(isnan(zrmaps)) = 0;                          % SET all nan valued elements to zeros
validDimsInds = sum(zrmaps==0,2)<zdims(2)/2;        % SELECT subspace {1/3 non-zero samples} 
zrmaps = zrmaps(validDimsInds,:);                   % REDUCE to selected subspace

% DECOMPOSE 
[LU,LR,FSr,VT] = erpPCA(zrmaps',5);                 % COMPUTE covariance-based PCA with Varimax rotation
V = LR;

% helper function to reshape eigenvectors
reshape_eigen_vector = @(V,pfd) reshape(V(:,1),pfd{1}.adata.binSizes')';



% LOAD Behavioral state contours
[stateContourMaps,stateContourHandles] = bhv_contours(sessionListName,                              ... sessionListName
                                                  'fet_HB_HPS',                                     ... featureSet
                                                  [1,2],                                            ... featureInd
                                                  {{linspace(-2,2,50),linspace(-2,2,50)}},          ... featureBin
                                                  'Ed05-20140529.ont.all',                          ... referenceTrial
                                                  {{'lloc+lpause&theta','hloc+hpause&theta',        ... states
                                                    'rear&theta'}},                                 ...
                                                  'wcr'                                             ... stateColors
);

hfig = figure(666003);clf();
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,20,6];
hfig.PaperPositionMode = 'auto';



nV = 5;
hax = gobjects([1,5]);
for i = 1:nV,
    hax(i) = subplot(1,nV,i);
    fpc = nan([zdims(1),1]);
    fpc(validDimsInds) = V(:,i);
    imagescnan({bins{:},abs(reshape_eigen_vector(fpc,pfd(:,pfindex)))},[-0.5,3.5],'linear',false,[0,0,0],1,1);       % PRINT eigenvectors
    colorbar();                                           
    %title(sprintf('PC%i Var:%3.2f',i,dS(i)));             % APPEND rank and variance as title
    axis('xy');
    axis('tight');
    hold('on');
    for s = 1:numel(stateContourHandles),                              % OVERLAY state Contours
        copyobj(stateContourHandles{s},hax(i));
    end
    caxis([-0.5,3.5])
end
af(@(h) set(h,'Units','centimeters'),            hax);    
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
ForAllSubplots('daspect([1,1,1])');

figName = ['erpPCA_HPITCHxVEL_',version,'_',sessionListName];
print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));

% ERPPCA HPITCHxSPEED END -------------------------------------------------------------------------


%% pfd matrix for anton %%

pfSizes = [pft{1}.adata.binSizes';...
           pfd{1,1}.adata.binSizes';...
           pfd{1,2}.adata.binSizes'];

pfsMat = nan([sum(cellfun(@numel,units)),max(pfSizes)]);

tpftMat = [];
for tind = 1:numTrials
    for uind = units{tind},
        tpftMat(end+1,:,:) = pft{tind}.plot(uind,'mean',false,[],false,0.9);
    end
end
tpfdHPBPMat = [];
for tind = 1:numTrials
    for uind = units{tind},
        tpfdHPBPMat(end+1,:,:) = pfd{tind,1}.plot(uind,'mean',false,[],false,0.9);
    end
end
tpfdHPBVMat = [];
for tind = 1:numTrials
    for uind = units{tind},
        tpfdHPBVMat(end+1,:,:) = pfd{tind,2}.plot(uind,'mean',false,[],false,0.9);
    end
end


map = [0,0,0,0,0];
uCumCount = cumsum(cellfun(@numel,units));
uCount = cellfun(@numel,units);
for tind = 1:numTrials
    map = cat(1,map,[[map(end,1)+1:uCumCount(tind)]',repmat(tind,[uCount(tind),1]),pft{tind}.data.clu',pft{tind}.data.el',pft{tind}.data.elClu']);
end
map(1,:) = [];


mapLabels = {'index','Trial','clu','el','elClu'};
pfsTheta.data = tpftMat;
pfsTheta.bins = pft{1}.adata.bins;
pfdHpitchBpitch.data = tpfdHPBPMat;
pfdHpitchBpitch.bins = pfd{1,1}.adata.bins;
pfdHpitchBspeed.data = tpfdHPBVMat;
pfdHpitchBspeed.bins = pfd{1,2}.adata.bins;


pfsData =  reshape(permute(pfsTheta.data,[2,3,1]),[],1,sum(uCount));
pfsData(1:prod(pfd{1,1}.adata.binSizes),2,:) = reshape(permute(pfdHpitchBpitch.data,[2,3,1]),[],1,sum(uCount));
pfsData(1:prod(pfd{1,2}.adata.binSizes),3,:) = reshape(permute(pfdHpitchBspeed.data,[2,3,1]),[],1,sum(uCount));
pfsDataLabels = {'rate',{'theta','HPITCHxBPITCH','HPITCHxBSPEED'},'units'};



save('/storage/share/Projects/BehaviorPlaceCode/placefields_behaviorfields_for_anton.mat',...
     'map','mapLabels','pfsTheta','pfdHpitchBpitch','pfdHpitchBspeed','pfsData','pfsDataLabels');


figure,imagesc(pfdHpitchBspeed.bins{:},sq(pfdHpitchBspeed.data(402,:,:))');



%% pfd plot parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARTS  

analDir = ['parts-v',version];
create_directory(fullfile(figDir,analDir));

for tind = 1:numTrials,
    hfig = figure(666004);
    hfig.Units = 'centimeters';
    hfig.Position = [0.5,0.5,20,12];
    hfig.PaperPositionMode = 'auto';
    ny = 12;
    hax = gobjects([1,6]);
    for u = 1:numel(units{tind}), 
        clf();    
        maxPfsRate = max([maxRate(pft{tind}  ,units{tind}(u),false,'mean'),...
                          maxRate(pfd{tind,1},units{tind}(u),false,'mean'),...
                          maxRate(pfd{tind,2},units{tind}(u),false,'mean')]);

% PLOT placefield RATE map
% PLOT placefield SNR map        
        hax(1) = subplot(231);  hold('on');  plot(pft{tind},units{tind}(u),'mean',true,maxPfsRate,false,0.99);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta Place Field, unit:',num2str(units{tind}(u))]);

        hax(2) = subplot(234);  
        hold('on');  
        plot(pft{tind},units{tind}(u),'snr',true,[],false,0.99);
        xlabel('mm');  xlim([-500,500]);
        ylabel('mm');  ylim([-500,500]);
        title(['Theta SNR Field, unit:',num2str(units{tind}(u))]);

        
% PLOT Rate map PITCH x HEIGHT | DRZ [-0.5,0.5]
% PLOT SNR  map PITCH x HEIGHT | DRZ [-0.5,0.5]        
        hax(3) = subplot(232);  
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-2,2]);
        ylabel('body pitch (rad)');  ylim([-2,2]);
        title('RateMap');

        hax(4) = subplot(235);  
        hold('on');  
        plot(pfd{tind,1},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-2,2]);
        ylabel('body pitch (rad)');  ylim([-2,2]);
        title('SNR Map');

% PLOT Rate map HPITCH x BSPEED | DRZ [-0.5,0.5]
% PLOT SNR  map HPITCH x BSPEED | DRZ [-0.5,0.5]        
        hax(5) = subplot(233);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'mean',true,maxPfsRate,false,0.85,false);
        colorbar();    
        xlabel('head-body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body pitch (rad)');  ylim([-2,2]);
        title('RateMap');

        hax(6) = subplot(236);  
        hold('on');  
        plot(pfd{tind,2},units{tind}(u),'snr',true,5,false,0.85,false);
        xlabel('head-body pitch (rad)');  xlim([-pi/2,2]);
        ylabel('body speed log(cm/s)');  ylim([-2,2]);
        title('SNR Map');
        
        
        
        
        % FORMAT figure
        af(@(h) set(h,'Units','centimeters'),            hax);    
        af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax(1:2));
        af(@(h) set(h,'Position',[h.Position(1:2),1.5,1.5]), hax(3:end));        
        af(@(h) set(h.Title,'Units','pixels'),           hax);
        af(@(h) set(h.Title,'Position',h.Title.Position+[0,20,0]),  hax);

        % SAVE figure
        drawnow();
        figName = ['rateMap_BHPITCHxBPITCHxBSPEED_',version,'_',Trials{tind}.filebase,'_unit-',num2str(units{tind}(u))];
        print(hfig,'-depsc2',fullfile(figDir,analDir,[figName,'.eps']));        
        print(hfig,'-dpng',  fullfile(figDir,analDir,[figName,'.png']));
    end%for u

end%for tind



