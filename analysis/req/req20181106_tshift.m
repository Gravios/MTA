
% ANALYSIS DESCRIPTION -----------------------------------------------------------------------------

% SAG - Secondary Analysis Goal
% TAG - Tertiary Analysis Goal


% PAG : Primary Analysis Goal
%    Deterimine the prefered theta phase for behavior state depedent rate changes associated with spatial
%    tuning curves
%
% VAR - Analysis variables
%    spike time (res)
%    head position (xyz)
%    local field potential (lfp)
%
% TRN - Transformations
%    FFT of lfp -> gmHMM of mean power within 6-12Hz band -> thetaPeriods (tper).
%    HILBERT of 6-12 Hz bandpass filtered lfp -> thetaPhase
%    DISCRETIZE head position at spike times durng theta states -> theta spike occupancy map.
%    DISCRETIZE head position durng theta states -> theta position occupancy map.
%    ELEMENT wise division of theta spike occupancy map and theta position occpancy map -> spatial rate map.
%    SELECTION of spike times within the 90th percentile within fixed radius around center of spatial rate
%        map -> place restricted selection of spike times during theta state.
%    CONVERT head markers positions from Cartesian to polar coordinate system.
%    
% MOD : Analysis Model 
%    meanFiringRate(thetaPhase, headPitch)
% CMD : Causal Model
%     {Sensory Information} EC3 -> CA1(theta[pi,2*pi])
%     {Auto Completion}     CA3 -> CA1(theta[0,pi])
%
% 
% PFSTAGS : hptp : head pitch, theta phase
%           hpef : head pitch, egocentric forward

% END ANALYSIS DESCRIPTION -------------------------------------------------------------------------


% ANALYSIS MAIN ------------------------------------------------------------------------------------
global MTA_PROJECT_PATH;

MjgER2016_load_data();

% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};

thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;

dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',...
                        ['req20181106_tshift-data-',...
                    DataHash({[sessionList.sessionName],sampleRate,states}),'.mat']);


%[pfd] = req20180123_ver5(Trial);


for t = [1:19,21:23];
%for t = 1:23
%for t = t
    %t = 20;    
    Trial = Trials{t}; 
    unitSubset = units{t};        
    subjectId = regexp(Trial.name,'^(\w*)-','tokens');
    subjectId = subjectId{1}{1};

% $$$     switch subjectId,
% $$$       case 'jg05'        
% $$$         pft = MTAApfs(Trial,unitSubset,[],[],'CA1thetaCA3inputPhase');
% $$$       otherwise  
        pft = pfs_2d_theta(Trial,unitSubset);
% $$$     end

    stc = Trial.load('stc','msnn_ppsvd_raux');
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
    mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
    mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                              atan2(diff(xyz(:,{'hcom','nose'},2),1,2),...
                                    diff(xyz(:,{'hcom','nose'},1),1,2)));
    try
        lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    catch err
        lfp = load(Trial,'lfp',sessionList(t).thetaRef);
    end
    
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    

% $$$     tvec = circshift(fxyz(:,'hcom',[1,2]),-50)-fxyz(:,'hcom',[1,2]);
% $$$     tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
% $$$     tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);

% $$$     tvecb = -circshift(fxyz(:,'hcom',[1,2]),50)-fxyz(:,'hcom',[1,2]);
% $$$     tvecb = sq(bsxfun(@rdivide,tvecb,sqrt(sum(tvecb.^2,3))));
% $$$     tvecb = cat(3,tvecb,sq(tvecb)*[0,-1;1,0]);
    
    hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
    hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
    hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
    [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);
    
% $$$     vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
% $$$     vxy.data(vxy.data<1e-3) = 1e-3;
% $$$     vxy.data = log10(vxy.data);

% $$$     headAngle = atan2(diff(fxyz(:,{'hcom','nose'},2),1,2),diff(fxyz(:,{'hcom','nose'},1),1,2));
% $$$     vang = [circ_dist(circshift(headAngle,-25),headAngle),circ_dist(circshift(headAngle,-75),headAngle),...
% $$$             circ_dist(circshift(headAngle,-125),headAngle),circ_dist(circshift(headAngle,-200),headAngle)];
% $$$     vangb = [-circ_dist(circshift(headAngle,25),headAngle),-circ_dist(circshift(headAngle,75),headAngle),...
% $$$              -circ_dist(circshift(headAngle,125),headAngle),-circ_dist(circshift(headAngle,200),headAngle)];
    
    stcm = stc2mat(stc,xyz,states);

    nq = get(Trial.load('nq'),'nq');
    edist = nq.eDist(unitSubset)';
    
    spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    
    
    [~,sind] = sort(pft.data.clu);
    si = subsref(pft.data.si(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));
    spar = subsref(pft.data.spar(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));    

% $$$     pch = fet_HB_pitch(Trial,sampleRate,false,'trb');
% $$$     pch.data  = pch.data(:,3);
% $$$     pch.name  = 'headPitch_thetaPhase';
% $$$     pch.label = 'fet_hptp';
% $$$     pch.key   = 'p';


    pchb = fet_HB_pitchB(Trial,sampleRate,false,'trb');        
    %pchb.data  = pch.data(:,1);
    


    aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);

    for shift = -250:50:250,
        
        pargs = struct('units',              unitSubset,                           ...
                       'states',             'theta',                              ...
                       'overwrite',          false,                                ...
                       'tag',                ['hbpptbpFS1v3_ts',num2str(shift)],   ...
                       'binDims',            [0.1,pi/8,0.1],                       ...
                       'SmoothingWeights',   [2,1.1,2],                            ...
                       'type',               'xyz',                                ...
                       'spkShuffle',         false,                                ...
                       'posShuffle',         false,                                ...
                       'numIter',            1,                                    ...
                       'xyzp',               [],                                   ...
                       'boundaryLimits',     [-1.8,1;-pi,pi;-0.8,2],               ...
                       'bootstrap',          false,                                ...
                       'halfsample',         false,                                ...
                       'compute_pfs',        @PlotPFCirc,                          ...
                       'autoSaveFlag',       false,                                ...
                       'spk',                spk                                   ...
                       );
        
        pfs = MTAApfs(Trial,'tag',pargs.tag);
        pfs.purge_savefile();
        pfs = Trial;
        for unit = unitSubset,
            [mxr,mxp] = pft.maxRate(unit);        
            pargs.xyzp = copy(pchb);
            pargs.xyzp.data = [circshift(pchb(:,1),shift),...
                               phz(:,spk.map(spk.map(:,1)==unit,2)),...
                               circshift(pchb(:,2),shift)];
            pargs.units = unit;
            pargs.states = MTADepoch([],                                                ...
                                     [],                                                ...
                                     ThreshCross(aper & abs(drz(:,unit==unitSubset))<0.8   ...
                                                 & abs(ddz(:,unit==unitSubset))<250,  ...
                                                 0.5,1),                           ...
                                     sampleRate,pargs.xyzp.sync.copy(),                       ...
                                     pargs.xyzp.origin,'TimePeriods','sts',[],'tdrz','d');
            
            pfsArgs = struct2varargin(pargs);
            pfs = MTAApfs(pfs,pfsArgs{:});    
            if unit==unitSubset(1),
                pfs.save();
            end
        end
        pfs.save();
        
    end

    
end


Trial = Trials{t}; 
unitSubset = units{t};        


tshifts = -250:50:250;
pfs = {};
for shift = -250:50:250,
    pfs{end+1} = MTAApfs(Trial,'tag',['hbpptbpFS1v3_ts',num2str(shift)]);
end

phzOrder = [9:16,1:8];
nanColor = []
hfig = figure();


nanColor = [0,0,0];
clf();
fax = axes('Position',[0,0,1,1],'Visible','off','FontSize',8,'LineWidth',1);
ylim([0,1]);xlim([0,1]);
hax = tight_subplot(numel(tshifts),numel(pfs{1}.adata.bins{2}),0.01,0.1);
u = unitSubset(1);
while u ~=-1,
    for s = 1:numel(tshifts),
        rmap = plot(pfs{s},u,1,'colorbar',[],false);
        if s == 1,
            rmax = prctile(rmap(nniz(rmap(:))),99.9);
        end
        if isnan(rmax), rmax = 1;end
        for i = 1:numel(phzOrder),
            if s==1, title(num2str(u));end
            hx = hax(i+(s-1).*numel(phzOrder));
            axes(hx);
            imagescnan({pfs{s}.adata.bins{[1,3]},sq(rmap(:,phzOrder(i),:))'},[0,rmax],...
                       'nanRGB',nanColor,'colorMap',@parula);
            axis('xy');
            hx.YTickLabels = {};
            hx.XTickLabels = {};            
            if s==6, 
                axes(fax);
                rectangle('Position',hx.Position+[-0.0005,-0.0005,0.001,0.001],'EdgeColor','m','LineWidth',2);
            end                        
        end
    end
    u = figure_controls(hfig,u,unitSubset,false,[],[]);
    ForAllSubplots('cla');
end










% $$$ 
% $$$ [pfd] = req20180123_ver5(Trial);
% $$$ 
% $$$ 
% $$$ t = 18;
% $$$ Trial = Trials{t}; 
% $$$ unitSubset = units{t};
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbp');
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbpF');
% $$$ pfs = MTAApfs(Trial,'tag','hbpptbpHS');
% $$$ pfe = MTAApfs(Trial,'tag','ptpefel');
% $$$ pfd = req20180123_ver5(Trial);
% $$$ pft = pfs_2d_theta(Trial,unitSubset);
% $$$ [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
% $$$ [ddz] = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);
% $$$ xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
% $$$ stc = Trial.load('stc','msnn_ppsvd_raux');    
% $$$ stcm = stc2mat(stc,xyz,states);    
% $$$ aper = stcm(:,1)==1 & ~any(stcm(:,[7,8]),2);
% $$$ 
% $$$ 
% $$$ pch = fet_HB_pitchB(Trial,[],false,'trb');        
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ hfig = figure();
% $$$ clf();
% $$$ hax = tight_subplot(3,13,0.01,0.1);
% $$$ u = unitSubset(1);
% $$$ while u ~=-1,
% $$$     rmap = plot(pfs,u,'mean','colorbar',[],false,0.8);
% $$$     %rmap = plot(pfs,u,1,'colorbar',[],false);    
% $$$     rmax = prctile(rmap(nniz(rmap(:))),98)+2;
% $$$     if isnan(rmax), rmax = 1;end
% $$$     axes(hax(1));
% $$$     plot(pft,u,1,'',[0,rmax],false,'colorMap',@jet);
% $$$     title(num2str(u));
% $$$ % $$$     axes(hax(14));
% $$$ % $$$     plot(pfta,u,1,'',[0,rmax],false,'colorMap',@jet);
% $$$     axes(hax(2));
% $$$     ind = aper & abs(drz(:,u==unitSubset))<0.8 & abs(ddz(:,u==unitSubset))<250;
% $$$     plot(xyz(ind,'hcom',1),xyz(ind,'hcom',2),'.r');
% $$$     xlim([-500,500]);
% $$$     ylim([-500,500]);
% $$$     axes(hax(3));        
% $$$     plot(pfd{1},u,1,'',[0,rmax],false,'colorMap',@jet);
% $$$     title('bhv space');
% $$$     
% $$$     
% $$$     for i = 1:10,
% $$$         axes(hax(i+3));
% $$$ 
% $$$         imagescnan({pfs.adata.bins{[1,3]},sq(rmap(:,mod(i+4,10)+1,:))'},[0,rmax],'colorMap',@jet);
% $$$         axis('xy');
% $$$         title(num2str(pfs.adata.bins{2}(mod(i+4,10)+1)));
% $$$     end
% $$$     %xlim([-1.5,0.5])
% $$$     
% $$$     if sum(rmap(:))~=0,
% $$$         rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))==10;
% $$$         zmap = reshape(permute(rmap,[1,3,2]),[],10);
% $$$         zmap(repmat(rmapNind(:),[1,10])&isnan(zmap)) = 0;    
% $$$         zmap(repmat(~rmapNind(:),[1,10])) = nan;
% $$$         zmap(zmap<0) = 0;
% $$$         zmap(isnan(zmap(:,1)),:) = [];
% $$$         
% $$$         [LU,LR,FSr,VT] = erpPCA(zmap',4);
% $$$         for v = 1:4,
% $$$             axes(hax(16+v));
% $$$             evec = nan([size(rmap,1),size(rmap,3)]);
% $$$             evec(rmapNind) = LR(:,v);
% $$$             imagesc(pfs.adata.bins{[1,3]},evec');
% $$$             axis('xy');        
% $$$             caxis(repmat(max(abs(caxis())),[1,2]).*[-1,1])
% $$$         end
% $$$         axes(hax(22));
% $$$         plot(VT(:,4),'-+');
% $$$         axes(hax(23));    
% $$$         imagesc(FSr([6:10,1:5],:)');
% $$$     end
% $$$ 
% $$$     rmap = plot(pfe,u,1,'colorbar',[],false);
% $$$     for i = 1:10,
% $$$         axes(hax(i+29));
% $$$ 
% $$$         imagescnan({pfe.adata.bins{[1,3]},sq(rmap(:,mod(i+4,10)+1,:))'},[0,rmax],'colorMap',@jet);
% $$$         axis('xy');
% $$$         title(num2str(pfs.adata.bins{2}(mod(i+4,10)+1)));
% $$$     end
% $$$ 
% $$$     axes(hax(27));
% $$$     colorbar();
% $$$     caxis([0,rmax]);
% $$$     colormap('jet');
% $$$     
% $$$     u = figure_controls(hfig,u,unitSubset,false,[],[]);
% $$$     ForAllSubplots('cla');
% $$$     
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for u = unitSubset
% $$$         subplot(1,3,3);    
% $$$     plot(pfd{1},u,1,'',[],false,'colorMap',@jet);
% $$$     title('bhv space');    
% $$$     
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ [pfd] = req20180123_ver5(Trial); 
% $$$ 
% $$$ 
% $$$ pfs = MTAApfs(Trial,'tag','hbptpefel');
% $$$ 
% $$$ hfig = figure();
% $$$ hax = tight_subplot(10,13,0.01,0.1);
% $$$ u = unitSubset(1);
% $$$ while 1,
% $$$ %for u = unitSubset
% $$$     try
% $$$     rmap = plot(pfs,u,1,'colorbar',[],false);
% $$$     rmax = prctile(rmap(nniz(rmap(:))),99.9)+4;
% $$$     axes(hax(1));
% $$$     plot(pft,u,1,'colorbar',[0,rmax],false,'colorMap',@jet);
% $$$     title(num2str(u));
% $$$     axes(hax(2));    
% $$$     plot(pfd{1},u,1,'colorbar',[0,rmax],false,'colorMap',@jet);
% $$$     title('bhv space');
% $$$     
% $$$ for i = 1:10,
% $$$     for j = 1:10
% $$$         ind = sub2ind([13,10],i+2,j);
% $$$         
% $$$     axes(hax(ind));
% $$$     imagescnan({pfs.adata.bins{[3,4]},sq(rmap(11-j,mod(i+4,10)+1,:,:))'},[0,rmax],'colorMap',@jet);
% $$$     axis('xy');
% $$$     title(num2str(pfs.adata.bins{2}(mod(i+4,10)+1)));
% $$$     Lines(0,[],'m');
% $$$     Lines([],0,'m');
% $$$ end    
% $$$ end
% $$$ end
% $$$ u = figure_controls(hfig,u,unitSubset,false,[],[])
% $$$ ForAllSubplots('cla');
% $$$ end
% $$$ 
% $$$ 
% $$$ pfa = MTAApfs(Trial,'tag','hptp');
% $$$ pfe = MTAApfs(Trial,'tag','hpef');
% $$$ 
% $$$ pft = pfs_2d_theta(Trial);
% $$$ 
% $$$ figure,
% $$$ for u = unitSubset,
% $$$     clf();
% $$$     subplot(241);
% $$$     plot(pft,u,1,'colorbar',[0,10],false,'colorMap',@jet);
% $$$     title(num2str(u));
% $$$     subplot(242);
% $$$     plot(pfd{1},u,1,'colorbar',[0,10],false,'colorMap',@jet);
% $$$     title('bhv space');
% $$$     subplot(243);
% $$$     plot(pfe,u,1,'colorbar',[0,10],false,'colorMap',@jet);
% $$$     ylabel('ego front/back')
% $$$     xlabel('head pitch');
% $$$     subplot(244);    
% $$$     plot(pfa,u,1,'colorbar',[0,10],false,'colorMap',@jet);    
% $$$     ylabel('theta phase')
% $$$     xlabel('head pitch');
% $$$     hax = subplot(248);    
% $$$     plot(pfa,u,1,'colorbar',[0,10],false,'colorMap',@jet);    
% $$$     ylabel('theta phase')
% $$$     xlabel('head pitch');
% $$$     drawnow();
% $$$     hax.Position(2) =     hax.Position(2) +0.135;
% $$$     drawnow();    
% $$$     waitforbuttonpress();
% $$$ end
% $$$ 
% $$$ 
% $$$ rmap = pfd.data.rateMap;
% $$$ rmap(~nniz(rmap(:)))=0;
% $$$ [U,S,V] = svd(rmap',0);
% $$$ 
% $$$ figure,
% $$$ for i = 1:10,
% $$$     subplot(2,5,i);
% $$$     imagesc(reshape(V(:,i),pfd.adata.binSizes')');
% $$$     axis('xy');
% $$$     caxis([-0.06,0.06]);
% $$$ end
% $$$ 
% $$$ 
% $$$ [~,V,FSr,VT] = erpPCA(rmap',10);
% $$$ 
% $$$ figure,
% $$$ for i = 1:10,
% $$$     subplot(2,5,i);
% $$$     imagesc(reshape(V(:,i),pfd.adata.binSizes')');
% $$$     axis('xy');
% $$$     caxis([-0.03,0.03]);
% $$$ end
% $$$ 
% $$$ 
% $$$ % SET interp parameters                
% $$$ interpParPfs = struct('bins',{{linspace(-500,500,50),...
% $$$                     linspace(-500,500,50),...
% $$$                     linspace(  -2,  2,50)}},...
% $$$                       'nanMaskThreshold', 0.1,...
% $$$                       'methodNanMap',     'cubic',...
% $$$                       'methodRateMap',    'cubic');
% --------------------------------------------------------------------------------------------------

% END ANALYSIS MAIN --------------------------------------------------------------------------------