


global MTA_PROJECT_PATH
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/fwdrebayesiandecodingtools/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/tempScripts/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/aux/

MjgER2016_load_data();

% $$$ if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
% $$$ numComp = size(eigVec{1},2);
% $$$ pfindex = 1;
% $$$ MjgER2016_load_bhv_erpPCA_scores();
% $$$ % output:
% $$$ %    fsrcz
% $$$ %    FSrC
% $$$ %    rmaps
% $$$ clear('pfd','tags','eigVec','eigVar','eigScore','validDims','zrmMean','zrmStd',...
% $$$       'clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','rmapsShuffledMean',...
% $$$       'rmapsShuffled','FSrM','FSrS','fsrsMean','fsrsStd','fsrsMean','fsrsStd');



for trialIndex =  [1:4,6,7,17,18,20:23];

Trial = Trials{trialIndex}; 
unitSubset = units{trialIndex};
decodingSampleRate = 10;

pfstype = 'xyhb';
%pfstype = 'xy';

% COMPUTE Ndim placefield
if strcmp(pfstype,'xyhb'),
    xyz = preproc_xyz(Trial,'SPINE_SPLINE_HEAD_EQI');
    xyz.resample(16);
    fet = fet_HB_pitchB(Trial,16);
    xyzp = copy(fet);
    xyzp.data = cat(2,sq(xyz(:,'nose',[1,2])),xyzp.data);
    pfsArgs = struct('units',              unitSubset,                           ...
                     'states',             'theta-groom-sit',                    ...
                     'overwrite',          false,                                ...
                     'tag',                '',                                   ...
                     'binDims',            [ 100, 100,0.4,0.4],                  ...
                     'SmoothingWeights',   [0.8,0.8,0.8,0.8],                    ...
                     'type',               'xyhb',                               ...
                     'spkShuffle',         false,                                ...
                     'posShuffle',         false,                                ...
                     'numIter',            1000,                                 ...
                     'xyzp',               xyzp,                                 ...
                     'boundaryLimits',     [-500,500;-500,500;-2,2;-2,2],        ...
                     'bootstrap',          false,                                ...
                     'halfsample',         true                                  ...
                     );
    pfsArgs = struct2varargin(pfsArgs);
    pfs = MTAApfs(Trial,pfsArgs{:});    
    interpParPfsNdim = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50),...
                                       linspace(  -2,  2,50),...
                                       linspace(  -2,  2,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
     
elseif strcmp(pfstype,'xy')
    pfs = pfs_2d_theta(Trial,unitSubset);
    interpParPfsNdim = struct('bins',{{linspace(-500,500,50),...
                                       linspace(-500,500,50)}},...
                              'nanMaskThreshold', 0.1,...
                              'methodNanMap',     'linear',...
                              'methodRateMap',    'linear');
    
end




% SET interpolation parameters


pfsBins = cell([1,4]);
if ~isempty(interpParPfsNdim),
    pfsBins(1:numel(interpParPfsNdim.bins)) = interpParPfsNdim.bins;
else,
    pfsBins = pfs.adata.bins;
end
pfsBinsDims = cellfun(@numel,pfsBins);
pfsBinsDims(pfsBinsDims==0) = 1;
pfsBins(cellfun(@isempty,pfsBins)) = [];

% $$$ figure,
% $$$ for u = 1:numel(pfs.data.clu),
% $$$ clf();
% $$$ srmap = pfs.plot(unitSubset(u),'mean',false,[],false,0.25,false,interpParPfsNdim,@jet,mazeMask);
% $$$ subplot2(8,6,1,[1,2]);
% $$$ plot(pft,unitSubset(u),'mean',true,[],false);
% $$$ title(num2str(u));
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ for i = 0:5,
% $$$ for j = 0:6,    
% $$$     subplot2(8,6,j+2,i+1);
% $$$     imagescnan({pfsBins{1},pfsBins{2},srmap(:,:,i*3+20,j*3+15)'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ end
% $$$ end
% $$$ waitforbuttonpress();
% $$$ end
% $$$ mrt = pft.maxRate(unitSubset);



stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = resample(preproc_xyz(Trial,'trb'),decodingSampleRate);
lfp = Trial.load('lfp',72);
fet = fet_HB_pitchB(Trial,decodingSampleRate);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,lfp,'',unitSubset,0.02,true,'count');
%ufr = ufr.create(Trial,xyz,'gper',unitSubset,0.025,true);


% TEST N-dimensional interpolation and masking of rate map
% $$$ 
% $$$ %srmap = srmap.*mask;
% $$$ figure,
% $$$ rmax = pft.maxRate(unitSubset(u));
% $$$ subplot(121);
% $$$ imagescnan({pfsBins{1},pfsBins{2},sq(srmap(:,:,25,25))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');
% $$$ subplot(122);
% $$$ imagescnan({pfsBins{3},pfsBins{4},sq(srmap(30,20,:,:))'},[0,rmax],[],true,'colorMap',@jet);axis('xy');

    
    
% CREATE spatial mask
width = numel(pfsBins{1});
height = numel(pfsBins{2});
radius = round(numel(pfsBins{1})/2)-find(pfsBins{1}<-440,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
circMask = repmat(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),...
                  [1,1,pfsBinsDims(3:4)]);

if strcmp(pfstype,'xyhb'),
    % CREATE behavioral mask
    % LOAD bhv analysis vars
    ds = load(fullfile(MTA_PROJECT_PATH,'analysis','req20180123_pfd_erpPCA-HBPITCHxBPITCH_v7.mat'));
    if ~exist('pfb','var'),  pfb = MTAApfs(Trial,'tag',['HBPITCHxBPITCH_v7']);  end
    interpParPfb = struct('bins',{{pfsBins{3},...
                                   pfsBins{4}}},...
                          'nanMaskThreshold', 0.01,...
                          'methodNanMap',     'linear',...
                          'methodRateMap',    'linear');
    bhvMask = zeros([cellfun(@numel,pfb.adata.bins)]);
    bhvMask(ds.vDims) = 1;
    interpGrids = cell([1,numel(interpParPfb.bins)]);
    [interpGrids{:}] = ndgrid(interpParPfb.bins{:});
    bhvMask = interpn(pfb.adata.bins{:},bhvMask,interpGrids{:},interpParPfb.methodRateMap);
    bhvMask(isnan(bhvMask)) = 0;
    sinterpGrids = interpGrids;
    SmoothingWeights = [0.1,0.1];
    for i = 1:ndims(sinterpGrids),
        sinterpGrids{i} = sinterpGrids{i}.^2/SmoothingWeights(i)^2/2;
    end
    Smoother = exp(sum(-cat(ndims(sinterpGrids)+1,sinterpGrids{:}),ndims(sinterpGrids)+1));
    Smoother = Smoother./sum(Smoother(:));
    bhvMask   = convn(bhvMask, Smoother,'same');
    bhvMask = double(bhvMask>0.05);
    bhvMask = repmat(permute(bhvMask,[3,4,1,2]),...
                     [pfsBinsDims(1:2),1,1]);
else
    bhvMask = 1;    
end


% CREATE rate map mask
mazeMask = circMask.*bhvMask;
mazeMask(mazeMask==0)=nan;


% ACCUMULATE rate maps for bayesian decoding
rateMap = zeros([cell2mat(cf(@numel,pfsBins)),0]);
for u = 1:numel(unitSubset),
    trm = pfs.plot(unitSubset(u),'mean',false,[],false,0.25,false,interpParPfsNdim,[],mazeMask);
    rateMap = cat(numel(pfsBins)+1,rateMap,trm);
end

% REDUCE prior dimensionality
mapSize = size(rateMap);
erows = {};
grows = {};
eeCount = [];
for i = 1:ndims(rateMap)-1;
    eeCount(:,i) = sum(reshape(permute(isnan(rateMap),[i,find(~ismember(1:ndims(rateMap),i))]),mapSize(i),[]),2);
    erows{i}(:) = eeCount(:,i)==prod(mapSize(~ismember(1:ndims(rateMap),i)));
    grows{i}(:) = eeCount(:,i)~=prod(mapSize(~ismember(1:ndims(rateMap),i)));
end
grows{end+1} = ones([1,numel(unitSubset)]);
erows{end+1} = ones([1,numel(unitSubset)]);

% MANUAL Reduction of a few dims Matrix size is too big
% jg05-20120310
% $$$ grows{1}([27:28]) = 0;
% $$$ erows{1}([27:28]) = 1;
% $$$ grows{2}([3,4,27:28]) = 0;
% $$$ erows{2}([3,4,27:28]) = 1;
% jg05-20120312
% $$$ grows{2}(27:28) = 0;
% $$$ erows{2}(27:28) = 1;
gpfsBins = {};
for i = 1:ndims(rateMap)-1;
    gpfsBins{i} = pfsBins{i}(grows{i});
end

smap = rateMap(grows{:});
smap = rateMap;
for i = 1:ndims(rateMap)-1,
    smap(erows{i},:,:,:,:) = [];
    smap = permute(smap,[2:ndims(rateMap),1]);
end
smap = permute(smap,[2:ndims(rateMap),1]);


% COMPUTE Posterior distribution base on ratemaps and unit firing rates
smap(isnan(smap)) = 0;
atE = [];
clear('tE');

bufferSize = 1000;
i = 1;
%tE = decode_bayesian_poisson(smap,ufr.data(1+(bufferSize*(i-1)):bufferSize*i,:)'-eps);


g = decode_bayesian_poisson(smap,ufr.data(1,:)');
atE = cat(ndims(tE),atE,tE);

% SELECT ripples
rper = stc{'g',ufr.sampleRate};
rufr = ufr(rper,:);
atE =  decode_bayesian_poisson(smap,rufr(1:5000,:)');

% ESTIMATE positions from posterior
tpos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
binsE = {};
for i = 1:numel(pfs.adata.bins),  binsE{i} = pfsBins{i}(grows{i});  end
gbins = cell([1,numel(pfsBins)]); 
[gbins{:}] = ndgrid(binsE{:});
gbins = cat(numel(gbins)+1,gbins{:});
ss = substruct('()',[repmat({':'},[1,ndims(gbins)-1]),{1}]);

for tind = 1:size(tE,ndims(tE)),
    ss.subs{end} = tind;
    tpos(tind,:) = sum(reshape(gbins.*repmat(subsref(tE,ss),[ones([1,ndims(gbins)-1]),...
                        size(gbins,ndims(gbins))]),...
                               [],size(gbins,ndims(gbins))));
end


%[~sort(reshape(atE(:,:,:,:,100),[],1));

figure()
subplot(311);  plot(tpos(:,[1,2]));  %Lines([0;cumsum(diff([rper.data,1,2))],[],'k');
subplot(312);  plot(tpos(:,[3,4]));  %Lines([0;cumsum(diff(rper.data,1,2))],[],'k');
subplot(313);  imagesc(rufr(1:tind,:)');  %Lines([0;cumsum(diff(rper.data,1,2))],[],'k');
linkaxes(findobj(gcf,'Type','Axes'),'x');



% CREATE video of decoding
% $$$ xbinsReal = discretize(xyz(:,'nose',1),binsE{1});
% $$$ ybinsReal = discretize(xyz(:,'nose',2),binsE{2});
% $$$ 
% $$$ xbinsReal = discretize(tpos(:,1),binsE{1});
% $$$ ybinsReal = discretize(tpos(:,2),binsE{2});
% $$$ hbinsReal = discretize(tpos(:,3),binsE{3});
% $$$ bbinsReal = discretize(tpos(:,4),binsE{4});
% $$$ 
% $$$ 
% $$$ [C,H] = bhv_contours();
% $$$ 
% $$$ H{1}.ZData(sqrt(sum(cat(3,H{1}.XData,H{1}.YData).^2,3))<0.2)=0;
% $$$ 
% $$$ 
% $$$ vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/decode/posterior_example_immobile.avi','Uncompressed AVI');
% $$$ open(vidObj);
% $$$ 
% $$$ hfig = figure();
% $$$ hfig.Units = 'centimeters';
% $$$ pause(0.1);
% $$$ hfig.Position(3:4) = [30,20];
% $$$ pause(0.1);
% $$$ axPosition = axes('Units','centimeters','Position',[2,4,6,6]);
% $$$ hold(axPosition,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axPosture = axes('Units','centimeters','Position',[2,12,6,6]);
% $$$ hold(axPosture,'on');
% $$$ daspect(axPosition,[1,1,1]);
% $$$ axBg = axes('Position',[0 0 1 1],'Visible','off');    
% $$$ 
% $$$ axPosX = axes('Units',         'centimeters',...
% $$$               'Position',      [13,16,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('X (cm)');
% $$$ haxLinesTime(1) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,13,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10,'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10,'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-50,50]);
% $$$ ylabel('Y (cm)');
% $$$ haxLinesTime(2) = animatedline([310,310],[-50,50],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,10,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Head','Pitch (rad)'});
% $$$ haxLinesTime(3) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ 
% $$$ axPosY = axes('Units',         'centimeters',...
% $$$               'Position',      [13,7,14,3],...
% $$$               'NextPlot',      'add',...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on',...
% $$$               'XTickLabel',    {});
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4),'m','LineWidth',2);
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2),'g','LineWidth',2);
% $$$ xlim([310,375]);
% $$$ ylim([-1.5,1.5]);
% $$$ ylabel({'Body','Pitch (rad)'});
% $$$ haxLinesTime(4) = animatedline([310,310],[-1.5,1.5],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ axStc  = axes('Units',         'centimeters',...
% $$$               'Position',      [13, 4,14,3],...
% $$$               'XMinorGrid',    'on',...
% $$$               'YMinorGrid',    'on');
% $$$ plotSTC(stc,1,[],fliplr({'rear','loc','pause'}),fliplr('rbk'));
% $$$ axStc.YTickLabel = {'pause','loc','rear'};
% $$$ xlim([310,375]);
% $$$ xlabel('Time (s)');
% $$$ haxLinesTime(5) = animatedline([310,310],[1,4],'Color','k','LineStyle','--','LineWidth',2);
% $$$ 
% $$$ limitsBound = [-50,50;-50,50;-1.5,1.5;-1.5,1.5,;1,4];
% $$$ 
% $$$ for t = 3100:1:3750
% $$$ 
% $$$     axes(axPosition);    
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-50,50,100),linspace(-50,50,100),zeros([100,100]));
% $$$         axPositionCirc = circle(0,0,42,'--r');        
% $$$     end
% $$$ 
% $$$     if ~isnan(hbinsReal(t))&~isnan(bbinsReal(t)),
% $$$         axPositionIm = imagesc(binsE{1}/decodingSampleRate,binsE{2}/decodingSampleRate,sq(tE(:,:,hbinsReal(t),bbinsReal(t),t))');
% $$$         axPositionScXyz = scatter(xyz(t,'nose',1)/decodingSampleRate,xyz(t,'nose',2)/decodingSampleRate,40,'g','filled');
% $$$         axPositionScTpos = scatter(tpos(t,1)/decodingSampleRate,tpos(t,2)/decodingSampleRate,40,'m','filled');
% $$$     end    
% $$$     uistack(axPositionCirc,'top')
% $$$ 
% $$$     if t == 3100,        
% $$$ 
% $$$     end
% $$$     
% $$$     xlabel('X Position (cm)');
% $$$     ylabel('Y Position (cm)');
% $$$     title('Decoded Position');
% $$$     xlim([-50,50]);    
% $$$     ylim([-50,50]);    
% $$$     daspect([1,1,1]);
% $$$ 
% $$$     axes(axPosture);
% $$$     if t == 3100,    
% $$$         imagesc(linspace(-1.65,1.4,100),linspace(-1,1.75,100),zeros([100,100]));
% $$$     end
% $$$     
% $$$     if ~isnan(xbinsReal(t))&~isnan(ybinsReal(t)),
% $$$         axPostureIm = imagesc(binsE{3},binsE{4},sq(tE(xbinsReal(t),ybinsReal(t),:,:,t))');
% $$$         axPostureScFet = scatter(fet(t,1),fet(t,2),40,'g','filled');
% $$$         axPostureScTpos = scatter(tpos(t,3),tpos(t,4),40,'m','filled');
% $$$     end
% $$$     xlabel('Head-Body Pitch (rad)');
% $$$     ylabel('Body Pitch (rad)')
% $$$     title('Decoded Posture');
% $$$     xlim([-1.65,1.4]);
% $$$     ylim([-1,1.75]);
% $$$     
% $$$     if t == 3100,
% $$$         HC = cf(@(h) copyobj(h,gca), H);
% $$$         h = plot(0,0,'w');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'c');
% $$$         h.Visible = 'off';
% $$$         h = plot(0,0,'r');
% $$$         h.Visible = 'off';        
% $$$         lax = legend({'real','decode','low','high','rear'},'Location','eastoutside');
% $$$         lax.Position = lax.Position+[0.1,0,0,0];
% $$$         lax.Color = [0.6,0.6,0.6];
% $$$     end
% $$$     cf(@(h) uistack(h,'top'), HC);
% $$$     
% $$$ 
% $$$     
% $$$     axes(axBg);
% $$$     cla(axBg);
% $$$     text(0.14,0.1,['time: ',num2str(t/decodingSampleRate),' s']);
% $$$     
% $$$     af(@(h) clearpoints(h), haxLinesTime);
% $$$     for h = 1:numel(haxLinesTime),
% $$$          addpoints(haxLinesTime(h),[t/decodingSampleRate,t/decodingSampleRate],limitsBound(h,:));
% $$$     end
% $$$ 
% $$$     drawnow('update');
% $$$ 
% $$$     writeVideo(vidObj,getframe());
% $$$     
% $$$     % Clear axes
% $$$     delete([axPositionIm,axPositionScXyz,axPositionScTpos]);    
% $$$     delete([axPostureIm,axPostureScFet,axPostureScTpos]);
% $$$     
% $$$ end
% $$$ close(vidObj);

% $$$ 
% $$$ % DIAGNOSTIC plot of decoded and real positions in physical and behavioral spaces
% $$$ figure();
% $$$ subplot(511);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,1)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',1)/10);
% $$$ title('x pos');
% $$$ subplot(512);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,2)/10);
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,xyz(1:tind,'nose',2)/10);
% $$$ title('y pos');
% $$$ subplot(513);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,3));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,1));
% $$$ title('head pitch');
% $$$ subplot(514);
% $$$ plot((1:tind)/decodingSampleRate,tpos(1:tind,4));
% $$$ hold('on');
% $$$ plot((1:tind)/decodingSampleRate,fet(1:tind,2));
% $$$ title('body pitch');
% $$$ subplot(515);
% $$$ plotSTC(stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause','theta'}),fliplr('rbcgkm'));
% $$$ xlabel('Time (s)');
% $$$ linkaxes(get(gcf,'Children'),'x');
% $$$ linkaxes(get(gcf,'Children'),'x');


% ENCAPSULATE decoded position in MTAData object

rpos = MTADfet.encapsulate(Trial,                                                             ... % MTATrial object
                           tpos,                                                              ... % Data
                           decodingSampleRate,                                                ... % Sample rate
                           ['bayesiandecoded_',pfs.parameters.type,'_',pfs.parameters.states],... % Name
                           ['bd',pfs.parameters.type],                                        ... % Label
                           'b'                                                                ... % key
);

ind = [stc{'a&t-m-s',xyz.sampleRate}];


% FIGURE - jpdf of physical and behavioral spaces
figure,
subplot(121);
out = hist2(rpos(ind,[1,2])/10,linspace([-50,50,50]),linspace([-50,50,50]));
imagesc(linspace([-50,50,50]),linspace([-50,50,50]),out'./sum(out(:)));
axis('xy');
caxis([0,max(caxis)/2]);
daspect([1,1,1]);
title({'Decoded Position',Trial.filebase});
xlabel('cm');
ylabel('cm');
if strcmp(pfstype,'xyhb');
    subplot(122);
    out = hist2(rpos(ind,[3,4]),linspace(-1.5,1,50),linspace(-1,1.8,50));
    imagesc(linspace(-1.5,1,50),linspace(-1,1.8,50),out'./sum(out(:)));
    axis('xy');
    caxis([0,max(caxis)/2]);
    daspect([1,1,1]);
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_jpdfs_pos_a_bhv_',pfstype,'_',Trial.filebase,'.png']]);



meanPopRate = sum(ufr(ind,:),2);
meanPopIncl = sum(ufr(ind,:)>1,2);

rposError = [sqrt((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2)];
spatialError = sqrt(sum((sq(xyz(ind,'nose',[1,2]))-rpos(ind,[1,2])).^2,2));
if strcmp(pfstype,'xyhb'),
    rposError = cat(2,rposError,sqrt((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2));  
    bhvError = sqrt(sum((sq(fet(ind,[1,2]))-rpos(ind,[3,4])).^2,2));
else
    bhvError = [];
    bhvErrorCondPos = [];
end






% FIGURE - jpdf of phisical and behavioral errors
if strcmp(pfstype,'xyhb'),
figure();
hist2([spatialError/10,bhvError],linspace(0,20,25),linspace(0,1.4,25));
title({'Spatial Vs Behavioral Error',[Trial.filebase,' ',num2str(numel(unitSubset)),' units']});
xlabel('cm');
ylabel('rad');
hax = gca();
hax.Units = 'centimeters';
hax.Position = [hax.Position([1,2]),4,4];
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_jpdf_pos_x_bhv_',pfstype,'_',Trial.filebase,'.png']]);
end


% COMPUTE mean error for decoded components individually
shuffle = @(x) circshift(x,randi(size(x,1)));
ind = [stc{'a&t-m-s',xyz.sampleRate}];
meanPopRate = sum(ufr(ind,:),2);
txss = sq(xyz(ind,'nose',[1,2]));
tfss = fet(ind,[1,2]);
trps = rpos(ind,:);
meanError = mean(rposError);
% COMPUTE shuffled error for decoded components
meanErrorShuffled = nan([1000,numel(pfsBins)]);
bhvErrorShuffled = [];
mbes = [];
for i = 1:1000,
    sptErrorShuffled = sqrt((shuffle(txss)-trps(:,[1,2])).^2);
    if strcmp(pfstype,'xyhb'),
        bhvErrorShuffled = sqrt((shuffle(tfss)-trps(:,[3,4])).^2);
        mbes = mean(bhvErrorShuffled(meanPopRate>0.6,:));
    end
    meanErrorShuffled(i,:) = [mean(sptErrorShuffled(meanPopRate>0.6,:)), mbes];
end

% COMPUTE zscore of
for e =1:numel(meanError),
zsMSError(e) = (meanError(e)-mean(meanErrorShuffled(:,e)))/std(meanErrorShuffled(:,e));
end



% FIGURE - mean error conditiond on position
posBins = linspace(-50,50,25);
xBins = discretize(xyz(ind,'nose',1)./10,posBins);
yBins = discretize(xyz(ind,'nose',2)./10,posBins);
aind = xBins>0&yBins>0;
figure();
subplot(121);
spatialErrorCondPos = accumarray([xBins(aind),yBins(aind)],spatialError(aind)/10,[numel(posBins),numel(posBins)],@median);
spatialErrorCondPos(spatialErrorCondPos==0) = nan;
imagescnan({posBins,posBins,spatialErrorCondPos'},[0,20],'linear',true,'colorMap',@jet);axis('xy');
title({'Mean Decoded Spatial Error','Conditioned on Position'});
if strcmp(pfstype,'xyhb'),
    subplot(122);
    bhvErrorCondPos = accumarray([xBins(aind),yBins(aind)],bhvError(aind),[numel(posBins),numel(posBins)],@median);
    bhvErrorCondPos(bhvErrorCondPos==0) = nan;
    imagescnan({posBins,posBins,bhvErrorCondPos'},[0,1],'linear',true,'colorMap',@jet);axis('xy');
    title({'Mean Decoded Pitch Error','Conditioned on Position'});
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
hax(1).Position = [hax(1).Position([1,2]),4,4];
hax(2).Position = [hax(2).Position([1,2]),0.5,4];
if strcmp(pfstype,'xyhb'),
    hax(3).Position = [hax(3).Position([1,2]),4,4];
    hax(4).Position = [hax(4).Position([1,2]),0.5,4];
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_pos_a_bhv_Cond_pos',pfstype,'_',Trial.filebase,'.png']]);



% FIGURE - Independence of errors
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,14,16];

subplot(423);hist2(mud([rposError(:,[1,2])]),50,50);title('x vs y');
if strcmp(pfstype,'xyhb'),
    subplot(421);hist2(mud([spatialError,bhvError]),50,50);title('xy vs hb');
    subplot(424);hist2(mud([rposError(:,[3,4])]),50,50);title('h vs b');
    subplot(425);hist2(mud([rposError(:,[1,3])]),50,50);title('x vs h');
    subplot(426);hist2(mud([rposError(:,[2,4])]),50,50);title('y vs b');
    subplot(427);hist2(mud([rposError(:,[1,4])]),50,50);title('x vs b');
    subplot(428);hist2(mud([rposError(:,[2,3])]),50,50);title('y vs h');
end
hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]),  hax);
axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
text(0.5,0.9,{Trial.filebase,['units: ',num2str(numel(unitSubset))]});
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_error_independence_',pfstype,'_',Trial.filebase,'.png']]);





save(fullfile(MTA_PROJECT_PATH,'analysis',['bhv_decode_',pfstype,'_',Trial.filebase,'.mat']),...
     'ind',                         ...
     'rposError',                   ...
     'spatialError',                ...
     'bhvError',                    ...
     'spatialErrorCondPos',         ...     
     'bhvErrorCondPos',             ...
     'posBins',                     ...
     'meanError',                   ...
     'meanErrorShuffled',           ...
     'zsMSError',                   ...
     'meanPopRate',                 ...
     'meanPopIncl',                 ...
     'unitSubset'                   ...
);

end


clear('ds');
trialSubset = [1:4,6,7,17,18,20:23];
pfstypes = {'xy','xyhb'};
for trialIndex =  1:numel(trialSubset);
    Trial = Trials{trialSubset(trialIndex)}; 
    for typeIndex = 1:numel(pfstypes)
        ds(trialIndex,typeIndex) = load(fullfile(MTA_PROJECT_PATH,'analysis',...
                              ['bhv_decode_',pfstypes{typeIndex},'_',Trial.filebase,'.mat']));
    end
end

me2d = reshape([ds(:,1).meanError],[2,12])';
me4d = reshape([ds(:,2).meanError],[4,12])';

me2d = [ds(:,1).spatialError']';
me4d = reshape([ds(:,2).meanError],[4,12])';



me2dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,1)));
me4dx = cell2mat(af(@(s) median(s.rposError(:,1)), ds(:,2)));

me2dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,1)));
me4dy = cell2mat(af(@(s) median(s.rposError(:,2)), ds(:,2)));

me2s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,1)));
me4s = cell2mat(af(@(s) median(s.spatialError(:,1)), ds(:,2)));

me4p = cell2mat(af(@(s) median(s.bhvError(:,1)), ds(:,2)));

nunits = cell2mat(af(@(s) numel(s.unitSubset(:)), ds(:,1)));

figure();
subplot
hold('on');
scatter(nunits,me2dx,20,'b','filled')
scatter(nunits,me4dx,20,'r','filled')

figure();
hold('on');
scatter(nunits,me2dy,20,'b','filled')
scatter(nunits,me4dy,20,'r','filled')



figure();
subplot(121);
hold('on');
scatter(nunits,me2s/10,20,'b','filled');
scatter(nunits,me4s/10,20,'r','filled');
xlabel('Number of Units');
ylabel('Distance (cm)');
title('Median Spatial Error');
legend({'2d place fields','4d place fields'},'location','eastoutside');
grid('on')
set(gca,'XTick',[0,25,50,75,100]);
set(gca,'YTick',[4,5,6,7,8,9,10]);

subplot(122);
scatter(nunits,me4p,20,'r','filled')
xlabel('Number of Units');
ylabel('Distance (rad)');
title('Median Behavioral Error');
legend({'4d place fields'},'location','eastoutside');
grid('on');
set(gca,'XTick',[0,25,50,75,100]);

hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),3,3]),  hax);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/decode/',...
       ['dc_regress_nunits_x_error.png']]);




% 20180725 ----------------------------------------------
trialIndex = 20;
Trial = Trials{trialIndex};
stc = Trial.load('stc','msnn_ppsvd_raux');
unitSubset = units{trialIndex};

xyz = preproc_xyz(Trial,'trb');
xyz.resample(250);

vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','head_back'});

states = {'rear','hloc','hpause','lloc','lpause','groom','sit','theta','spw'};

drz = xyz.copy();
drz.data = compute_drz(Trial,unitSubset,pfs,'feature',xyz.copy());

lfp = Trial.load('lfp',[68,72,76,82]);

phz = phase(resample(copy(lfp),xyz),[5,12]);

spk = Trial.load('spk',xyz.sampleRate,[],unitSubset);

ufrAll = Trial.load('ufr',xyz,[],[],0.04);

unitsInt = select_units(Trial,'int');
ufrInt = Trial.load('ufr',xyz,[],unitsInt,0.04);

int = Trial.load('spk',xyz.sampleRate,'',unitsInt);

figure
for unit = 1:numel(unitsInt),
    subplotfit(unit,20);
    rose(phz(int(unitsInt(unit)),1));
end



ufrIntObj = ufrInt.copy();
ufrIntObj.data = mean(ufrInt(:,[5,6,10]),2);
uiphz = ufrIntObj.phase([5,12]);

ufrAllObj = ufrAll.copy();
ufrAllObj.data = mean(ufrAll.data,2);


figure,
plot(mean(ufrInt(:,[5,6,10]),2));
hold('on');
plot(phz(:,2)*10+50);
plot(uiphz(:)*10+50);


ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'',unitSubset,0.01,true,'gauss');

xyn = xyz.copy();
xyn.data = sq(xyz(:,'nose',[1,2]));
xyl = xyn.copy();;
xyl.resample(lfp);

specArgs = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,40]);
   
[ys,fs,ts] = fet_spec(Trial,lfp,[],[],[],specArgs);

figure,
imagesc(ts,fs,log10(ys(:,:,1)));

lys = ys.copy;
lys.data = mean(log10(ys(:,fs<15,2)),2);
lys.resample(xyz);

lds = ys.copy;
lds.data = mean(log10(ys(:,5<fs & fs<12,4)),2)./mean(log10(ys(:,12<fs&fs<18,4)),2);
lds.resample(xyz);

lrs = ys.copy;
lrs.data = mean(log10(ys(:,5<fs&fs<12,3)),2)./mean(log10(ys(:,5<fs&fs<12,4)),2);
lrs.resample(xyz);

uind = sum(ufr.data-eps,2)~=0&nniz(xyz);

bins = ones([1,30000]).*0.004;
tEn = decode_bayesian_poisson(smap,(ufr(1:30000,:)'./70),'bins',bins,'baseline',0.01);

mtEn = max(reshape(tEn,[],30000));


% ESTIMATE positions from posterior
eposn = nan([size(tEn,ndims(tEn)),ndims(tEn)-1]);
binsE = {};
for i = 1:numel(pfs.adata.bins),  binsE{i} = pfsBins{i}(grows{i});  end
gbins = cell([1,numel(pfsBins)]); 
[gbins{:}] = ndgrid(binsE{:});
gbins = cat(numel(gbins)+1,gbins{:});
ss = substruct('()',[repmat({':'},[1,ndims(gbins)-1]),{1}]);

for tind = 1:size(tEn,ndims(tEn)),
    ss.subs{end} = tind;
    eposn(tind,:) = sum(reshape(gbins.*repmat(subsref(tEn,ss),[ones([1,ndims(gbins)-1]),...
                        size(gbins,ndims(gbins))]),...
                               [],size(gbins,ndims(gbins))));
end



figure,
hold('on');

for t = 1:30000,    
    if mtEn(t)<0.002||unitInclusion(t)<3, continue, end
    %if unitInclusion(t)<3, continue, end
    i = imagesc(gpfsBins{:},tEn(:,:,t)');
    p = plot(xyz(t,'hcom',1),xyz(t,'hcom',2),'*m');
    g = plot(eposn(t,1),eposn(t,2),'*g');
    drawnow();
    pause(0.05);
    delete(i);
    delete(p);
    delete(g);
end



tE = decode_bayesian_poisson(smap,(ufr(:,:)'));%,'bins',bins,'baseline',0.001);






% ESTIMATE positions from posterior
epos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
binsE = {};
for i = 1:numel(pfs.adata.bins),  binsE{i} = pfsBins{i}(grows{i});  end
gbins = cell([1,numel(pfsBins)]); 
[gbins{:}] = ndgrid(binsE{:});
gbins = cat(numel(gbins)+1,gbins{:});
ss = substruct('()',[repmat({':'},[1,ndims(gbins)-1]),{1}]);

for tind = 1:size(tE,ndims(tE)),
    ss.subs{end} = tind;
    epos(tind,:) = sum(reshape(gbins.*repmat(subsref(tE,ss),[ones([1,ndims(gbins)-1]),...
                        size(gbins,ndims(gbins))]),...
                               [],size(gbins,ndims(gbins))));
end

% $$$ tE = RectFilter(tE,3,1);
% $$$ tE = permute(RectFilter(permute(tE,[2,1,3]),3,1),[2,1,3]);

tpos = nan([size(tE,ndims(tE)),ndims(tE)-1]);
apos = nan([size(tE,ndims(tE)),1]);
for tind = 1:size(tE,ndims(tE)),
    tbin = LocalMinima2(-tE(:,:,tind),0,30);
    if ~isempty(tbin),
        tpos(tind,:) = [gpfsBins{1}(tbin(1)),gpfsBins{2}(tbin(2))];
        apos(tind)   = tE(tbin(1),tbin(2),tind);
    end
end



stcm = stc2mat(stc,xyz,states);

timeSpec = [1:size(ys,1)]./ys.sampleRate;

cmapLims = [1,2.5;1,2.5;1,3;1,3.5];

figure();

sp = [];
for s = 1:4,
    sp(s+1) = subplot(6,1,s+1);
    imagesc(ts,fs,log10(ys(:,:,s))');
    axis('xy');
    caxis(cmapLims(s,:));
    colormap(gca,'jet');
end
sp(6) = subplot(6,1,6);
haxSTS = plotSTC(Trial.stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');

linkaxes(get(gcf,'Children'),'x');

sp(1) = subplot(611);
hold('on');
position = animatedline();
position.Marker = '*';
position.MarkerSize = 20;
position.MarkerFaceColor  = 'm';
position.MarkerEdgeColor  = 'm';

predicted = animatedline();
predicted.Marker = '*';
predicted.MarkerSize = 20;
predicted.MarkerFaceColor  = 'm';
predicted.MarkerEdgeColor  = 'm';

for i = 5266163:2:5286163,
    if apos(i)<0.001, continue;end
    
    posterior = imagesc(pfsBins{1}(3:end-2),...
                pfsBins{2}(3:end-2),...
                tE(:,:,i)');
    
    position.clearpoints();
    position.addpoints(xyn(i,1),xyn(i,2),10);
    predicted.clearpoints();
    predicted.addpoints(tpos(i,1),tpos(i,2),10);
    daspect([1,1,1]);
    
    t = NearestNeighbour(timeSpec,stper(i));

    xlim(sp(2),[t-5,t+5]);
    drawnow();
    
    delete(posterior);
end


unitInclusion = sum(logical(ufr(:,:)-eps),2);
statePeriods = any(stcm==7|stcm==7,2);
statePeriods = any(stcm==6|stcm==6,2);

%positionError = sqrt(sum((tpos-xyls).^2,2));
positionError = sqrt(sum((epos-sq(xyz(:,'nose',[1,2]))).^2,2));
lowFreqMeanPower =lys(:);
thetaRatioRadLM =lrs(:);




ind = unitInclusion>=1&round(tpos(:,1)/5)~=0&round(tpos(:,2)/5)~=0&statePeriods&nniz(tpos)&apos>0.001;
figure,
subplot(121);
hist2([lowFreqMeanPower(ind),positionError(ind)],linspace(1,2.5,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])
subplot(122);
hist2([thetaRatioRadLM(ind),positionError(ind)],linspace(0.6,1.1,50),linspace(0,800,50));
colormap('jet');
caxis([0,500])



%% segmentation of immobility sequences ------------------------------------------------------------

mufr = unity(mean(ufr.data,2));



hfig = figure();
hfig.CurrentCharacter = 's';
sp = gobjects([1,2]);;
sp(1) = subplot(221);
plot([1:numel(mufr)]./xyz.sampleRate,mufr);
hold('on');
plot([1:numel(mufr)]./xyz.sampleRate,apos);
plot([1:numel(mufr)]./xyz.sampleRate,unity(ufrAllObj.data));
sp(2) = subplot(223);
haxSTS = plotSTC(stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');


sp(3) = subplot(222);
while ~strcmp(hfig.CurrentCharacter,'q')
    if diff(sp(1).XLim).*xyz.sampleRate<10000,
        cla(sp(3));
        segInd = round(sp(1).XLim.*xyz.sampleRate);
        segInd(segInd<1) = 1;
        segInd(segInd>numel(apos)) = numel(apos);    
        segInd = segInd(1):segInd(2);
        segInd = segInd(apos(segInd)>0.1);
        if ~isempty(segInd),
            axes(sp(3));
            imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
            hold('on');
            plot(tpos(segInd,1),tpos(segInd,2),'-g');
            plot(epos(segInd,1),epos(segInd,2),'-c');            
            scatter(tpos(segInd(1),1),tpos(segInd(1),2),10,'m','filled');
            plot(xyz(segInd(1):segInd(end),'nose',1),xyz(segInd(1):segInd(end),'nose',2),'.','MarkerSize',10);
            axis('xy');
        end
        
        drawnow();
        waitforbuttonpress();
    else
        waitforbuttonpress();
        continue    
    end
end


% Pfeiffer, Foster 2013
% Criteria for spw decoding events
% normalized mean Population must > 1 Hz
% no sequence may jump more than 10cm
%
% output vars
% mean observed speed of each segment
% mean speed of segment


tper = [stc{'theta',ufr.sampleRate}];
%sequencePeriods = ThreshCross(RectFilter(apos,11,3),0.2,5);
sequencePeriods = ThreshCross(unitInclusion+apos,3.2,5); 
% REMOVE theta periods
%sequencePeriods(WithinRanges(mean(sequencePeriods,2),tper.data),:) = [];


lastPosition = [0,0];
newSeqPer = [];
seqStart = 1;
seqEnd = 1;
trmSeqPer = sequencePeriods;
ipd = [];
for p = 1:size(sequencePeriods,1),
    segInd = sequencePeriods(p,1):sequencePeriods(p,2);
    segInd = segInd(apos(segInd)>0.1);
    if isempty(segInd),
        trmSeqPer(p,:) = nan;
        continue;
    elseif 100 < sqrt(sum((epos(segInd(1),:)-lastPosition).^2,2)),
        newSeqPer = cat(1,newSeqPer,[sequencePeriods(seqStart,1),sequencePeriods(seqEnd,2)]);
        seqStart = p;        
    end
    seqEnd = p;            
    lastPosition = epos(segInd(end),:);
    ipd(p) =sqrt(sum((epos(segInd(1),:)-lastPosition).^2,2));
end
trmSeqPer(isnan(trmSeqPer(:,1)),:) = [];
newSeqPer(1,:) = [];



ysg.resample(xyz);


seq.sampleRate = xyz.sampleRate;
%seq.periods = newSeqPer;
seq.periods = sequencePeriods;

seq.duration = nan([size(seq.periods,1),1]);
seq.obsDistance = nan([size(seq.periods,1),1]);
seq.decDistance = nan([size(seq.periods,1),1]);
seq.obsDiffMean = nan([size(seq.periods,1),1]);
seq.decDiffMean = nan([size(seq.periods,1),1]);
seq.obsDiffStd = nan([size(seq.periods,1),1]);
seq.decDiffStd = nan([size(seq.periods,1),1]);
seq.obsDiffMax = nan([size(seq.periods,1),1]);
seq.decDiffMax = nan([size(seq.periods,1),1]);
seq.meanError = nan([size(seq.periods,1),1]);

seq.obsPathDistMean = nan([size(seq.periods,1),1]);
seq.decPathDistMean = nan([size(seq.periods,1),1]);
seq.obsPathDistStd = nan([size(seq.periods,1),1]);
seq.decPathDistStd = nan([size(seq.periods,1),1]);
seq.obsPathAngMean = nan([size(seq.periods,1),1]);
seq.decPathAngMean = nan([size(seq.periods,1),1]);
seq.obsPathAngStd = nan([size(seq.periods,1),1]);
seq.decPathAngStd = nan([size(seq.periods,1),1]);

seq.obsPathAngPPC = nan([size(seq.periods,1),1]);
seq.decPathAngPPC = nan([size(seq.periods,1),1]);
seq.errorMean = nan([size(seq.periods,1),1]);
seq.errorMin = nan([size(seq.periods,1),1]);
seq.lysMean = nan([size(seq.periods,1),1]);
seq.lrsMean = nan([size(seq.periods,1),1]);
seq.ldsMean = nan([size(seq.periods,1),1]);
seq.ysgMean = nan([size(seq.periods,1),15,4]);



dvxy = vxy.data;
dxyz = xyz(:,'nose',[1,2]);
depos = permute(epos,[1,3,2]);
for p = 1:size(seq.periods,1),
    segInd = seq.periods(p,1):seq.periods(p,2);
    segInd = segInd(apos(segInd)>0.1);

    if numel(segInd)>5,
        seq.duration(p) = diff(segInd([1,end]));
        seq.obsDistance(p) = sqrt(sum(diff(dxyz(segInd([1,end]),1,:)).^2,3));
        seq.decDistance(p) = sqrt(sum(diff(depos(segInd([1,end]),1,:)).^2,3));

        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);
        mdiffXYZ = repmat(dxyz(segInd,1,:),[1,numel(segInd),1])-repmat(permute(dxyz(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);
        mdistXYZ = sqrt(sum(mdiffXYZ.^2,3));

        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);
        mdiffEpos = repmat(depos(segInd,1,:),[1,numel(segInd),1])-repmat(permute(depos(segInd,1,:),[2,1,3]),[numel(segInd),1,1]);
        mdistEpos = sqrt(sum(mdiffEpos.^2,3));    

        seq.obsPathDistMean(p) = mean(diag(mdistXYZ,1));
        seq.decPathDistMean(p) = mean(diag(mdistEpos,1));

        seq.obsPathDistStd(p) = std(diag(mdistXYZ,1));
        seq.decPathDistStd(p) = std(diag(mdistEpos,1));
        

        obsAng = atan2(diag(mdiffXYZ(:,:,1),1),diag(mdiffXYZ(:,:,2),1));
        decAng = atan2(diag(mdiffEpos(:,:,1),1),diag(mdiffEpos(:,:,2),1));        
        
        seq.obsPathAngMean(p) = circ_mean(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngMean(p) = circ_mean(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        seq.obsPathAngStd(p) = circ_std(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngStd(p) = circ_std(circ_dist(decAng(1:end-1),decAng(2:end)));

        seq.obsPathAngPPC(p) = PPC(circ_dist(obsAng(1:end-1),obsAng(2:end)));
        seq.decPathAngPPC(p) = PPC(circ_dist(decAng(1:end-1),decAng(2:end)));
        
        
        seq.obsDiffMean(p) = mean(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffMean(p) = mean(mdistEpos(mdistEpos(:)~=0));
        
        seq.obsDiffStd(p) = std(mdistXYZ(mdistXYZ(:)~=0  ));
        seq.decDiffStd(p) = std(mdistEpos(mdistEpos(:)~=0));

        if sum(mdistXYZ(:)~=0)~=0,
            seq.obsDiffMax(p) = max(mdistXYZ(mdistXYZ(:)~=0  ));
        else
            seq.obsDiffMax(p) = 0;
        end
        seq.decDiffMax(p) = max(mdistEpos(mdistEpos(:)~=0));
        
        seq.lysMean(p) = mean(lys(segInd));
        seq.lrsMean(p) = mean(lrs(segInd));        
        seq.ldsMean(p) = mean(lds(segInd));                
        
        seq.ysgMean(p,:,:) = mean(log10(ysg(segInd,:,:)));
        
        seq.errorMean(p)   = mean(sqrt(sum([depos(segInd,1,:)-dxyz(segInd,1,:)].^2,3)));
        seq.errorMin(p)   = min(sqrt(sum([depos(segInd,1,:)-dxyz(segInd,1,:)].^2,3)));        
    end
    
end



figure();plot(seq.lysMean,log10(seq.obsDiffMean),'.')
figure();plot(seq.lrsMean,log10(seq.obsDiffMean),'.')


notTheta = ~WithinRanges(mean(sequencePeriods,2),tper.data)&nniz(abs(seq.ysgMean(:,1,3)))&abs(seq.lysMean(:))>1.9;
notTheta = ~WithinRanges(mean(sequencePeriods,2),tper.data)&abs(seq.lysMean(:))>1.9;

%plot(log10(seq.decPathDistMean(notTheta)),seq.ysgMean(notTheta,1,3),'.')
hist2([(seq.decPathAngPPC(notTheta)),abs(seq.ysgMean(notTheta,6,4))],linspace(0.2,1,20),linspace(2.7,4.4,20));
hist2([log10(seq.decDiffStd(notTheta)),abs(seq.ysgMean(notTheta,6,4))],linspace(0.2,3,20),linspace(2.7,4.4,20));
hist2([log10(seq.decPathDistStd(notTheta)),abs(seq.ysgMean(notTheta,1,3))],linspace(0,3,20),linspace(2.7,4.1,20));
hist2([log10(seq.errorMin(notTheta)),abs(seq.lysMean(notTheta))],linspace(1,3,25),linspace(1,2.6,25));
hist2([log10(seq.errorMean(notTheta)),abs(seq.lysMean(notTheta))],linspace(1.5,3,25),linspace(1,2.6,25));

figure();
hold('on');
plot(seq.obsPathAngPPC,log10(seq.obsPathDistMean),'.')
plot(seq.decPathAngPPC,log10(seq.decPathDistMean),'.')


figure();
subplot(121);
sp = plot(seq.obsPathAngPPC,seq.obsPathAngMean,'.');



subplot(122);
plot(seq.decPathAngPPC,seq.decPathAngMean,'.')



hfig = figure();
sp = gobjects([0,1]);
sp(end+1) = subplot(221);
hold('on');
scatter(log10(seq.errorMean),log10(seq.errorMin),3,log10(seq.duration),'filled')
%scatter(seq.decPathAngPPC,log10(seq.decPathDistMean),8,log10(seq.duration),'filled')
sp(end+1) = subplot(222);
hold('on');
scatter(seq.lysMean,log10(seq.obsDiffMean),3,log10(seq.duration),'filled')
%linkaxes(get(gcf,'Children'),'xy');
drawnow();
c = gobjects([1,numel(sp)]);;
xy = [0,0];
while isempty(hfig.CurrentCharacter)||hfig.CurrentCharacter~='q',
    waitforbuttonpress();
% REMOVE old unit marker
    delete(c);     
% GET current subplot index    
    axind = find(arrayfun(@isequal,sp,repmat([gca],size(sp))));
% GET xy position of currsor on left mouse button down within current subplot    
    xy = sp(axind).CurrentPoint(1,1:2);
% GET axes Data
    axData = [sp(axind).Children.XData',sp(axind).Children.YData'];
% FIND closest point to currsor on left mouse button down
    dist =sqrt(sum(bsxfun(@minus,xy,axData).^2,2));
    [~,mind]=min(dist);
% PLOT Posterior
    segInd = seq.periods(mind,1):seq.periods(mind,2);
    segInd = segInd(apos(segInd)>0.1);
    subplot(223);
    imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
    caxis([0,0.75]);
    colormap('hot');
    axis('xy');
    hold('on');    
    plot(epos(segInd,1),epos(segInd,2),'m','LineWidth',2);
    plot(epos(segInd(1),1),epos(segInd(1),2),'c*','MarkerSize',10);    
    plot(xyn(segInd,1), xyn(segInd,2),'g*','MarkerSize',3);
    quiver(xyn(segInd,1), xyn(segInd,2),...
             cos(ang(segInd,'head_back','head_front',1))*100,...
             sin(ang(segInd,'head_back','head_front',1))*100,0,'Color',[1,1,1]);
    
    subplot(224);
    
    pd = prctile(dist,5);
    
    scatter(log10(seq.errorMean(dist<pd)),log10(seq.errorMin(dist<pd)),3,log10(seq.duration(dist<pd)),'filled');

    axes(sp(axind));
% HIGHLIGHT selected unit on current subplot
    c(axind) = circle(axData(mind,1),axData(mind,2),0.25,'g');
    for a = find(~ismember(1:numel(sp),axind))
        axes(sp(a));        
        c(a) = circle(sp(a).Children.XData(mind),sp(a).Children.YData(mind),0.25,'k');
    end
end




figure();
plot(seq.obsPathAngPPC,log10(seq.decDiff),'.')
plot(seq.decPathAngPPC,log10(seq.decDiffMax),'.')

figure,
hist2(log10(abs([seq.decMeanDiff,seq.obsDistance])+1e-4),100,100)
Lines([],log10(200),'m');

% Sequence analysis 

% Types of sequences
% THETA 
% PROSPECTIVE 
% SPW 
% 
% Theta Sequence: 
% The hippocampus represents the environment as a squence of action potentials of place cells arranged
% temoporally by the spatial order of the placefields with respect to the rats trajectory. Each sequence
% can last between 
% 
% Filter the sequences based on continuity 



slfp = filter(resample(Trial.load('lfp',[65:96]),xyz),'ButFilter',4,[1,35],'bandpass');
rhm = resample(fet_rhm(Trial),xyz);


stc.states(stc.gsi('events')) = [];
stc.addState([],[],seq.periods,ufr.sampleRate,ufr.sync,ufr.origin,'events','d');

hfig = figure();
sp = gobjects([0,1]);
sp(end+1) = subplot2(6,4,1,[1:3]);
imagesc([1:size(vxy,1)]'./vxy.sampleRate,1:size(ufr,2),ufr.data');
% $$$ for s = 2:4,
% $$$     sp(s) = subplot2(6,4,s,[1:3]);
% $$$     imagesc(ts,fs,log10(ys(:,:,s))');    
% $$$     axis('xy');
% $$$     %caxis(cmapLimsGamma(s,:));
% $$$     colormap(gca,'jet');
% $$$     hold('on');
% $$$     plot([1:size(lfp,1)]./lfp.sampleRate,nunity(lfp.data(:,s))*2+10);    
% $$$ end
sp(end+1) = subplot2(6,4,2:4,[1:3]);
hold('on');
imagesc([1:size(slfp,1)]./slfp.sampleRate,1:32,slfp.data');
plot([1:size(rhm,1)]'./rhm.sampleRate,nunity(rhm.data),'k');
for chan = [7,12,20,29];
    plot([1:size(slfp,1)]'./slfp.sampleRate,-nunity(slfp.data(:,chan))+chan,'m');
end
colormap(gca,'jet');
caxis([-20000,20000]);
axis('ij');

sp(end+1) = subplot2(6,4,5,[1:3]);cla();
hold('on');
plot([[1:size(vxy,1)]'./vxy.sampleRate,[1:size(vxy,1)]'./vxy.sampleRate],vxy.data)
% $$$ for s = 1:numel(sunits)
% $$$     plot([1:size(vxy,1)]'./xyz.sampleRate,ones([numel(spk(sunits(s))),1])./s,'*m');
% $$$ end
plot([1:size(vxy,1)]'./vxy.sampleRate,apos,'-m');
plot([1:size(vxy,1)]'./vxy.sampleRate,lds.data,'-b');
plot([1:size(vxy,1)]'./vxy.sampleRate,lrs.data,'-r');
plot([1:size(vxy,1)]'./vxy.sampleRate,lys.data,'-m');
legend({'body speed','head speed','max posterior','lm theta-delta','rad/lm ratio','pyr <20hz mean PSD'})
sp(end+1) = subplot2(6,4,6,[1:3]);
haxSTS = plotSTC(stc,1,[],[{'events'},fliplr(states)],[,'k',fliplr('rbcgkmy')]);
haxSTS.YTickLabel = [{'events'},fliplr(states)];
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');


c = gobjects([1,numel(sp)]);;
xy = [0,0];
axData = mean([stc{'events',1}.data],2);
while isempty(hfig.CurrentCharacter)||hfig.CurrentCharacter~='q',
    waitforbuttonpress();
% GET current subplot index    
    axind = find(arrayfun(@isequal,sp,repmat([gca],size(sp))));
    if axind == numel(sp),    
% GET time of currsor on left mouse button down within current subplot    
        t = sp(axind).CurrentPoint(1,1);
% FIND closest point to currsor on left mouse button down
        dist =sqrt(sum(bsxfun(@minus,t,axData).^2,2));
        [~,mind]=min(dist);
% PLOT Posterior
        segInd = seq.periods(mind,1):seq.periods(mind,2);
        segInd = segInd(unitInclusion(segInd)>1);
        if numel(segInd)>5,
            subplot2(6,4,[1,2],4);
            imagesc(gpfsBins{:},sum(tE(:,:,segInd),3)');
            caxis([0,0.75]);
            colormap(gca,'hot');
            axis('xy');
            hold('on');    
            plot(epos(segInd,1),epos(segInd,2),'m','LineWidth',2);
            plot(epos(segInd(1),1),epos(segInd(1),2),'c*','MarkerSize',10);    
            plot(xyn(segInd,1), xyn(segInd,2),'g*','MarkerSize',3);
            quiver(xyn(segInd,1), xyn(segInd,2),...
                   cos(ang(segInd,'head_back','head_front',1))*100,...
                   sin(ang(segInd,'head_back','head_front',1))*100,0,'Color',[1,1,1]);
        end
        axes(sp(axind));
    end    
end
