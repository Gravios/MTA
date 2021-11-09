
% LOAD Session data --------------------------------------------------------------------------------

configure_default_args();

MjgER2016_load_data();

% LOAD behavioral data
% $$$ [pfd ,tags ,eigVec, eigVar, eigScore, validDims,...
% $$$  unitSubsets, unitIntersection, zrmMean, zrmStd] = req20180123_ver5(Trials);
% $$$ numComp = size(eigVec{1},2);
% $$$ pfindex = 1;

% LOAD behavioral scores
% $$$ MjgER2016_load_bhv_erpPCA_scores();

trialId = 20;

% RASTER STUFF
Trial = Trials{trialId};    % jg05-20120312.cof.all
unitsSubset = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsPyr = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsInt = Trial.spk.get_unit_set(Trial,'interneurons');

pft = pfs_2d_theta(Trial,unitsPyr);
[mxr,mxp]= pft.get_max_rate(unitsPyr);
mxd = sqrt(sum(mxp.^2,2));
mxn = mxp./mxd;
mxa = atan2(mxn(:,2),mxn(:,1));
mxg = nan(size(mxd));
mxg(mxd<=175) = 1;
abins = linspace(-pi,pi,9);
for a = 1:numel(abins)-1
    mxg(mxd>175& mxa>abins(a) & mxa<=abins(a+1)) = a+1;
end    
[mxv,mxi] = sort(mxg,'descend');
unitsPyr = unitsPyr(mxi);
ucolors = hsv(8);
unitsPyrColor = zeros([size(mxg,1),3]);
unitsPyrColor(mxv>1,:) = ucolors(mxv(mxv>1)-1,:);
unitsIntColor = repmat([0,0,0],[numel(unitsInt),1]);

Trial.load('nq');
% $$$ unitsNon = 1:size(Trial.nq.eDist,1);
% $$$ unitsNon = unitsNon(~ismember(unitsNon,[unitsPyr,unitsInt]));


sampleRate = 250; % Hertz
spikeWindow = 0.025; % Seconds
halfSpkWindow = 0.012; % Seconds

xyz = preproc_xyz(Trial,'trb',sampleRate);

fxyz = filter(xyz.copy(),'ButFilter',3,15,'low')

vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();

ang = create(MTADang,Trial,fxyz);

pch = fet_HB_pitchB(Trial,sampleRate);

bfs      = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);

% COMPUTE bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units, [], [], false);
numComp = size(eigVecs,2);


channels = [66,69,75,82,85];

Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',channels);

phz = copy(lfp);
phz.data = lfp.data(:,1);
phz = phase(phz);
phz.data = phz.data+phzCorrection(trialId);
phz.data(phz.data<-pi) = phz.data(phz.data<-pi)+2*pi;


% $$$ tbp = filter(lfp.copy,'ButFilter',4,[6,12],'bandpass');
% $$$ mrl = zeros([size(lfp,1),1]);
% $$$ mrl(nniz(lfp),:) = abs(hilbert(tbp.data(nniz(lfp),4)));


% $$$ pft = pfs_2d_theta(Trial,1:111,false,true,struct('numIter',1,'halfsample','false'));

% $$$ figure,
% $$$ for u = 1:111,
% $$$     plot(pft,u,1,'text',[],true);
% $$$     title(num2str(u));
% $$$     waitforbuttonpress();
% $$$ end


spk = Trial.load('spk',Trial.lfp.sampleRate,'',[unitsPyr,unitsInt]);
pyr = Trial.load('spk',Trial.lfp.sampleRate,'',unitsPyr);
int = Trial.load('spk',Trial.lfp.sampleRate,'theta',unitsInt);
% $$$ non = Trial.load('spk',Trial.lfp.sampleRate,'',unitsNon);

% Decoding stuff
ufr = Trial.load('ufr',xyz,[],unitsPyr,spikeWindow,'boxcar');
pfs = pfs_2d_states(Trial,unitsPyr,[],{'loc&theta'});
ratemaps = [];
for unit = unitsPyr
    ratemaps(:,end+1) = pfs{1}.data.rateMap(:,pfs{1}.data.clu==unit);
end
mask = create_tensor_mask(pfs{1}.adata.bins);
ratemaps(~mask(:),:) ==nan;

bhvRatemaps = [];
for unit = unitsPyr
    bhvRatemaps(:,end+1) = bfs{trialId}.data.rateMap(:,bfs{trialId}.data.clu==unit);
end
bhvRatemaps(~validDims,:) = nan;


intPhzPref = zeros([numel(unitsInt),1]);
for ii = 1:numel(unitsInt),
    intPhzPref(ii) = circ_mean(phz(int(unitsInt(ii))));
end
intPhzPref(intPhzPref<0) = intPhzPref(intPhzPref<0)+2*pi;
[mpv,mpi] = sort(intPhzPref,'descend');
unitsInt = unitsInt(mpi);





% $$$ 
% $$$ ufrInt = Trial.load('ufr',lfp,spk,unitsInt,0.075,'boxcar');
% $$$ ufrPyr = Trial.load('ufr',lfp,spk,unitsPyr,0.075,'boxcar');
% $$$ 
% $$$ ufrPyrCnt = copy(ufrPyr);
% $$$ ufrPyrCnt.data = sum(double(logical(ufrPyrCnt(:,:)>0.0001)),2);
% $$$ ufrPyrCnt.filter('ButFilter',4,10,'low');
% $$$ 
% $$$ plot([1:size(ufrPyrCnt,1)]./lfp.sampleRate,ufrPyrCnt(:,1))
% $$$ 
% $$$ plot([1:size(ufrPyrFilt,1)]./lfp.sampleRate,sum(double(logical(ufrPyrFilt(:,:)>0.0001)),2))
% $$$ plot([1:size(ufrPyr,1)]./lfp.sampleRate,sum(double(logical(ufrPyr(:,:)>0.0001)),2))


% COMPUTE lfp spectra
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,25]);
   
[lys,lfs,lts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);

thp = mean(lys(:,lfs>5&lfs<10,4),2)./mean(lys(:,lfs<3|(lfs<13&lfs>10),4),2);
figure,plot(thp)
thp(nniz(thp)) = ButFilter(thp(nniz(thp)),4,[0.5]./(lys.sampleRate./2),'low');
hold('on');plot(thp)

% COMPUTE lfp spectra
% $$$ glfp = lfp.copy();
% $$$ glfp.data = glfp.data(:,2);
% $$$ specArgsTheta = struct('nFFT',2^8,...
% $$$                   'Fs',  lfp.sampleRate,...
% $$$                   'WinLength',2^7,...
% $$$                   'nOverlap',2^7*0.875,...
% $$$                   'NW',3,...                  'Detrend',[],...
% $$$                   'nTapers',[],...
% $$$                   'FreqRange',[30,200]);
% $$$    
% $$$ [gys,gfs,gts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);

% $$$ 
% $$$ figure,
% $$$ for c = 1:5
% $$$ subplot(5,1,c)
% $$$ imagesc(gts,gfs,log10(gys(:,:,c))')
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ caxis([2.25,4]);
% $$$ end
% $$$ ForAllSubplots('xlim([4250,4250+5]);')

cmapLimsTheta = [0.8,2.6;0.8,2.6;1,3.2;1,3.5;1,3.5];
states = {'theta','sit','groom','lpause','lloc','hpause','hloc','rear'};
stateColors = 'kymbbggr';

%% START FIG

[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'1080p','landscape',[],1100,50,10,10);
globalXOffset = 0;
globalYOffset = 0;


    
iax = gobjects([0,1]);
vline = gobjects([0,1]);
timeWindow = 12;
step = 0.008;
markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_left','head_front','head_right'};

for ii = 1:numel(xyz.model.Connections)
    xyz.model.Connections{ii}.color = [0,0,1];
    fxyz.model.Connections{ii}.color = [0,0,1];
end

%%%<<< PLOT unit raster timeseries
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(8, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height*10],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
spkRasterHax = sax(end);
unitSet = [fliplr(unitsInt),unitsPyr];
unitClr = [unitsIntColor;unitsPyrColor];
sWidth = 0.4./lfp.sampleRate;
xlim(sax(end),[0,lts(end)]);
ylim(sax(end),[-2*pi,numel(unitSet)+1]);
sphz = copy(phz);
sphz.data(sphz.data<0) = sphz.data(sphz.data<0)+2*pi;
plot([1:size(phz,1)]./phz.sampleRate,sphz(:)-2*pi);
Lines([],1:numel(unitSet),'w');
spkRasterHax.Color = [1,1,1];
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());
[yind, yOffSet, xind, xOffSet] = deal(8, 1, 2, 0);

[yind, yOffSet, xind, xOffSet] = deal(8, 1, 2, 0);
% CREATE subplot axes
thetaHax = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              (fig.page.ypos(yind)+yOffSet+globalYOffset)+fig.subplot.height*10*((2*pi+1)/(numel(unitsInt)+numel(unitsPyr)+2*pi)),...
                              20,                        ...
                              fig.subplot.height*10.*(numel(unitsInt)/(numel(unitsInt)+numel(unitsPyr)+2*pi))],...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
plot(cos(linspace(-pi,pi,100)),linspace(-pi,pi,100));
thetaHax.Visible = 'off';
%thetaHax.Visible = 'on';
thetaHax.YAxisLocation = 'right';
axis(thetaHax,'tight');

% Connect Interneuron raster to theta phase

delete(findobj(fax,'Type','line'))
for intInd = 1:numel(unitsInt)
    line(fax,                                                                        ... axes
         [sum(spkRasterHax.Position([1,3])),                                         ... x1 right edge of spike raster
          sum(spkRasterHax.Position([1,3])) +                                        ... x2 right edge of spike raster
          fig.subplot.horizontalPadding +                                            ... x2 horizontalPadding
          abs(((cos(mpv(intInd)-pi)+1)./2).*thetaHax.Position(3))],                  ... x2 cos(phz) to pixels
         [(fig.page.ypos(yind)+yOffSet+globalYOffset) +                              ... y1 right edge of spike raster
          fig.subplot.height*10*((numel(unitsInt)+1-intInd+2*pi)/(numel(unitsInt)+numel(unitsPyr)+2*pi)),... y1 bot offset from spike raster
          thetaHax.Position(2) +                                                     ... y2 bottom edge of thetaHax
          thetaHax.Position(4)*(mpv(intInd)./diff(ylim(thetaHax)))],         ... y2 offset from bottom spike raster
         'Color','k');
end

%%%<<< PLOT states
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(10, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','Pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height*2],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plotSTC(Trial.stc,1,'text',states,stateColors);
ylim(sax(end),[1,9]);
sax(end).XTickLabel = {};
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/7+1,'r');
plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/7+1,'g');
box(sax(end),'on');
sax(end).Color = [0.9,0.9,0.9];
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT ori layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(11, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(1,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'ori','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT pyramidal layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(12, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(2,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'pyr','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT RAD layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(13, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(3,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'rad','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT LM layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(14, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(4,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
sax(end).XTickLabel = {};
box(sax(end),'on');
ylabel({'lm','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());

%%%<<< PLOT LM layer timeseries of low freq spectra
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(15, 1, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width,                        ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
iax{end+1} = imagesc();
axis(sax(end),'xy');
caxis(sax(end),cmapLimsTheta(4,:));
colormap(sax(end),'jet');
axis(sax(end),'tight');
xlabel(sax(end),'Seconds');
box(sax(end),'on');
ylabel({'dg','Hz'});
ylim(lfs([1,end]));
vline{end+1} = line(sax(end),[timeWindow/2].*[1,1],ylim());
%% LINK spectra and raster
linkaxes(sax,'x');
%%%>>>




%%%<<< PLOT 3D stick-skelleton behaving rat
[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 80);
%[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 40);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
skeletonHax = sax(end);
grid(skeletonHax,'on');
%%%>>>

%%%<<< PLOT behaving rat maze view
%[yind, yOffSet, xind, xOffSet] = deal(11, 20-fig.subplot.height.*5, 2, 150+ fig.subplot.width./5);
[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
%[yind, yOffSet, xind, xOffSet] = deal(9, 20+fig.subplot.height.*5, 2, 150+ fig.subplot.width./5);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeHax = sax(end);
grid(mazeHax,'on');
% inner circle
patch(mazeHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
paxOH = pcolor(mazeHax,                                            ...
               pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
               pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
               reshape(posterior,pfs{1}.adata.binSizes')');
paxOH.EdgeColor = 'none';
caxis(mazeHax,[1.0e-5,0.2]);
axis(mazeHax,'xy');
daspect(mazeHax,[1,1,1]);
%%%>>>


%%%<<< PLOT behaving rate with bayes decoded position (300ms)
[yind, yOffSet, xind, xOffSet] = deal(15, 60+fig.subplot.height.*5, 2, 80);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeBayesHax = sax(end);
% inner circle
patch(mazeBayesHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeBayesHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
pax = pcolor(pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
             pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
             reshape(posterior,pfs{1}.adata.binSizes')');
pax.EdgeColor = 'none';
caxis(mazeBayesHax,[1.0e-5,0.2]);
axis(mazeBayesHax,'xy');
xlim(sax(end),[-500,500]);
ylim(sax(end),[-500,500]);
sax(end).XTick = [-400:200:400];
sax(end).YTick = [-400:200:400];
sax(end).XTickLabel = [];
sax(end).YTickLabel = [];
hold(sax(end),'on');
daspect(sax(end),[1,1,1]);
title(mazeBayesHax,'Bayesian Decoding: Window 300ms');
%%%>>>


%%%<<< PLOT behaving rat with bayes decoded position (40ms)
[yind, yOffSet, xind, xOffSet] = deal(15, 60+fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
% CREATE subplot axes
sax(end+1) = axes('Units','pixels',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                        ...
                              fig.subplot.height.*5],                    ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
mazeBayesFineHax = sax(end);
grid(sax(end),'on');
% inner circle
patch(mazeBayesFineHax,                   ... a
      175.*cos(linspace(-pi,pi,100)), ... x
      175.*sin(linspace(-pi,pi,100)), ... y
      -1.*ones([1,100]),              ... z
      'FaceColor','k',                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
abins = linspace(-pi,pi,9);
% outer sectors
for p = 1:8
patch(mazeBayesFineHax,                   ... a
      [175.*cos(linspace(abins(p),abins(p+1),100)),fliplr(450.*cos(linspace(abins(p),abins(p+1),100)))], ... x
      [175.*sin(linspace(abins(p),abins(p+1),100)),fliplr(450.*sin(linspace(abins(p),abins(p+1),100)))], ... y
      -1.*ones([1,200]),              ... z
      'FaceColor',ucolors(p,:),                ... c
      'EdgeColor','none',             ...
      'FaceAlpha',0.3);
end
halfSpkWindow = 0.04; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1000,sampleRate,ufr,ratemaps,halfSpkWindow);
posterior(~mask(:)) = nan;
posterior(posterior<0.001) = nan;
paxFine = pcolor(mazeBayesFineHax,...
                 pfs{1}.adata.bins{1}-pfs{1}.parameters.binDims(1)/2,...
             pfs{1}.adata.bins{2}-pfs{1}.parameters.binDims(2)/2,...
             reshape(posterior,pfs{1}.adata.binSizes')');
paxFine.EdgeColor = 'none';
caxis(mazeBayesFineHax,[1.0e-5,0.2]);
axis(mazeBayesFineHax,'xy');
xlim(sax(end),[-500,500]);
ylim(sax(end),[-500,500]);
sax(end).XTick = [-400:200:400];
sax(end).YTick = [-400:200:400];
sax(end).XTickLabel = [];
sax(end).YTickLabel = [];
hold(sax(end),'on');
daspect(sax(end),[1,1,1]);
title(mazeBayesFineHax,'Bayesian Decoding: Window 40ms');
%%%>>>


%%%<<< PLOT behavior decoding 300ms
[yind, yOffSet, xind, xOffSet] = deal(15, 280-fig.subplot.height.*5, 2, 80);
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                     ...
                              fig.subplot.height.*5],                   ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
bhvHax = sax(end);
grid(bhvHax,'on');
halfSpkWindow = 0.300; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1500,sampleRate,ufr,bhvRatemaps,halfSpkWindow);
posterior(~validDims(:)) = nan;
posterior(posterior<0.001) = nan;
bhvMask = double(validDims);
pchGrid = cell([2,1]);
[pchGrid{:}] = ndgrid(bfs{1}.adata.bins{:});
paxBhvBkgr = contour(bhvHax,pchGrid{:},reshape(bhvMask,bfs{1}.adata.binSizes'),[0.5,0.5],'-k');
paxBhv = pcolor(bhvHax,                                             ... Axes handle
                bfs{1}.adata.bins{1}-bfs{1}.parameters.binDims(1)/2,... x bins
                bfs{1}.adata.bins{2}-bfs{1}.parameters.binDims(2)/2,... y bins
                reshape(posterior,bfs{1}.adata.binSizes')');          % bhv posterior
paxBhv.EdgeColor = 'none';
xlim([-1.8,0.5]);
ylim([-0.5,1.8]);
pitchMarker = scatter3(bhvHax,pch(1500,1),pch(1500,2),1,20,'r','Filled');
xlabel(bhvHax,'Head-Body Pitch');
ylabel(bhvHax,'Body Pitch');
daspect(bhvHax,[1,1,1]);
%%%>>>


%%%<<< PLOT behavior decoding 40ms
[yind, yOffSet, xind, xOffSet] = deal(15, 280-fig.subplot.height.*5, 2, 190+ fig.subplot.width./5);
sax(end+1) = axes('Units','pixels',                                     ...
                  'Position',[fig.page.xpos(xind)+xOffSet+globalXOffset,...
                              fig.page.ypos(yind)+yOffSet+globalYOffset,...
                              fig.subplot.width./5,                     ...
                              fig.subplot.height.*5],                   ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
bhvFineHax = sax(end);
grid(bhvFineHax,'on');
halfSpkWindow = 0.040; % Seconds
[posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,1500,sampleRate,ufr,bhvRatemaps,halfSpkWindow);
posterior(~validDims(:)) = nan;
posterior(posterior<0.001) = nan;
bhvMask = double(validDims);
pchGrid = cell([2,1]);
[pchGrid{:}] = ndgrid(bfs{1}.adata.bins{:});
paxBhvFineBkgr = contour(bhvFineHax,pchGrid{:},reshape(bhvMask,bfs{1}.adata.binSizes'),[0.5,0.5],'-k');
paxFineBhv = pcolor(bhvFineHax,                                             ... Axes handle
                bfs{1}.adata.bins{1}-bfs{1}.parameters.binDims(1)/2,... x bins
                bfs{1}.adata.bins{2}-bfs{1}.parameters.binDims(2)/2,... y bins
                reshape(posterior,bfs{1}.adata.binSizes')');          % bhv posterior
paxFineBhv.EdgeColor = 'none';
xlim([-1.8,0.5]);
ylim([-0.5,1.8]);
xlabel(bhvFineHax,'Head-Body Pitch');
ylabel(bhvFineHax,'Body Pitch');
daspect(bhvFineHax,[1,1,1]);
pitchFineMarker = scatter3(bhvFineHax,pch(1500,1),pch(1500,2),1,20,'r','Filled');
%%%>>>


% RECORD video while incrementing time
% $$$ vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/explore_spk_lfp_state','Archival');
vidPicDir = create_directory('/storage/share/Projects/BehaviorPlaceCode/explore_spk_lfp_state_rest');
% $$$ vidObj.Height = fig.page.height;
% $$$ vidObj.Width = fig.page.width;

% $$$ open(vidObj);
start = 4200;
x = start;
% $$$ x = 4304-timeWindow-5;
% $$$ x = 4230;
fstart = x;
fend = x+timeWindow;
xlim(spkRasterHax,[fstart,fend]);
lh = gobjects([0,1]);
for u = unitSet
    uind = find(u==unitSet);
    res = spk(u);
    res = res(res./spk.sampleRate>(fstart)&res./spk.sampleRate<fend);    
    %res = res(res./spk.sampleRate<timeWindow);
    if ~isempty(res)
        for r = 1:numel(res)
            lh(end+1) = line(spkRasterHax,repmat(res(r),[1,2])./spk.sampleRate,...
                             [uind,uind+1],'Color',unitClr(uind,:),'LineWidth',1);
        end
    end
end

figFrame = getframe(hfig);
figFrame = repmat(figFrame,[50,1]);

segment = 1;
frame = 1;

% $$$ for x = 0.1:step:10%round(lts(end))  
for x = start:step:start+60%round(lts(end))      
    tic
    fstart = x;
    fend = x+timeWindow;
    xlim(spkRasterHax,[fstart,fend]);
    oldLH = arrayfun(@(l) l.XData(1)<fstart,lh);
    delete(lh(oldLH));
    lh(oldLH) = [];
    for u = unitSet
        uind = find(u==unitSet);
        res = spk(u);
        res = res(res./spk.sampleRate>(fend-step) & res./spk.sampleRate<fend);
        if ~isempty(res)
            for r = 1:numel(res)
            lh(end+1) = line(spkRasterHax,repmat(res(r),[1,2])./spk.sampleRate,...
                 [uind,uind+1],'Color',unitClr(uind,:),'LineWidth',1);
            end
        end
    end

    tindex = round((x+timeWindow/2).*xyz.sampleRate);
    
% UPDATE side perspective
    [mkrs,stks,mkrc] = plotSkeletonLine(Trial,fxyz,tindex,'line',ang,markers,[],[],skeletonHax);
    view(skeletonHax,circ_rad2ang(ang(tindex,'spine_lower','spine_upper',1)+pi/4),20);
    [xp,yp] = deal(mkrs{4}.XData(1),mkrs{4}.YData(1));
    xlim(skeletonHax,xp(1)+[-180,180]);
    ylim(skeletonHax,yp(1)+[-180,180]);
    zlim(skeletonHax,[0,300]);      
    
% UPDATE top perspective fixed
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.3);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxOH.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        paxOH.ZData = ones(size(paxOH.ZData));        
    else
        paxOH.CData = nan(pfs{1}.adata.binSizes');
    end    
    plotSkeletonLine(Trial,xyz,tindex,'line',ang,markers,[],[],mazeHax);
    xlim(mazeHax,xp(1)+[-150,150]);
    ylim(mazeHax,yp(1)+[-150,150]);
    view(0,90);
    
    
% UPDATE spectra
    trange = xlim(iax{1}.Parent);
    trange = round(trange(1).*lys.sampleRate):round(trange(2).*lys.sampleRate);
    for cc = 1:5,
        iax{cc}.XData = lts(trange);
        iax{cc}.YData = lfs;
        iax{cc}.CData = log10(lys(trange,:,cc))';
    end

% UPDATE position decoding overhead
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.3);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeBayesHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        pax.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        pax.ZData = ones(size(pax.ZData));        
        pax.Visible = 'on';        
    else
        pax.CData = nan(pfs{1}.adata.binSizes');
        pax.Visible = 'off';
    end    
    plotSkeletonLine(Trial,fxyz,index,'line',ang,markers,[],[],mazeBayesHax);

% UPDATE position decoding overhead
    index = round((x+timeWindow/2)*ufr.sampleRate);
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,0.04);
    posterior(~mask(:)) = nan;
    posterior(posterior<prctile(posterior,95)) = nan;
    caxis(mazeBayesFineHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxFine.CData = reshape(posterior,pfs{1}.adata.binSizes')';
        paxFine.ZData = ones(size(paxFine.ZData));        
        paxFine.Visible = 'on';
    else
        paxFine.CData = nan(pfs{1}.adata.binSizes');
        paxFine.Visible = 'off';        
    end    
    plotSkeletonLine(Trial,fxyz,index,'line',ang,markers,[],[],mazeBayesFineHax);

% UPDATE Behavior Decoding
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,bhvRatemaps,0.3);
    posterior(~validDims(:)) = nan;
    caxis(bhvHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxBhv.CData = reshape(posterior,bfs{1}.adata.binSizes')';
        paxBhv.ZData = ones(size(paxBhv.ZData));        
        paxBhv.Visible = 'on';
    else
        paxBhv.CData = nan(pfs{1}.adata.binSizes');
        paxBhv.Visible = 'off';
    end    
    pitchMarker.XData = pch(index,1);
    pitchMarker.YData = pch(index,2);    


% UPDATE Fine Behavior Decoding
    [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,bhvRatemaps,0.04);
    posterior(~validDims(:)) = nan;
    caxis(bhvFineHax,[min(posterior),max(posterior)]);
    if ~isempty(posterior)
        paxFineBhv.CData = reshape(posterior,bfs{1}.adata.binSizes')';
        paxFineBhv.ZData = ones(size(paxFineBhv.ZData));        
        paxFineBhv.Visible = 'on';
    else
        paxFineBhv.CData = nan(pfs{1}.adata.binSizes');
        paxFineBhv.Visible = 'off';        
    end    
    pitchFineMarker.XData = pch(index,1);
    pitchFineMarker.YData = pch(index,2);    
    
% SET vertical line time position
    set([vline{:}],'XData',[1,1].*(x+timeWindow/2));
    
    drawnow();

    figFrame(frame) = getframe(hfig);
    %imwrite(frame2im(getframe(hfig)),fullfile(vidPicDir,['explore_spk_lfp_',num2str(frame,'%08d'),'.tiff']),'tiff');    

    delete(findobj(skeletonHax,'Type','line'));
    delete(findobj(mazeHax,'Type','line'));    
    delete(findobj(mazeBayesHax,'Type','line'));
    delete(findobj(mazeBayesFineHax,'Type','line'));    
    
    disp(num2str([frame,toc],'[INFO] frame %d ... %d'))
    frame = frame+1;

    if frame==51
        tic
        cd(vidPicDir);
        for f = 1:50,
            imwrite(frame2im(figFrame(f)),fullfile(vidPicDir,['explore_spk_lfp_',num2str(f,'%08d'),'.tiff']),'tiff');
        end

        disp(num2str([toc],'[INFO] writing images  ... %d'));
        tic
        system(['ffmpeg -framerate 30 -i ''explore_spk_lfp_%08d.tiff''  rest_jg05-20120312_seg',num2str(segment,'%08d'),'.mp4 && rm ./*.tiff']);
        disp(num2str([toc],'[INFO] writing video segment  ... %d'));
        segment = segment+1;
        frame = 1;
    end
    
end

system('ffmpeg -safe 0 -f concat -i <(find . -type f -name ''*'' -printf "file ''$PWD/%p''\n" | sort) -c copy active_example.mp4')
!mv active_example.mp4 ../videos/
system(['rm ./*_seg*.mp4']);

%ffmpeg -framerate 30 -i 'explore_spk_lfp_%08d.tiff'  rest_jg05-20120312.mp4


% $$$ close(vidObj);
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'s'},lfs<4|(lfs>12&lfs<15),3)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'s'},lfs<4|(lfs>12&lfs<15),4)),2),...,...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'w'},lfs<4|(lfs>12&lfs<15),3)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'w'},lfs<4|(lfs>12&lfs<15),4)),2),...,...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ 
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'s'},lfs<4,2)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'s'},lfs<4,4)),2),...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,3)),2)./mean(log10(lys(stc{'w'},lfs<4,2)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2)./mean(log10(lys(stc{'w'},lfs<4,4)),2),...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ 
% $$$ 
% $$$ figure()
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'s'},lfs>5&lfs<10,4)),2),...
% $$$      mean(log10(lys(stc{'s'},lfs>5&lfs<10,5)),2),...
% $$$      '.');
% $$$ hold('on');
% $$$ plot(mean(log10(lys(stc{'w'},lfs>5&lfs<10,4)),2),...
% $$$      mean(log10(lys(stc{'w'},lfs>5&lfs<10,5)),2),...
% $$$      '.g');
% $$$ caxisLims = cat(1,xlim(),ylim());
% $$$ line([1,4],[1,4],'Color','k');
% $$$ xlim(caxisLims(1,:));
% $$$ ylim(caxisLims(2,:));
% $$$ 
% $$$ figure()
% $$$ hist2([mean(log10(lys(stc{'w+s'},lfs>5&lfs<10,3)),2),...
% $$$      mean(log10(lys(stc{'w+s'},lfs>5&lfs<10,4)),2)],...
% $$$      linspace(caxisLims(1,1),caxisLims(1,2),100),...
% $$$      linspace(caxisLims(2,1),caxisLims(2,2),100));
% $$$ 
% $$$ 
