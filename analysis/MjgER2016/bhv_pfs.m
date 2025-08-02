OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';    
Trial= MTATrial('jg05-20120310');    
Trial.load('stc','nn0317');
states = Trial.stc.list_state_attrib;
states(~cellfun(@isempty,regexp(states,'^spw$')))=[]; % drop the spw    

binDims = [20,20];    
units = [];
overwrite = false;
numIter = 1001;

smoothingWeights = [2.2,2.2];

nNearestNeighbors = 300;
distThreshold = 125;
ufrShufBlockSize = 1;
sampleRate = 30;



pfk = {};
pfs = {};    
for s = 1:numel(states)
    pfk{s} = MTAAknnpfs_bs(Trial,units,states{s},overwrite, ...
                           'binDims',binDims,...
                           'nNearestNeighbors',nNearestNeighbors,...
                           'ufrShufBlockSize',ufrShufBlockSize,...
                           'distThreshold',distThreshold,...
                           'numIter',numIter);
% $$$         pfs{s} = MTAApfs(Trial,units,...
% $$$                          states{s},...
% $$$                          overwrite, ...
% $$$                          'binDims',binDims,...
% $$$                          'SmoothingWeights',smoothingWeights,...
% $$$                          'numIter',numIter);
    fprintf('pfk %s: complete\n',states{s});
end
units = pfk{1}.data.clu;    
mRate = pfk{7}.maxRate(units);

[accg,tbin] = autoccg(Trial,units,'theta');



width = pfk{1}.adata.binSizes(1);
height = pfk{1}.adata.binSizes(2);
radius = round(pfk{1}.adata.binSizes(1)/2)-find(pfk{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;


%FigDir = '/jg05-20120310_pfk_bhv';
FigDir = '/jg05-20120310_pfs_bhv';
FigPrefix = 'pfs_';
try,mkdir(fullfile(OwnDir,FigDir));end



autoincr = true;
%autoincr = false;    
hfig = figure(849274);clf
unit = units(1);
set(hfig,'units','centimeters')
set(hfig,'Position',[26,2,25,24])
set(hfig,'PaperPositionMode','auto');

while unit~=-1,
    for s = 1:numel(states)
        subplot(4,4,s)
        hold('on')            
        pf = pfk{s};
        %pf = pfs{s};

        % Correct color of nans and plot place field

        % PFS plot
        %ratemap = pf.plot(unit,'mazeMaskFlag',false);

        % PFK plot
        %ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,1),fliplr(pf.adata.binSizes'));
        ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,2:end),...
                          [fliplr(pf.adata.binSizes'),1000]);
        
        if size(ratemap,3)>1
            ratemap = mean(ratemap,3)'.*mask;
        else
            ratemap = ratemap.*mask;
        end

        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);

        text(pf.adata.bins{1}(end)-250,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        colormap([0,0,0;parula]);
        caxis([-1,mRate(unit==units)]);        

        title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);
    end   

    subplot(4,4,16)
    bar(tbin,accg(:,unit));
    xlim([min(tbin),max(tbin)]);
    title([' AutoCCG: Unit ',num2str(unit)]);

    print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigPrefix,num2str(unit),'.eps']));
    print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigPrefix,num2str(unit),'.png']));

    
    unit = figure_controls(hfig,unit,units,autoincr);
    
end



xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
eds = linspace([40,190,100]);

figure,hold on
ind = Trial.stc{'q'};    
hs = bar(eds,histc(xyz(ind,5,3),eds),'histc');
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
ind = Trial.stc{'j'};    
hs = bar(eds,histc(xyz(ind,5,3),eds),'histc');
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
hs.FaceColor = 'r';
hs.EdgeColor = 'r';

eds = linspace([-pi/2,pi/2,100]);    
figure,hold on
ind = Trial.stc{'q'};    
hs = bar(eds,histc(ang(ind,5,7,2),eds),'histc');
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
hs.FaceColor = 'c';
hs.EdgeColor = 'c';
ind = Trial.stc{'j'};    
hs = bar(eds,histc(ang(ind,5,7,2),eds),'histc');
hs.FaceAlpha = 0.5;
hs.EdgeAlpha = 0.5;
hs.FaceColor = 'r';
hs.EdgeColor = 'r';

unit = 34;
ratemap = reshape(pfk{1}.data.rateMap(:,unit==pfk{1}.data.clu,2:end),...
                  [fliplr(pfk{1}.adata.binSizes'),1000]);
oratemap = reshape(pfk{1}.data.rateMap(:,unit==pfk{1}.data.clu,1),...
                   [fliplr(pfk{1}.adata.binSizes')]);
figure,
subplot(131);
imagesc((mean(ratemap,3)'./std(ratemap,[],3)').*mask),axis xy,
subplot(132);
imagesc(mean(ratemap,3)'.*mask),axis xy,
subplot(133);
imagesc(std(ratemap,[],3)'.*mask),axis xy,

figure
subplot(131);
imagesc(oratemap'.*mask),axis xy
subplot(132);
imagesc(mean(ratemap,3)'.*mask),axis xy,
subplot(133);
imagesc((oratemap'-mean(ratemap,3)').*mask),axis xy,

