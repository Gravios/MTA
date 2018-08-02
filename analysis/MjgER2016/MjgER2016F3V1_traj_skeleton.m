

% LOAD Session data --------------------------------------------------------------------------------
MjgER2016_load_data();

% LOAD behavioral data
[pfd ,tags ,eigVec, eigVar, eigScore, validDims,...
 unitSubsets, unitIntersection, zrmMean, zrmStd] = req20180123_ver5(Trials);
numComp = size(eigVec{1},2);
pfindex = 1;

% LOAD behavioral scores
MjgER2016_load_bhv_erpPCA_scores();


% RASTER STUFF
Trial = Trials{20};    % jg05-20120312.cof.all
unitsPyr  = units{20};     % good units
unitsInt = select_units(Trial,'int');
xyz = preproc_xyz(Trial,'trb');


vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();

pft = pfs_2d_theta(Trial);

pyr = Trial.load('spk',Trial.lfp.sampleRate,'',unitsPyr);
int = Trial.load('spk',Trial.lfp.sampleRate,'',unitsInt);

lfp = Trial.load('lfp',[68,72,76,82]);

ufrInt = Trial.load('ufr',[],[],unitsInt,0.5);
ufrPyr = Trial.load('ufr',xyz,[],unitsPyr,0.75);

drz = xyz.copy();
drz.data = compute_drz(Trial,unitsPyr,pft);

xyl = xyz.copy();
xyl.data = sq(xyz(:,'nose',[1,2]));
xyl.resample(Trial.lfp.sampleRate);


bhvs = {'rear','hloc+hpause','lloc+lpause'};
ufrMeanBhvRate = zeros([numel(unitSet),numel(bhvs)]);
for s = 1:numel(bhvs),
    ufrMeanBhvRate(s,:) = mean(ufrPyr(stc(bhvs{s}),:));
end
[pyrSortedBhv,indPyrSortedBhv] = sortrows(bsxfun(@rdivide,ufrMeanBhvRate(:,[1,3]),max(ufrMeanBhvRate(:,[1,3]),[],2)));


cluSessionMapSubset = cluSessionMap(unitSubsets{1},:);

unitSubset = cluSessionMapSubset(cluSessionMapSubset(:,1)==20,2);
unitBhvScore = FSrC(cluSessionMapSubset(:,1)==20,1:3);
unitBhvSig  = fsrcz(cluSessionMapSubset(:,1)==20,1:3);

% END LOAD DATA ------------------------------------------------------------------------------------




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
cmapLims = [1,2.5;1,2.5;1,3;1,3.5];
for s = 1:4,
subplot(7,1,s);
imagesc(ts,fs,log10(ys(:,:,s))');
axis('xy');
caxis(cmapLims(s,:));
colormap(gca,'jet');
end


subplot(7,1,5);cla
hold('on');
plot(ts,mean(log10(ys(:,[1:5],1)),2));
plot(ts,mean(log10(ys(:,[1:30],1)),2));
%plot(ts,mean(log10(ys(:,[1:5],2)),2));
%plot(ts,mean(log10(ys(:,[1:5],3)),2));
%plot(ts,mean(log10(ys(:,[10:20],4)),2));
%plot(ts,mean(log10(ys(:,[1:5],1)),2)./mean(log10(ys(:,[1:5],4)),2));
plot(ts,mean(log10(ys(:,[10:20],3)),2)./mean(log10(ys(:,[10:20],4)),2),'m');
axis('xy');
Lines([],1.75,'k');
Lines([],0.75,'g');

subplot(716);cla();
v = vxy.data;
plot([1:size(vxy,1)]./xyz.sampleRate,v)
% $$$ v(v<=1e-3)=1e-3;
% $$$ plot([1:size(vxy,1)]./xyz.sampleRate,log10(v))

subplot(7,1,7);cla();
haxSTS = plotSTC(Trial.stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause','groom','sit'}),fliplr('rbcgkmy'));
xlim(indLFP([1,end])./lfp.sampleRate)
haxSTS.YTickLabel = fliplr({'rear','hloc','hpause','lloc','lpause','groom','sit'});
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');

figure();
subplot(8,1,[1:3]);
hold('on');
unitSet = unitsPyr(indPyrSortedBhv');
unitClr = repmat('b',[1,numel(unts)]);
for u = unitSet
    uind = find(u==unitSet);
    res = pyr(u);
    cts = 1:max(res);
    sWidth = 0.4./lfp.sampleRate;
    sti = ismember(cts,res);
    patch(reshape(repmat(reshape(bsxfun(@plus,cts(sti)/lfp.sampleRate,[-sWidth/2;sWidth/2]),[],1)',2,1),[],1),...
          repmat([0,1,1,0],[1,numel(res)])'+uind,...
          unitClr(uind),'EdgeColor',unitClr(uind));    
end
Lines([],1:numel(unitSet),'w');

subplot(8,1,4);
plot([[1:size(vxy,1)]'./vxy.sampleRate,[1:size(vxy,1)]'./vxy.sampleRate],vxy.data)
ylim([0,80]);

sp = subplot(8,1,5);
plotSTC(Trial.stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause','groom','sit'}),fliplr('rbcgkmy'));
ylim([0,8])


cmapLims = [1,3;1,3.5;1,4];
for s = 1:3,
subplot(8,1,s+5);
imagesc(ts,fs,log10(ys(:,:,s))');
axis('xy');
caxis(cmapLims(s,:));
colormap(gca,'jet');
end


linkaxes(get(gcf,'Children'),'x');

axBg = axes('Position',[0 0 1 1],'Visible','off');    
text(0.1,0.95,Trial.filebase);
Lines(0.5,[],'k');



% $$$ xlim(sp,[0,40]);
% $$$ for x = 0:0.25:round(ts(end))
% $$$     xlim(sp,[x,x+20]);
% $$$     drawnow();
% $$$ end
% $$$ 
% $$$ axes(sp);


vidObj = VideoWriter('/storage/share/Projects/BehaviorPlaceCode/decode/explore_spk_lfp_state.avi','Archival');
open(vidObj);
for x = 0:0.1:round(ts(end))
    xlim(sp,[x,x+40]);
    drawnow();
    writeVideo(vidObj,getframe());
end
close(vidObj);







% FIGURE Skeleton example -------------------------------------------------------------------------
%
% 
%
%

ang = create(MTADang,Trial,xyz);

% vars

indLFP = [935,956]*lfp.sampleRate

%indLFP = [1390,1425]*lfp.sampleRate ------
examplePoints = round([indLFP(6251),...
                    indLFP(round(numel(indLFP)./3)),...
                    indLFP(end-12000)]./lfp.sampleRate.*xyz.sampleRate);
% NOPES 
indLFP = [5675,5715]*lfp.sampleRate

% MAYBES 
indLFP = [2570,2600]*lfp.sampleRate
indLFP = [500,575]*lfp.sampleRate;
indLFP = [264,335]*lfp.sampleRate;

% YEPS 
indLFP = [6176,6237]*lfp.sampleRate;

% NOT YEPS
%indLFP = [1422,1465]*lfp.sampleRate
indLFP = [320,355]*lfp.sampleRate;
indLFP = [350,380]*lfp.sampleRate;


indLFP = [6882,6925]*lfp.sampleRate;

indLFP = [500,575]*lfp.sampleRate;

indLFP = [550,640]*lfp.sampleRate;
indLFP = [768,825]*lfp.sampleRate;
indLFP = [1010,1060]*lfp.sampleRate;
indLFP = [1350,1432]*lfp.sampleRate

indLFP = [1336,1390]*lfp.sampleRate
indLFP = [2575,2636]*lfp.sampleRate
indLFP = [3860,3935]*lfp.sampleRate
indLFP = [6806,6880]*lfp.sampleRate;

indLFP = indLFP(1):indLFP(2);
examplePoints = round([indLFP(6251),...
                    indLFP(round(numel(indLFP)./3)),...
                    indLFP(end-12000)]./lfp.sampleRate.*xyz.sampleRate);

indXYZ = round(indLFP(1)./lfp.sampleRate.*xyz.sampleRate:4:indLFP(end)./lfp.sampleRate.*xyz.sampleRate);
indXYZSub = indXYZ;
gridPft = cell(size(interpParPfs.bins));
[gridPft{:}] = ndgrid(interpParPfs.bins{:});
colorLimits = [min(xyl(indLFP,1)),max(xyl(indLFP,1));min(xyl(indLFP,2)),max(xyl(indLFP,2));0,1];



hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [2.5, 2.5, 24, 28];
% PLOT Trajectory in 3D space
subplot(8,1,[1,2,3]);
hold('on');
% $$$ scatter3(xyz(indXYZSub,'nose',1),...
% $$$          xyz(indXYZSub,'nose',2),...
% $$$          xyz(indXYZSub,'nose',3),10,[0.8,0.8,0.8],'filled');
scatter(xyz(indXYZSub,'nose',1),xyz(indXYZSub,'nose',2),10,...
        [(xyz(indXYZSub,'nose',1)-colorLimits(1,1))./(abs(diff(colorLimits(1,:)))),...
         (xyz(indXYZSub,'nose',2)-colorLimits(2,1))./(abs(diff(colorLimits(2,:)))),...
         zeros([numel(indXYZSub),1])],'filled');


%examplePoints = [indXYZ(220),indXYZ(990),indXYZ(600)];
for e = examplePoints, plotSkeleton(Trial,xyz,e,'surface',ang); end
camlight('left');
daspect([1,1,1]);



% SELECT units within example timeseries
unitPyrSpkCnt = [];
for unit = unitsPyr,
    unitPyrSpkCnt(end+1) = numel(SelectPeriods(pyr(unit),indLFP([1,end])));
end
% REORDER units from mean bhv ufr (rear to lloc)
unitTSub = unitsPyr(find(unitPyrSpkCnt>30));
unitTSub = unitTSub(ismember(unitTSub,unitSubset'));

% SORTROWS by behavior scores
tubs = unitBhvScore(:,[2,3,1]);
%tubs(tubs<-0.25) = -0.25;
%[bhvScr,bhvOrd] = sortrows(tubs,[1,2,-3]);
[~,stsOrd] = max(tubs,[],2);
[bhvMax,bhvOrd] = sort(stsOrd,'descend');
stubs = tubs(bhvOrd,:);
for s = 1:3, 
    [~,sortedMaxSubsetInd] = sort(stubs(bhvMax==s,s),'ascend');
    bhvOrdMaxSubset = bhvOrd(bhvMax==s);
    bhvOrd(bhvMax==s) = bhvOrdMaxSubset(sortedMaxSubsetInd);
end

utsordInd = ismember(unitSubset(bhvOrd)',unitTSub);
utsorded = unitSubset(bhvOrd);
unitTSub = utsorded(utsordInd)';

% MANUAL arrangement
% $$$ unitTSub = [unitTSub(1:3),unitTSub(18),unitTSub(4:17),unitTSub(19:end)];
% $$$ unitTSub = [unitTSub([16,1,3,2,4:6,13,11,9,10,12,14,15]),unitTSub(19:end),unitTSub(17:18)];
% $$$ % removed unit 8&7
% $$$ %extraUnits = unitSubset(~ismember(unitSubset,unitTSub));
% $$$ unitTSub = [102,118,136,73,unitTSub,145];



% GET DRZ subset
[~,sind] = sort(unitTSub);
[~,sind] = sort(sind);
upind = find(ismember(unitsPyr,unitTSub));
mdrz = drz.copy();
mdrz.data = drz(:,upind(sind));


% $$$ figure();
% $$$ for u = 1:16,
% $$$     hold('on');
% $$$     plot(pft,unitTSub(u));
% $$$     plot(xyz(indXYZ,'nose',1),xyz(indXYZ,'nose',2),'.w');
% $$$     plot(mrp(u,1),mrp(u,2),'*m');
% $$$     title(num2str(u));
% $$$     waitforbuttonpress();
% $$$ end
% $$$ close(gcf);

%listOfUnitsReqAltPos = [1,3,4,6,9,15];

% PLOT unit raster 
hax = subplot(8,1,[4:7]);
hold('on');
hax.Units = 'centimeters';
hax.Position(3:4) = [15, 10];
ylabel('Units');

rightTickLabel = [];
for unit = unitTSub,
    uind = find(unit==unitTSub);
    res = SelectPeriods(pyr(unit),indLFP([1,end]),'d',1,0);
    sWidth = 0.4;
    %sti = ismember(indLFP,res);
    patch([indXYZ(1),indXYZ,indXYZ(end)]./xyz.sampleRate,[0,1-abs(mdrz(indXYZ,uind))',0]+uind-1,...
          [0.75,0.75,0.75],'EdgeColor','none');
    for r = res',
        patch(reshape(repmat(reshape(bsxfun(@plus,r,[-sWidth/2;sWidth/2]),[],1)',2,1),[],1)./lfp.sampleRate,...
              [0.1,0.9,0.9,0.1]'+uind-1,[0,0,0],...
              'EdgeColor',([xyl(r,:),0]-colorLimits(:,1)')./abs(diff(colorLimits,1,2)'), ...
              'FaceColor','none');
    end    

    rightTickLabel(end+1,:) = round(unitBhvScore(unitSubset==unit,[2,3,1]),2);
end

hax.YTick = 0.5:[numel(unitTSub)+0.5];
hax.YTickLabel = 1:numel(unitTSub);


Lines([],1:numel(unitTSub),'k');
Lines(examplePoints./xyz.sampleRate,[],'k');
xlim(indLFP([1,end])./lfp.sampleRate)
grid('on');
sax = axes('Units','centimeters','Position',[hax.Position(1,2)+[12.5,0],4,hax.Position(4)]);
sax.Visible = 'off';

% SET grey scale limits for behavior f-score
colorLimitsBHV = [min(rightTickLabel);max(rightTickLabel)]

for unit = 1:size(rightTickLabel,1)
    for s = 1:size(rightTickLabel,2)
        bhvsColor = repmat((rightTickLabel(unit,s)-colorLimitsBHV(1,s))./(diff(colorLimitsBHV(:,s))*1.4),[1,3]);
        patch([s-1,s-1,s,s],...
              [unit-1,unit,unit,unit-1],...
              bhvsColor);
        if  bhvsColor(1)<0.6
            bhvColor = [0.1,0.1,0.1];
        elseif bhvsColor(1)>=0.6
            bhvColor = [0.9,0.9,0.9];
        end
        
        text(s-0.7,unit-0.5,sprintf('% 1.2f',rightTickLabel(unit,s)),'Color',1-bhvColor);
    end
end
xlim([0,3]);    

linkaxes([hax,sax],'y');
ylim([0,numel(unitTSub)]);



% PLOT BHV
subplot(8,1,[8]);
haxSTS = plotSTC(Trial.stc,1,[],fliplr({'rear','hloc','hpause','lloc','lpause'}),fliplr('rbcgk'));
xlim(indLFP([1,end])./lfp.sampleRate)
haxSTS.YTickLabel = fliplr({'rear','hloc','hpause','lloc','lpause'});
xlabel('Time (s)');
haxSTS.Units = 'centimeters';
haxSTS.Position(3) = [15];
grid('on');
Lines(examplePoints./xyz.sampleRate,[],'k');



% END FIGURE ---------------------------------------------------------------------------------------





%ufrIntSubset = ufrInt(stc{'x+r'},:);



[LU, LR, FSr, VT] = erpPCA(unity(ufrIntSubset(1:10:end,:)),15);
[U,S,V] = svd(unity(ufrIntSubset(1:10:end,:)),0);

binDims = [25,25];
bins = af(@(b) linspace(-500,500,b),  binDims);


binInds = cf(@(b,i)  discretize(xyl(:,i),b),  bins,mat2cell(1:size(xyl,2),1,ones([1,size(xyl,2)])));

stcm = sum(stc2mat(stc,xyl,{'walk','turn','pause','rear','groom','sit'}),2);

G = unique(int.clu);
pccg = zeros([binDims,numel(G),numel(G)]);
sccg = zeros([binDims,numel(G),numel(G),6,2]);
for indX = 1:binDims(1),
    for indY = 1:binDims(2),
        for sts = 1:6,
            pind = binInds{1}==indX&binInds{2}==indY&stcm==sts;
            
            myres = int.res(pind(int.res)==1);
            myclu = int.clu(pind(int.res)==1);

            [mccg,tbin] = CCG(myres,myclu,20,20,Trial.lfp.sampleRate,G,'hz');
            for j = 1:numel(G),
                for k = 1:numel(G),
                    try
                        [sccg(indX,indY,j,k,sts,1),sccg(indX,indY,j,k,sts,2)] = ...
                            max(spline(tbin,RectFilter(mccg(:,j,k),5,5),tbin));
                        sccg(indX,indY,j,k,sts,2) = tbin(sccg(indX,indY,j,k,sts,2));
                    end
                end
            end
            
            %pccg(indX,indY,:,:) = sq(max(mccg));
        end
    end
end

s = 4;
p = [5,15];
figure();
subplot(231);
imagesc(sccg(:,:,p(1),p(1),s,1)');  axis('xy');
subplot(232);
imagesc(sccg(:,:,p(1),p(2),s,1)');  axis('xy');
subplot(233);
imagesc(sccg(:,:,p(2),p(2),s,1)');  axis('xy');
sp = findobj(gcf,'Type','Axes');
af(@(a) caxis(a,[0,max(reshape(cell2mat(af(@(a) caxis(a), sp)),[],1))]), sp)
subplot(234);
imagesc(sccg(:,:,p(1),p(1),s,2))');  axis('xy');  caxis([-80,80]);
subplot(235);
imagesc(sccg(:,:,p(1),p(2),s,2)');  axis('xy');  caxis([-80,80]);
subplot(236);
imagesc(sccg(:,:,p(2),p(2),s,2)');  axis('xy');  caxis([-80,80]);


figure();imagesc(sccg(:,:,5,8,s,1)');  axis('xy');

