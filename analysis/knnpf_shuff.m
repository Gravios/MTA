

%% Primary Goal - knn estimation of unit firing rate
s = MTASession('jg05-20120317');

Trial = MTATrial(s,{{'CluRes',s.xyzSampleRate}},'all');

Trial = Trial.filter();

klen = 32;
kern = ones(klen,1);
overlap = 8;
nnn = 81;
dist_thresh = 75;
niter = 5000;
nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';


myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1,
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);

newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));

rper = round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-3;
wper = round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate)-3;

myrx = SelectPeriods(myxyz,rper,'c',1,1);
mywx = SelectPeriods(myxyz,wper,'c',1,1);



pfknnmrw = {};
pfknnmdw = {};
pfknnmrr = {};
pfknnmdr = {};

pfknnbsrw = {};
pfknnbsrr = {};
pfknnbsdr = {};
pfknnbsdw = {};

pfknnmrs = {};
pfknnmds = {};


% unit = 57;
ufrs = {};
% $$$ %%iu = [12,14,24];
% $$$ while unit ~=-1
for unit = gu(3:end),

myres = Trial.res(Trial.clu==unit);
myufr = zeros(size(Trial.xyz,1),1);
myufr(1:myres(end)) = accumarray(myres,ones(size(myres),1));
t =         permute(reshape(myufr(1:size(myufr,1)-mod(size(myufr,1),klen)),klen,[]),[3,1,2]);
for shift = 1:klen/overlap-1
t = cat(1,t,permute(reshape(circshift(myufr(1:size(myufr,1)-mod(size(myufr,1),klen)),-overlap*shift),klen,[]),[3,1,2]));
end
myufr = t; 
myufr = reshape(sq(sum(repmat(permute(repmat(kern,1,size(myufr,3)),[3,1,2]),size(myufr,1),1).*myufr,2))/(klen/Trial.xyzSampleRate),[],1);

try,ufrs{unit} = unique(myufr);end

myru = SelectPeriods(myufr,rper,'c',1,1);
mywu = SelectPeriods(myufr,wper,'c',1,1);

pfknnmrw{unit} = nan(nxbins,nybins,nnn);
pfknnmdw{unit} = nan(nxbins,nybins,nnn);
pfknnmrr{unit} = nan(nxbins,nybins,nnn);
pfknnmdr{unit} = nan(nxbins,nybins,nnn);

tic
for y = 1:length(ybins),
    for x = 1:length(xbins),
        distw = sqrt(sum((mywx-repmat([xbins(x),ybins(y)],size(mywx,1),1)).^2,2));
        distr = sqrt(sum((myrx-repmat([xbins(x),ybins(y)],size(myrx,1),1)).^2,2));
        [~,distIndw ] = sort(distw);        
        [~,distIndr] = sort(distr);        
        pfknnmdw{unit}(y,x) = median(distw(distIndw(1:nnn)));
        pfknnmdr{unit}(y,x) = median(distr(distIndr(1:nnn)));
        if pfknnmdw{unit}(y,x)<dist_thresh,
            %pfknnmrw{unit}(y,x) = min(mywu(distIndw(1:nnn)));
            %pfknnmrw{unit}(y,x) = max(mywu(distIndw(1:nnn)));
            pfknnmrw{unit}(y,x,:) = mywu(distIndw(1:nnn));            
        end
        if pfknnmdr{unit}(y,x)<dist_thresh,
            %pfknnmrr{unit}(y,x) = min(myru(distIndr(1:nnn)));
            %pfknnmrr{unit}(y,x) = max(myru(distIndr(1:nnn)));
            pfknnmrr{unit}(y,x,:) = myru(distIndr(1:nnn));
        end
    end
end
toc

%% Bootstrap Analysis


pfknnbsrw{unit} = nan(nxbins,nybins,niter);
pfknnbsrr{unit} = nan(nxbins,nybins,niter);
pfknnbsdr = nan(nxbins,nybins);
pfknnbsdw = nan(nxbins,nybins);


tic
for y = 1:length(ybins),
    for x = 1:length(xbins),
        distw = sqrt(sum((mywx-repmat([xbins(x),ybins(y)],size(mywx,1),1)).^2,2));
        distr = sqrt(sum((myrx-repmat([xbins(x),ybins(y)],size(myrx,1),1)).^2,2));
        [~,distIndw ] = sort(distw);        
        [~,distIndr] = sort(distr);        
        pfknnbsdw = median(distw(distIndw(1:nnn)));
        pfknnbsdr = median(distr(distIndr(1:nnn)));
        for i = 1:niter
            nnnind = randi(nnn,[nnn,1]);
            if pfknnbsdw<dist_thresh,
                %pfknnmrw{unit}(y,x) = min(mywu(distIndw(1:nnn)));
                %pfknnmrw{unit}(y,x) = max(mywu(distIndw(1:nnn)));
                pfknnbsrw{unit}(y,x,i) = mean(mywu(distIndw(nnnind)));            
            end
            if pfknnbsdr<dist_thresh,
                %pfknnmrr{unit}(y,x) = min(myru(distIndr(1:nnn)));
                %pfknnmrr{unit}(y,x) = max(myru(distIndr(1:nnn)));
                pfknnbsrr{unit}(y,x,:) = mean(myru(distIndr(nnnind)));
            end
        end
    end
end
toc


% $$$ %% Permutation Analysis
% $$$ 
% $$$ %% initialize rate maps
% $$$ pfknnmrs{unit} = nan(nxbins,nybins,niter);
% $$$ pfknnmds{unit} = nan(nxbins,nybins,niter);
% $$$ 
% $$$ wSampSize = length(mywx);
% $$$ rSampSize = length(myrx);
% $$$ wSubSampSize = round(wSampSize/(wSampSize+rSampSize)*rSampSize);
% $$$ rSubSampSize = round(rSampSize/(wSampSize+rSampSize)*rSampSize);
% $$$ 
% $$$ wind = randi(wSampSize,[wSubSampSize,niter]);
% $$$ rind = randi(rSampSize,[rSubSampSize,niter]);
% $$$ 
% $$$ tic
% $$$ for y = 1:length(ybins),
% $$$     for x = 1:length(xbins),
% $$$         for i = 1:niter,
% $$$             %% generate permutation
% $$$             mysu = cat(1,mywu(wind(:,i)),myru(rind(:,i)));
% $$$             mysx = cat(1,mywx(wind(:,i),:),myrx(rind(:,i),:));
% $$$ 
% $$$             dists = sqrt(sum((mysx-repmat([xbins(x),ybins(y)],size(mysx,1),1)).^2,2));
% $$$             [~,distInds ] = sort(dists);        
% $$$             pfknnmds{unit}(y,x,i) = median(dists(distInds(1:nnn)));
% $$$             if pfknnmds{unit}(y,x,i)<dist_thresh,
% $$$                 pfknnmrs{unit}(y,x,i) = mean(mysu(distInds(1:nnn)));
% $$$             end
% $$$         end
% $$$     end
% $$$ end
% $$$ toc

end

% $$$ subplot(231);
% $$$ %imagescnan({xbins,ybins,pfknnmrs},[],[],1,[0,0,0]);,axis xy
% $$$ hist(mywu,15);
% $$$ xl = xlim;
% $$$ xlim([-0.05,xl(2)]);
% $$$ subplot(232);
% $$$ imagescnan({xbins,ybins,pfknnmrw},[],[],1,[0,0,0]);,axis xy
% $$$ title(['walk unit ' num2str(Trial.map(unit,:))]);
% $$$ subplot(233);
% $$$ imagescnan({xbins,ybins,pfknnmrr},[],[],1,[0,0,0]);,axis xy
% $$$ title(['rear unit ' num2str(Trial.map(unit,:))]);
% $$$ subplot(234);
% $$$ bar(tbin,accg(:,unit)),axis tight,
% $$$ subplot(235);
% $$$ pfw.plot(unit,1);,
% $$$ subplot(236);
% $$$ pfr.plot(unit,1);,
% $$$ 
% $$$ Trial.printFig(gcf,'png',unit,'knnpf_100ms_50msol_1nn_walk_rear');
% $$$ unit = figure_controls(gcf,unit)
% $$$ %unit = unit+1;




% $$$ 
% $$$ xyoccw = xyocc(Trial,SelectPeriods(Trial.xyz(:,7,[1,2]),Trial.Bhv.getState('walk').state,'c',1,1),50);
% $$$ xyoccr = xyocc(Trial,SelectPeriods(Trial.xyz(:,7,[1,2]),Trial.Bhv.getState('rear').state,'c',1,1),50);
% $$$ 
% $$$ 
% $$$ [accg,tbin] = autoccg(Trial);
% $$$ 
% $$$ figure,imagescnan({xbins,ybins,pfknnmrs(:,:,1)},[],[],1,[0,0,0]);,axis xy
% $$$ 
% $$$ 
% $$$ nnn = 81
% $$$ pfknnmdmw = zeros(50,50);
% $$$ pfknnmrr = nan(nxbins,nybins,nnn);
% $$$ pfknnmrw = nan(nxbins,nybins,nnn);
% $$$ tic
% $$$ for y = 1:length(ybins),
% $$$     for x = 1:length(xbins),
% $$$         distw = sqrt(sum((mywx-repmat([xbins(x),ybins(y)],size(mywx,1),1)).^2,2));
% $$$         distr = sqrt(sum((myrx-repmat([xbins(x),ybins(y)],size(myrx,1),1)).^2,2));
% $$$         [~,distIndw ] = sort(distw);        
% $$$         [~,distIndr] = sort(distr);        
% $$$         pfknnmdw(y,x) = median(distw(distIndw(1:nnn)));
% $$$         pfknnmdr(y,x) = median(distr(distIndr(1:nnn)));
% $$$         pfknnmdmw(y,x) = max(distw(distIndw(1:nnn)));
% $$$         pfknnmdmr(y,x) = max(distr(distIndr(1:nnn)));
% $$$         if pfknnmdw(y,x)<150,
% $$$             %pfknnmrw(y,x) = min(mywu(distIndw(1:nnn)));
% $$$             %pfknnmrw(y,x) = max(mywu(distIndw(1:nnn)));
% $$$             pfknnmrw(y,x,:) = mywu(distIndw(1:nnn));
% $$$         end
% $$$         if pfknnmdr(y,x)<150,
% $$$             %pfknnmrr(y,x) = min(myru(distIndr(1:nnn)));
% $$$             %pfknnmrr(y,x) = max(myru(distIndr(1:nnn)));
% $$$             pfknnmrr(y,x,:) = myru(distIndr(1:nnn));
% $$$         end
% $$$     end
% $$$ end
% $$$ toc
% $$$ 
% $$$ 
% $$$ figure
% $$$ subplot(231);
% $$$ %imagescnan({xbins,ybins,pfknnmrs},[],[],1,[0,0,0]);,axis xy
% $$$ hist(mywu,15);
% $$$ xl = xlim;
% $$$ xlim([-0.05,xl(2)]);
% $$$ subplot(232);
% $$$ imagescnan({xbins,ybins,mean(pfknnmrw,3)},[],[],1,[0,0,0]);,axis xy
% $$$ title(['walk unit ' num2str(Trial.map(unit,:))]);
% $$$ subplot(233);
% $$$ imagescnan({xbins,ybins,mean(pfknnmrr,3)},[],[],1,[0,0,0]);,axis xy
% $$$ title(['rear unit ' num2str(Trial.map(unit,:))]);
% $$$ subplot(234);
% $$$ bar(tbin,accg(:,unit)),axis tight,
% $$$ subplot(235);
% $$$ pfw.plot(unit,1);,
% $$$ subplot(236);
% $$$ pfr.plot(unit,1);,
% $$$ 
% $$$ 
% $$$ 
% $$$ mp = mean(pfknnmrw,3),
% $$$ mp(isnan(mp))=0; 
% $$$ mind = LocalMinima2(-mp,-5,5);
% $$$ figure,hist(sq(pfknnmrw(mind(1,1),mind(1,2),:)),100)
% $$$ figure,hist(sq(pfknnmrw(mind(2,1),mind(2,2),:)),100)
% $$$ figure,hist(sq(pfknnmrw(mind(3,1),mind(3,2),:)),100)
% $$$ 
% $$$ 
% $$$ 
% $$$ mp = mean(pfknnmrr,3);
% $$$ mp(isnan(mp))=0; 
% $$$ mind = LocalMinima2(-mp,-1,5);
% $$$ figure,hist(sq(pfknnmrr(mind(1,1),mind(1,2),:)),100)
% $$$ figure,hist(sq(pfknnmrr(mind(2,1),mind(2,2),:)),100)
% $$$ figure,hist(sq(pfknnmrr(mind(3,1),mind(3,2),:)),100)
% $$$ 
% $$$ 
% $$$ 
% $$$ rearsegs = GetSegs(sq(Trial.xyz(:,7,:)),Trial.Bhv.getState('rear').state(:,1),80);
% $$$ srs = rearsegs-repmat(rearsegs(1,:,:),size(rearsegs,1),1);
% $$$ c = jet(size(rearsegs,2));
% $$$ 
% $$$ figure,
% $$$ hold on
% $$$ for i = 1:size(rearsegs,2),
% $$$ pa =plot3(srs(:,i,1),srs(:,i,2),srs(:,i,3),'color',c(i,:));
% $$$ end




% $$$ pfw = MTAPlaceField(Trial,[],'walk');
% $$$ pfr = MTAPlaceField(Trial,[],'rear');
% $$$ gu = [];
% $$$ for index =1:76,
% $$$ subplot(211),
% $$$ pfw.plot(index,1),
% $$$ subplot(212),
% $$$ pfr.plot(index,1),
% $$$ waitforbuttonpress
% $$$ whatkey = get(gcf,'CurrentCharacter');
% $$$ switch double(whatkey)
% $$$   case double('n')
% $$$     index = index+1;
% $$$   case double('m')
% $$$     gu(end+1) = index;
% $$$   case double('i')
% $$$     index = input('Enter index #: ');
% $$$   case double('p')
% $$$     index=index-1;
% $$$   case double('q')
% $$$     index = -1;
% $$$ end
% $$$ end



pval = 1./sum((permute(repmat(permute(mean(pfknnmrw{unit},3),circshift([1:ndims(pfknnmrw{unit})]',1)'),niter,1),circshift([1:ndims(pfknnmrw{unit})]',-1)')>pfknnbsrw{unit}),3);
pval(isinf(pval)) = nan;

figure
imagescnan({xbins,ybins,pval},[],[],1,[0,0,0]);
figure
imagescnan({xbins,ybins,mean(pfknnmrw{unit},3)},[],[],1,[0,0,0]);
