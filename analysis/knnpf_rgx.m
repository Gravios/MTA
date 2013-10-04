

%% Primary Goal - knn estimation of unit firing rate
s = MTASession('jg05-20120317');

Trial = MTATrial(s,{{'CluRes',s.xyzSampleRate}},'all');

Trial = Trial.filter();

klen = 32;
kern = ones(klen,1);
overlap = 8;
nnn = 81;
nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';


myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);

newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));

myrx = SelectPeriods(myxyz,round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
mywx = SelectPeriods(myxyz,round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);



%unit = 1;
ufrs = {};
%%iu = [12,14,24];
%while unit ~=-1
for unit = 1:76,

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

myru = SelectPeriods(myufr,round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
mywu = SelectPeriods(myufr,round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);

pfknnmrw{unit} = nan(nxbins,nybins,length(ufrs{unit}));
pfknnmdw{unit} = nan(nxbins,nybins);
pfknnmrr{unit} = nan(nxbins,nybins,length(ufrs{unit}));
pfknnmdr{unit} = nan(nxbins,nybins);

xyoccw = xyocc(Trial,mywx,50)';
xyoccwt = sum(xyoccw(:));
xyoccr = xyocc(Trial,myrx,50)';
xyoccrt = sum(xyoccr(:));

tic
for u = 1:length(ufrs{unit})
puw = sum(mywu==ufrs{unit}(u))/length(mywu);
pur = sum(myru==ufrs{unit}(u))/length(myru);
for y = 1:length(ybins),
    for x = 1:length(xbins),
        distw = sqrt(sum((mywx-repmat([xbins(x),ybins(y)],size(mywx,1),1)).^2,2));
        distr = sqrt(sum((myrx-repmat([xbins(x),ybins(y)],size(myrx,1),1)).^2,2));
        [~,distIndw ] = sort(distw);        
        [~,distIndr] = sort(distr);        
        pfknnmdw{unit}(y,x) = median(distw(distIndw(1:nnn)));
        pfknnmdr{unit}(y,x) = median(distr(distIndr(1:nnn)));
        if pfknnmdw{unit}(y,x)<50,
            %pfknnmrw(y,x) = min(mywu(distIndw(1:nnn)));
            %pfknnmrw(y,x) = max(mywu(distIndw(1:nnn)));
            %pfknnmrw(y,x) = mean(mywu(distIndw(1:nnn)));
            tu = mywu(distIndw(1:nnn));
            %pfknnmrw(y,x) = mean(tu(tu<40&tu>1));            
            pfknnmrw{unit}(y,x,u) = (sum(tu==ufrs{unit}(u))/nnn)*((xyoccw(y,x)/xyoccwt)/puw);
        end
        if pfknnmdr{unit}(y,x)<50,
            %pfknnmrr(y,x) = min(myru(distIndr(1:nnn)));
            %pfknnmrr(y,x) = max(myru(distIndr(1:nnn)));
            %pfknnmrr(y,x) = mean(myru(distIndr(1:nnn)));
            tu = myru(distIndr(1:nnn));
            %pfknnmrr(y,x) = mean(tu(tu<40&tu>1));
            pfknnmrr{unit}(y,x,u) = (sum(tu==ufrs{unit}(u))/nnn)*((xyoccw(y,x)/xyoccwt)/puw);
        end
    end
end
end
toc

end
save([Trial.spath.analysis Trial.filebase '.bayes2.mat'],'pfknnmrw','pfknnmdw','pfknnmrr','pfknnmdr','-v7.3');


%pfknnmrw = pfknnmrw./length(myru);
%pfknnmrr = pfknnmrr./length(myru);

% $$$ 
% $$$ %% Permutation Analysis
% $$$ niter = 1000;
% $$$ %% initialize rate maps
% $$$ pfknnmrs = nan(nxbins,nybins,niter);
% $$$ pfknnmds = nan(nxbins,nybins,niter);
% $$$ 
% $$$ wSampSize = length(mywx);
% $$$ rSampSize = length(myrx);
% $$$ wSubSampSize = round(wSampSize/(wSampSize+rSampSize)*rSampSize);
% $$$ rSubSampSize = round(rSampSize/(wSampSize+rSampSize)*rSampSize);
% $$$ 
% $$$ for i = 1:niter,
% $$$ fprintf('%i\n',i)
% $$$ tic
% $$$     %% generate permutation
% $$$     wind = randi(wSampSize,[wSubSampSize,1]);
% $$$     rind = randi(rSampSize,[rSubSampSize,1]);
% $$$     mysu = cat(1,mywu(wind),myru(rind));
% $$$     mysx = cat(1,mywx(wind,:),myrx(rind,:));
% $$$     for y = 1:length(ybins),
% $$$         for x = 1:length(xbins),
% $$$             dists = sqrt(sum((mysx-repmat([xbins(x),ybins(y)],size(mysx,1),1)).^2,2));
% $$$             [~,distInds ] = sort(dists);        
% $$$             pfknnmds(y,x,i) = median(dists(distInds(1:nnn)));
% $$$             if pfknnmds(y,x,i)<50,
% $$$                 pfknnmrs(y,x,i) = mean(mysu(distInds(1:nnn)));
% $$$             end
% $$$         end
% $$$     end
% $$$ toc
% $$$ end

b=6;
subplot(231);
%imagescnan({xbins,ybins,pfknnmrs},[],[],1,[0,0,0]);,axis xy
hist(mywu,15);
xl = xlim;
xlim([-0.05,xl(2)]);
subplot(232);
imagescnan({xbins,ybins,pfknnmrw(:,:,b)},[],[],1,[0,0,0]);,axis xy
title(['walk unit ' num2str(Trial.map(unit,:))]);
subplot(233);
imagescnan({xbins,ybins,pfknnmrr(:,:,b)},[],[],1,[0,0,0]);,axis xy
title(['rear unit ' num2str(Trial.map(unit,:))]);
subplot(234);
bar(tbin,accg(:,unit)),axis tight,
subplot(235);
pfw.plot(unit,1);,
subplot(236);
pfr.plot(unit,1);,

Trial.printFig(gcf,'png',unit,'knnpf_100ms_50msol_1nn_walk_rear');
unit = figure_controls(gcf,unit)
%unit = unit+1;
end



xyoccw = xyocc(Trial,mywx,50);
xyoccr = xyocc(Trial,myrx,50);


[accg,tbin] = autoccg(Trial);
