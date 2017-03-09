function [uccg,tbins] = unit_ccg(Trial,units,states,spkmode)
%[uccg,tbins] = unit_ccg(Trial,units,states,spkmode)
Trial = 'jg05-20120310.cof.all';
Trial = MTATrial.validate(Trial);
Trial.load('stc','nn0317');


Trial.load('nq');
[units,states] = DefaultArgs(varargin,{find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19),'theta'});
nu = numel(units);



units = Trial.spk.map(:,1);
states = 'spw';

Spk = Trial.spk.copy;
Spk.create(Trial,Trial.sampleRate,states,units,'deburst');

spwPer = [Trial.stc{'spw'}];
spwUnitRate = [];
for u = units'
    spwUnitRate(u) = numel(Spk(u))/(sum(diff(spwPer.data,1,2)+1)/spwPer.sampleRate);
end
spwUnitRate(spwUnitRate==0)=nan;
figure,hist(log10(spwUnitRate),50)


tSpk = Trial.spk.copy;
tSpk.create(Trial,Trial.sampleRate,'theta',units,'deburst');
thetaPer = [Trial.stc{'theta'}];
thetaUnitRate = [];
for u = units'
    thetaUnitRate(u) = numel(tSpk(u))/(sum(diff(thetaPer.data,1,2)+1)/thetaPer.sampleRate);
end
thetaUnitRate(spwUnitRate==0)=nan;

figure,hist(log10(thetaUnitRate),50)

figure,plot(log10(thetaUnitRate),log10(spwUnitRate),'.');




halfbins = 25;
binsize = 160;
fwin = 11;

units = Spk.map(:,1);


[fccg,tbins] = CCG(Spk.res,Spk.clu,binsize,halfbins,Spk.sampleRate,units,'hz');

uccg.ccg = fccg;
uccg.fccg = reshape(Filter0(gausswin(fwin)./sum(gausswin(fwin)),fccg),size(fccg));
[uccg.cpow,uccg.tshf] = max(uccg.ccg);
uccg.cpow = sq(uccg.cpow);
uccg.tshf = tbins(sq(uccg.tshf));



%% get placefields
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';    
states = Trial.stc.list_state_attrib;
states(~cellfun(@isempty,regexp(states,'^spw$')))=[]; % drop the spw
states(~cellfun(@isempty,regexp(states,'^sit$')))=[]; % drop the spw
states(~cellfun(@isempty,regexp(states,'^groom$')))=[]; % drop the spw
nsts = numel(states);

binDims = [20,20];
overwrite = false;
numIter = 1;

%MTAApfs args
smoothingWeights = [2.2,2.2];
%MTAAknnpfs_bs args
nNearestNeighbors = 150;
distThreshold = 125;
ufrShufBlockSize = 1;
sampleRate = 15;

pfk = {};
pfs = {};    
for s = 1:nsts
    pfs{s} = MTAAknnpfs_bs(Trial,units,[states{s},'&theta'],overwrite, ...
                           'binDims',binDims,...
                           'nNearestNeighbors',nNearestNeighbors,...
                           'ufrShufBlockSize',ufrShufBlockSize,...
                           'distThreshold',distThreshold,...
                           'numIter',numIter);
    
% $$$     pfs{s} = MTAApfs(Trial,units,...
% $$$                      states{s},...
% $$$                      overwrite, ...
% $$$                      'binDims',binDims,...
% $$$                      'SmoothingWeights',smoothingWeights,...
% $$$                      'numIter',numIter);

    fprintf('pfk %s: complete\n',states{s});
end

mRate = pfs{7}.maxRate(units,true);

[accg,tbin] = autoccg(Trial,units,'theta');
[sccg,tbin] = autoccg(Trial,[],'spw');

mxr = [];
mxp = [];
for s = 1:nsts
    [mxr(s,:),mxp(s,:,:)] = pfs{s}.maxRate(units,true);
end

width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;

mxd = sqrt(sum(bsxfun(@minus,...
                      permute(mxp,[1,2,4,3]),...
                      permute(mxp,[1,4,2,3])).^2,4));


u = 45;
uind = mxd(7,:,u)<100&mxr(7,:)>2;

round(sum(uccg.ccg(:,uind,u)))
round(uccg.cpow(uind,uind))


unit = u;
hfig = figure(2039320);clf
hfig.Units = 'centimeters';
hfig.Position = [0,5,55,10.5];
for s = 1:nsts
    subplot(3,nsts,s); hold('on');
    pf = pfs{s};
    % Correct color of nans and plot place field
    ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,1),fliplr(pf.adata.binSizes')).*mask;
    ratemap(isnan(ratemap)) = -1;
    imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
    text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
         sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
    %colormap([0,0,0;parula]);
    caxis([-1,sq(mRate(unit==units)).*1.5]);    
    %caxis([-1,8]);        
    title(pf.parameters.states);    
end


cunits = find(uind);

for t = 1:numel(cunits),
    subplot(3,nsts,nsts+t)
    bar(tbins,fccg(:,cunits(t),cunits(t)))
    xlim([tbin([1,end])])    
    title(num2str(cunits(t)));
end

for t = 1:numel(cunits),
    subplot(3,nsts,nsts.*2+t)
    bar(tbins,fccg(:,u,cunits(t)))
    xlim([tbin([1,end])])    
    title(num2str(cunits(t)));
end


sunits = cunits;


hfig = figure(2039321);clf
hfig.Units = 'centimeters';
hfig.Position = [0,0,55,3.5*numel(sunits)];
for t = 1:numel(sunits)
    for s = 1:nsts
        unit = sunits(t);
        subplot(numel(sunits),nsts,s+nsts*(t-1)); hold('on');
        pf = pfs{s};
        % Correct color of nans and plot place field
        ratemap = reshape(pf.data.rateMap(:,unit==pf.data.clu,1),fliplr(pf.adata.binSizes')).*mask;
        ratemap(isnan(ratemap)) = -1;
        imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap);    
        text(pf.adata.bins{1}(end)-350,pf.adata.bins{2}(end)-50,...
             sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
        %colormap([0,0,0;parula]);
        caxis([-1,sq(mRate(unit==units)).*1.5]);    
        %caxis([-1,8]);        
        title(pf.parameters.states);    
        if s==1,ylabel(num2str(unit));end
    end
end




% $$$ figure
% $$$ subplot(121),imagesc(log10(uccg.cpow)')
% $$$ subplot(122),imagesc(mut')
% $$$ 
% $$$ states = {'theta','rear&theta','walk&theta'};
% $$$ numsts = numel(states);
% $$$ pfs={};
% $$$ for i = 1:numsts,
% $$$     pfs{i}  =  MTAAknnpfs(Trial,[],states{i},0,'numIter',1, ...
% $$$                           'ufrShufBlockSize',0,'binDims',[20, ...
% $$$                         20],'distThreshold',70); ...
% $$$                                                                                                    
% $$$ end




% $$$ unt = [18,91];
% $$$ figure,
% $$$ for s = 1:3,
% $$$ subplot2(2,3,1,s),pfs{s}.plot(pfs{s}.data.clu(unit(1)));
% $$$ subplot2(2,3,2,s),pfs{s}.plot(pfs{s}.data.clu(unit(2)));
% $$$ end
% $$$ uic = mat2cell(unit,[1],[1,1]);
% $$$ muc(uic{:})
% $$$ mut(uic{:})
% $$$ 
% $$$ xyc = round(get(gca,'CurrentPoint'));
% $$$ unit = [18,91];
% $$$ figure,
% $$$ for s = 1:3,
% $$$ subplot2(2,3,1,s),pfs{s}.plot(88);
% $$$ subplot2(2,3,2,s),pfs{s}.plot(35);
% $$$ end
% $$$ uic = mat2cell(unit,[1],[1,1]);
% $$$ muc(uic{:})
% $$$ mut(uic{:})
% $$$ 
% $$$ 
% $$$ figure
% $$$ nu = 10;
% $$$ unit = 22:31;
% $$$ for u = 1:nu,
% $$$ for o = 1:nu,
% $$$ subplot2(nu,nu,u,o);
% $$$ bar(t,uccg(:,unit(u),unit(o)));
% $$$ end
% $$$ end




% $$$ %% Diagnostics checking the status of xyz head position X pitch relationship
% $$$ xyz = Trial.load('xyz','seh');
% $$$ xyz = Trial.load('xyz');
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ 
% $$$ 
% $$$ edx = linspace(1.5,2.5,100);
% $$$ edy = linspace(-pi/2,pi/2,100);
% $$$ 
% $$$ figure
% $$$ ind = [Trial.stc{'hwalk'}];
% $$$ hist2([log10(xyz(ind,'head_back',3)),ang(ind,'head_back','head_front',2)],edx,edy);
% $$$ figure
% $$$ ind = [Trial.stc{'lwalk'}];
% $$$ hist2([log10(xyz(ind,'head_back',3)),ang(ind,'head_back','head_front',2)],edx,edy);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,hold on,
% $$$ eds = linspace(1.5,2.5,100);
% $$$ ind = [Trial.stc{'hwalk'}];
% $$$ hs = bar(eds,histc(log10(xyz(ind,'head_back',3)),eds),'histc');
% $$$ hs.FaceAlpha =0.5;
% $$$ hs.EdgeAlpha =0.5;
% $$$ hs.FaceColor ='c';
% $$$ hs.EdgeColor ='c';
% $$$ ind = [Trial.stc{'lwalk'}];
% $$$ hs = bar(eds,histc(log10(xyz(ind,'head_back',3)),eds),'histc');
% $$$ hs.FaceAlpha =0.5;
% $$$ hs.EdgeAlpha =0.5;
% $$$ hs.FaceColor ='r';
% $$$ hs.EdgeColor ='r';
