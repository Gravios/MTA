;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

% While the rat sits, a theta state exists where only CA1lm theta is present\


% load depth profile
Trial = MTATrial.validate('jg05-20120312.cof.all');
stc = Trial.load('stc','msnn_ppsvd_raux');
channels = [65:96];
%lfp = Trial.load('lfp',

% $$$ sampleRate = 250;

try, lfp = Trial.load('lfp',channels);
catch, lfp = Trial.load('lfp',channels);
end
% $$$ phz = lfp.phase([5,13]);    
% $$$ phz.data = unwrap(phz.data);
% $$$ phz.resample(xyz);    
% $$$ phz.data = mod(phz.data+pi,2*pi)-pi;
% $$$ lfp.resample(xyz);

flfp = lfp.copy();
flfp.filter('ButFilter',4,[6,15],'bandpass');
flfp.resample(ys);

specArgs = struct('nFFT',2^8,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^7,...
                  'nOverlap',2^7*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[100,300]);
   
[ys,fs,ts] = fet_spec(Trial,lfp,[],[],[],specArgs);




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
ufr = ufr.create(Trial,xyz,'',unitSubset,0.04,true,'gauss');

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
imagesc(ts,1:32,log10(sq(ys(:,21,:)))');
caxis([0,4])


% Identification of the pyramidal layer of the CA1-3
% Emperical Identification of probe center offset from pyramidal layer.
% Compute spectral power in the ripple band
% Smooth along time and space;
% Angle of incidence relative to the surface (plane) of the pyramidal layer.
% The top site of the probe ( assume the approach is from the cortex to the dorsal hippocampus )


probeDepth = 7
targetLayerChan = 7;
% spectral features
% non-whitened spectral power

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
    imagesc(ts,1:32,imgaussfilt(log10(sq(ys(:,20,:))),[10,1])');
    colormap('jet');
    caxis([0,5]);
figure,imagesc(ts,1:32,bsxfun(@rdivide,sq(ys(:,20,1:15)),sum(sq(ys(:,20,1:15)),2))');



targetLayerChan = 7;

chanElGrp = 1:32;
chanTargetLayer = 7; %CA1pyr
chanInterval = 50;%um
chanSpan = chanTargetLayer-8:chanTargetLayer+8;
chanSpan(chanSpan<=0|chanSpan>numel(chanElGrp)) = [];

tHalfWindow = 8;
fitWindows = {-tHalfWindow+1:tHalfWindow,chanSpan-chanTargetLayer}; %samples,chans
fitGrid = cell([1,numel(fitWindows)]);
[fitGrid{:}] = ndgrid(fitWindows{:});


fys = imgaussfilt(log10(sq(ys(:,20,chanSpan))),[5,1]);

figure,
scatter(fitGrid{1}(:),fitGrid{2}(:),40,reshape(fys(t:t+127,:),[],1)-min(reshape(fys(t:t+127,:),[],1)),'filled');



g = fittype( @(A,xa,ya,xya,xo,yo,x,y)                                   ...
             A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
             'independent',{'x', 'y'},...
             'dependent', 'z' ); 

t = 118950;
[fo,fg] = fit([fitGrid{1}(:),fitGrid{2}(:)],...
             reshape(fys(t:t+numel(fitWindows{1})-1,:),[],1) - ...
              min(reshape(fys(t:t+numel(fitWindows{1})-1,:),[],1)),...
             g,...
             'StartPoint',[1,1e-2,1e-2,0,0,0],...
              'Lower',[0.25,1e-5,1e-5,1e-5,-100,-100],...
              'Upper',[5,1e-1,1e-1,1e-1,100,100]);

figure
imagesc(fitWindows{:},...
        fys(t:t+numel(fitWindows{1})-1,:)' - ...
            min(reshape(fys(t:t+numel(fitWindows{1})-1,:),[],1)))
hold('on');
plot(fo.xo,fo.yo,'.')


[rmins] = LocalMinimaN(-fys,-2,10);


figure();
imagesc(1:size(fys,1),fitWindows{2},fys');
colormap('jet');
hold('on');
plot(rmins(:,1),fitWindows{2}(rmins(:,2)),'.');

rmins = sortrows(rmins,1);

fitWeights = exp(0.1.*-sqrt(sum(reshape(cat(ndims(fitGrid)+1,fitGrid{:}),[],ndims(fitGrid)).^2,2)));
fitWeights = fitWeights./sum(fitWeights(:));

figure
scatter(fitGrid{1}(:),fitGrid{2}(:),40,reshape(fitWeights,[],1),'filled');

figure,imagesc(fitWindows{:},fitWeights');

figure();
for r = 1:size(rmins,1),


    clf();
    
    t = rmins(r,1)-tHalfWindow:rmins(r,1)+tHalfWindow-1;
    
    hold('on');
    imagesc(fitWindows{:},fys(t,:)');
    plot(fitWindows{1}(tHalfWindow),fitWindows{2}(rmins(r,2)),'.');
    colorbar();    


    aMax(r)         = fys(rmins(r,1),rmins(r,2));
    aTemporalVar(r) = mean(diff(fys(t,rmins(r,2))).^2);
    aDepthVar(r)    = mean(diff(fys(rmins(r,1),:)).^2);
    aDist(r)        = fitWindows{2}(rmins(r,2));
    
    %    plot_ellipse([0,aDist(r)],5e-4.*sqrt(1./[aTemporalVar(r).^2,abs(aDepthVar(r).^2)]),0.1,0.1);
% $$$     plot_ellipse([fo.xo,fo.yo],sqrt(1./[fo.xa,fo.ya]),0,0.1);    

% $$$     fitWeights = exp(1e-1.*-sqrt(sum(bsxfun(@minus,reshape(cat(ndims(fitGrid)+1,...
% $$$                                                       fitGrid{:}),[],ndims(fitGrid)),...
% $$$                                            [0,fitWindows{2}(rmins(r,2))]).^2,2)));
% $$$     fitWeights = fitWeights./max(fitWeights(:));
% $$$     
% $$$     
% $$$     
% $$$     [fo,fg] = fit([fitGrid{1}(:),fitGrid{2}(:)],                                ...
% $$$                   rSfys.*fitWeights,                                                        ...
% $$$                   g,                                                            ...
% $$$                   'StartPoint',[1,1e-1,1e-1,0,0,fitWindows{2}(rmins(r,2))],     ...
% $$$                   'Lower',[0,1e-2,1e-2,0,-tHalfWindow,-numel(chanSpan)],   ...
% $$$                   'Upper',[1,1,1,1,tHalfWindow,numel(chanSpan)],       ...
% $$$                   'Normalize','on',...
% $$$                   'Robust','off');
% $$$     plot(fo.xo,fo.yo,'.r');
% $$$     plot_ellipse([fo.xo,fo.yo],sqrt(1./[fo.xa,fo.ya]),0,0.1);

    
    waitforbuttonpress();    
end

rmins = sortrows(rmins,1);
rmins = rmins(2:end-1,:);
aPar = nan([size(rmins,1),4]);
%Depth,Peak,svar,tvar

for r = 1:size(rmins,1),
    t = rmins(r,1)-tHalfWindow:rmins(r,1)+tHalfWindow-1;    
    aPar(r,:) = [fitWindows{2}(rmins(r,2)),...
                 fys(rmins(r,1),rmins(r,2)),...
                 mean(diff(fys(rmins(r,1),:)).^2),...
                 mean(diff(fys(t,rmins(r,2))).^2)];
end

figure();
plot3(log10(aPar(:,3)),log10(aPar(:,4)),aPar(:,2),'.');
figure();
plot3(aPar(:,3),uinc(rmins(:,1)),aPar(:,1),'.');


ind = abs(aPar(:,1))<1;
figure();
plot3(log10(aPar(ind,3)),aPar(ind,2),log10(aPar(ind,4)),'.')


%aMax = GetSegs(fys,rmins(:,1)


ufr = Trial.load('ufr',ys,[],units{20},0.05,true,'gauss');
uinc = sum(double(ufr.data>0.5),2);

ind = abs(aPar(:,1))<1;
figure();
plot3(log10(aPar(ind,3)),aPar(ind,2),uinc(rmins(ind,1)),'.')

figure();
plot3(log10(aPar(ind,3)),aPar(ind,2),log10(ys(rmins(ind,1),8)),'.')



figure();
ind = find(abs(aPar(:,1))<1);
plot3(log10(aPar(ind,3)),aPar(ind,2),uinc(rmins(ind,1)),'.')

flfp = lfp.copy();
flfp.filter('ButFilter',4,[6,15],'bandpass');
flfp.resample(ys);


lPar = nan([size(rmins,1),4]);
for r = 1:size(rmins,1),    
    t = rmins(r,1)-tHalfWindow:rmins(r,1)+tHalfWindow-1;
    lPar(r,1) = [sum(flfp(rmins(r,1),chanTargetLayer:15)),...
                 flfp(rmins(r,1),chanTargetLayer),...
                 mean(diff(flfp(rmins(r,1),1:15)).^2),...
                 mean(diff(flfp(t,chanTargetLayer+5)).^2)];
end

figure();
plot3(log10(aPar(ind,3)),aPar(ind,2),log10(lPar(ind,3)),'.')

figure
plot3(log10(aPar(ind,3)),aPar(ind,2),log10(lPar(ind,3)),'.')

[rmins] = LocalMinimaN(-fys,-2.5,10);
rmins = sortrows(rmins,1);
rmins = rmins(2:end-1,:);


figure();
sp = subplot(4,2,[7,8]);        
plot_stc(Trial.stc,ys.sampleRate,[],{'walk','rear','turn','pause','groom','sit'},'brgcmy');
colormap('jet');
for r = 1000:size(rmins,1),
    if ~(abs(fitWindows{2}(rmins(r,2)))<=1),continue,end;
    
    t = (rmins(r,1)-tHalfWindow):(rmins(r,1)+tHalfWindow-1);
    subplot(421);cla();    
        imagesc(flfp(t,:)');
        caxis([-5e3,5e3]);
    subplot(422);cla();
        imagesc(fys(t,:)');
        caxis([2.5,4]);
    subplot(4,2,[3,4]);cla();
        t = (rmins(r,1)-1000):(rmins(r,1)+1000-1);        
        imagesc(t,fitWindows{2},fys(t,:)');
        hold('on');
        plot(rmins(r,1),fitWindows{2}(rmins(r,2)),'.');
        xlim([t([1,end])]);
    subplot(4,2,[5,6]);cla();
        plot(uinc(t))
        ylim([0,15]);
        xlim(sp,[t([1,end])]);    
    waitforbuttonpress();
end


ind = abs(aPar(:,1))<1;

flfps = reshape(permute(GetSegs(flfp.data,rmins(ind,1)-8,16,0),[2,1,3]),sum(ind),[]);
ripss = reshape(permute(GetSegs(fys,rmins(ind,1)-8,16,0),[2,1,3]),sum(ind),[]);


sper = Trial.stc{'s+w+n',ys.sampleRate};
sind = WithinRanges(rmins(ind),sper.data);

% $$$ figure();
% $$$ plot3(log10(aPar(sind,3)),log10(aPar(sind,4)),aPar(sind,2),'.');


[U,S,V] = svd(nunity([flfps,ripss]),0);
[S,V,U,VT] = erpPCA(nunity([flfps(sind,:),ripss(sind,:)]));

figure,
for j = 1:10*6;
tpart = reshape(V(1:size(flfps,2),j),16,32);
subplot(10,6,j);
imagesc(tpart');
end
figure
for j = 1:10*6;
rpart = reshape(V(size(flfps,2)+1:size(flfps,2)+size(ripss,2),j),16,15);
subplot(10,6,j);
imagesc(rpart');
end

sper = Trial.stc{'s',ys.sampleRate};
sind = WithinRanges(rmins(ind),sper.data);


figure();
hold('on');
plot3(U(:,1),U(:,2),U(:,3),'.');
plot3(U(sind,1),U(sind,2),U(sind,3),'.');



%Depth,Peak,svar,tvar
figure();
plot3(U(:,1),U(:,2),(aPar(ind,2)),'.');
figure();
plot3(U(:,1),U(:,2),aPar(ind,3),'.');
figure();
plot3(U(:,1),aPar(ind,2),aPar(ind,3),'.');


sper = Trial.stc{'s+w+n',ys.sampleRate};
sind = WithinRanges(rmins(:,1),sper.data)&abs(aPar(:,1))<1;

aParSub = aPar(sind,:);
lParSub = lPar(sind,:);
figure();
j = 1:size(U,1);
mapped = tsne([U(j,1),U(j,2),aParSub(j,2),log10(aParSub(j,3)),log10(lParSub(j,3))],[],2,4,80);


figure,
plot(aParSub(j,2),log10(aParSub(j,3)),'.');

figure,
plot(mapped(:,1),mapped(:,2),'.');
cpnts = ClusterPP(gcf());

aParSub = aPar(ind,:);
pind = cpnts==2;
figure();
plot3(U(:,2),aPar(ind,2),aPar(ind,3),'.');
hold('on');
plot3(U(pind,2),aParSub(pind,2),aParSub(pind,3),'.','MarkerSize',10);