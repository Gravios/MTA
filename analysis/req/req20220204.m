

mrhm = rhmPhz([stc{'lloc'}]);
mphz = phz([stc{'lloc'}]);

rind =  randsample(size(mphz,1),100000);
figure,hist2([mphz(rind),mrhm(rind)],linspace(0,2*pi,11),linspace(-pi,pi,11),'xprob')

unitSubset = units{18};
meta = sessionList(18);
dc = accumulate_decoding_vars(Trial,                               ...
                              unitSubset,                          ...
                              meta.subject.channelGroup.theta,     ...
                              meta.subject.correction.thetaPhase,  ...
                              meta.subject.correction.headYaw,     ...
                              meta.subject.correction.headBody);
 

rhm = fet_rhm(Trial,sampleRate);
Trial.lfp.filename = [Trial.name,'.lfp'];    
tlfp = Trial.load('lfp',sessionList(tind).subject.channelGroup.theta);
tlfp.resample(rhm);
x = copy(tlfp);
x.data = cat(2,x.data,rhm.data);
mode = 'mtcsdglong';
specArgsTheta = struct('nFFT',2^9,...
                       'Fs',  x.sampleRate,...
                       'WinLength',2^8,...
                       'nOverlap',2^8*0.5,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[1,20]);

[spex,fs,ts ] = fet_spec(Trial,x,mode,false,[],specArgsTheta,[],true);

coh = abs(spex(:,:,1,2))./sqrt(spex(:,:,1,1).*spex(:,:,2,2));
mcoh = copy(spex);
mcoh.data = mean(coh(:,fs>5&fs<12),2);
mcoh.resample(rhm);

mang = copy(spex);
mang.data = circ_mean(angle(spex(:,fs>5&fs<12,1,2)),[],2);
mang.data = unwrap(mang.data);
mang.resample(rhm);
mang.data = mod(mang.data+2*pi,2*pi);
 

derr = nan([size(rhm,1),1]);
derr(dc.ind(dc.ucnt>2)) = sqrt(sum(dc.esax(dc.ucnt>2,:).^2,2));


mind = dc.stcm(:,1)==1 ...
       & (dc.stcm(:,3)==3|dc.stcm(:,5)==5)...
       & dc.ucnt>2 ...
       & abs(dc.hvfl(:,2))<10 ...
       & sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))<350;
clear('xcomp','ycomp','zcomp');
xcomp.label = 'theta phase';      xcomp.data = dc.phz(mind);        xcomp.edgs = linspace(0,2*pi,8); 
% $$$ ycomp.label = 'Rhm Theta Coh';    ycomp.data = mcoh(dc.ind(mind));  ycomp.edgs = linspace(0.2,0.9,21);
%ycomp.label = 'head-body ang';    ycomp.data = dc.hbang(mind);      ycomp.edgs = linspace(-1.2,1.2,8);  
%ycomp.label = 'head pitch';    ycomp.data = dc.fet(mind,1);      ycomp.edgs = linspace(-1.8,0,30);  
ycomp.label = 'velxy';    ycomp.data = dc.lvxy(mind,1);      ycomp.edgs = linspace(0,2,10);  
%ycomp.label = 'rhm phase';    ycomp.data = rhmPhz(dc.ind(mind));      ycomp.edgs = linspace(-pi,pi,8);  
% $$$ ycomp.label = 'rhm-tht dang';     ycomp.data = mang(dc.ind(mind));  ycomp.edgs = linspace(0,2*pi,17);
%vlabel = 'ego lat'; vData = dc.ecom(mind,2); vClim = [-100,100];
vlabel = 'ego fwd'; vData = dc.ecom(mind,1); vClim = [-40,130];
% $$$ vData = mang(dc.ind(mind));
% $$$ vClim = [0,2*pi];
% $$$ vlabel = 'ego';
% $$$ vData = dc.ucnt(mind);
% $$$ vClim = [2,5];
% $$$ vlabel = 'ucount';
[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vData);
zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
set(figure(),'Units','centimeters','Position',[0,-3,8,28]);
subplot(511); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
    cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(512); imagesc(xcomp.ctrs,ycomp.ctrs,zmedian'); axis('xy');
    cax = colorbar(); ylabel(cax,['Median ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(513); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',vlabel]); ylabel(ycomp.label);
    colormap('jet');
    caxis([70,150]);
subplot(514); 
imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
    xlabel(xcomp.label);
subplot(515); 
    imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2,'omitnan'))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
    xlabel(xcomp.label);








anatLoc = 'CA1';
tinds = [17,18,19,20,21,22,23,29];
% $$$ anatLoc = 'CA3';
% $$$ tinds = [6,7,30,26,27,28,29,1,2];
figure
for sts = 1:3;
    cPhi = [];
    cRes = [];
    bind =   sqrt(sum(spkv.ego.^2,2)) < 300 ...
             & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
             & spkv.vxy(:,2)<10 ...
             & spkv.stc(:,1)==1 & any(ismember(spkv.stc(:,stateInds{sts}),stateInds{sts}),2) ...
             & abs(spkv.ego(:,2))       < 150;
    mp = [];
    for t = 1:numel(tinds)
        ind = spkv.map(:,1)==tinds(t) ...
              & bind;
        mp = histcounts(spkv.rhp(ind),rEds);
        mpn = bsxfun(@rdivide,mp,sum(mp));
        if sum(mp)>100
            cPhi(t,1) = angle(sum(mpn.*exp(i*tCtr)));
            cRes(t,1) = abs(sum(mpn.*exp(i*tCtr)));
        else
            cPhi(t,1) = nan;
            cRes(t,1) = nan;
        end
    end
    bind =   sqrt(sum(spkv.ego.^2,2)) < 300 ...
             & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
             & spkv.vxy(:,2)>20 ...
             & spkv.stc(:,1)==1 & any(ismember(spkv.stc(:,stateInds{sts}),stateInds{sts}),2) ...         
             & abs(spkv.ego(:,2))       < 150;
    for t = 1:numel(tinds)
        ind = spkv.map(:,1)==tinds(t) ...
              & bind;
        mp = histcounts(spkv.rhp(ind),rEds);
        mpn = bsxfun(@rdivide,mp,sum(mp));
        if sum(mp)>100        
            cPhi(t,2) = angle(sum(mpn.*exp(i*tCtr)));
            cRes(t,2) = abs(sum(mpn.*exp(i*tCtr)));
        else
            cPhi(t,2) = nan;
            cRes(t,2) = nan;
        end
        
    end
    if sts == 1
        hold('on'),plot(cPhi(:,1)',cRes(:,1)','.r')
    elseif  sts == 2
        hold('on'),plot(cPhi',cRes','.-g')
    elseif  sts == 3
        hold('on'),plot(cPhi',cRes','.-b')
    end
end
xlim([-pi,pi]);
ylim([0,0.3]);


figure,
hist2([spkv.avl(ind),spkv.vxy(ind,2)],linspace(-30,30,20),linspace(0,60,20),'xprob');

figure,
hist2([spkv.avl(ind),log10(spkv.rhw(ind,1))],linspace(-30,30,20),linspace(-4.5,-1.5,20),'xprob');


modes = {'count','prob'};

% $$$ vEds = linspace(-4.5,-1.5,11);
% $$$ vCtr = linspace(-4.5,-1.5,10);
% $$$ vInd = discretize(log10(spkv.rhw(:,1)),vEds);

vEds = linspace(1,60,7);
vCtr = mean([vEds(1:end-1);vEds(2:end)]);
vInd = discretize(spkv.vxy(:,2),vEds);
%vCtr = linspace(-30,30,10);
%vInd = discretize(spkv.avl(:,2),vEds);




mpn = bsxfun(@rdivide,mp,sum(mp));
vlabels = round(vEds);
figure
hold('on');
cmap = cool(size(mpn,2));
for v = 1:size(mpn,1)
    plot(tCtr,mpn(:,v)','Color',cmap(v,:));
    legendLabel{v}  = [num2str(vlabels(v)),'-',num2str(vlabels(v+1)),' cm/s ', ...
                        '\phi: ', num2str(round(angle(sum(mpn(:,v)'.*exp(i*tCtr))),2)),...
                        '  r:', num2str(round(abs(sum(mpn(:,v)'.*exp(i*tCtr))),2))];
end
legend(legendLabel);

mpnCr = abs(sum(mpn(:,3)'.*exp(i*tCtr)));
mpnCm = angle(sum(mpn(:,1)'.*exp(i*tCtr)));
Lines(angle(sum(mpn(:,1)'.*exp(i*tCtr))),[],'r');
Lines(angle(sum(mpn(:,end)'.*exp(i*tCtr))),[],'k');


vEds = linspace(1,60,15);
vCtr = mean([vEds(1:end-1);vEds(2:end)]);
vInd = discretize(spkv.vxy(:,1),vEds);

tEds = linspace(0,2*pi,11);
tCtr = mean([tEds(1:end-1);tEds(2:end)]);

rEds = linspace(-pi,pi,11);
tCtr = mean([tEds(1:end-1);tEds(2:end)]);


mode = modes{2};
stateInds = {4,[5],[7]};
stateLabels = {'Rear','high','low'};
[hfig,fig,fax,sax] = set_figure_layout(figure(6660010),'A4','portrait',[],1.5,1.5,0.05,0.05);
for sts = 1:numel(stateInds)
anatLoc = 'CA1';
ind = ( spkv.map(:,1)==20 ...
       |spkv.map(:,1)==18 ...
       |spkv.map(:,1)==19 ...        
        |spkv.map(:,1)==21 ...
        |spkv.map(:,1)==22 ...
        |spkv.map(:,1)==23 ...
        |spkv.map(:,1)==29) ...
      & spkv.stc(:,1)==1 & any(ismember(spkv.stc(:,stateInds{sts}),stateInds{sts}),2) ...
      & sqrt(sum(spkv.ego.^2,2)) < 300 ...
      & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
      & abs(spkv.ego(:,2)) < 150;
subplot2(5,numel(stateInds),2,sts); 
    mp =[];
    for v = 1:numel(vCtr),
        mp(:,v) = histcounts(spkv.rhp(ind&v==vInd),rEds);
    end
    switch mode
      case 'prob'
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'speed cm/s'});
    end
    xlabel('\phi_{rhm}');
    title(stateLabels{sts});
subplot2(5,numel(stateInds),4,sts);
    mp =[];
    for v = 1:numel(vCtr)
        mp(:,v) = histcounts(spkv.phz(ind&v==vInd),tEds);
    end
    switch mode
      case 'prob'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'speed cm/s'});
    end
    xlabel('\phi_{theta}');
anatLoc = 'CA3';
ind = ( spkv.map(:,1)==6 ...
       |spkv.map(:,1)==7 ...
       |spkv.map(:,1)==30 ...        
        |spkv.map(:,1)==26 ...
        |spkv.map(:,1)==27 ...
        |spkv.map(:,1)==28 ...
        |spkv.map(:,1)==29 ...        
        |spkv.map(:,1)==1 ...        
        |spkv.map(:,1)==2) ...
      & spkv.stc(:,1)==1 & any(ismember(spkv.stc(:,stateInds{sts}),stateInds{sts}),2) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...   
      & abs(spkv.ego(:,2)) < 150;
subplot2(5,numel(stateInds),3,sts);
    mp =[];
    for v = 1:numel(vCtr)
        mp(:,v) = histcounts(spkv.rhp(ind&v==vInd),rEds);
    end
    switch mode
      case 'prob'    
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'speed cm/s'});
    end
    xlabel('\phi_{rhm}');
subplot2(5,numel(stateInds),5,sts);
    mp =[];
    for v = 1:numel(vCtr)
        mp(:,v) = histcounts(spkv.phz(ind&v==vInd),tEds);
    end
    switch mode
      case 'prob'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'speed cm/s'});
    end
    xlabel('\phi_{theta}');
end
cax = colorbar(gca());
sax = gca();
cax.Position(1) = sum(sax.Position([1,3]))+0.01;
pause(0.1);
cax.Position(1) = sum(sax.Position([1,3]))+0.01;
switch mode
  case 'prob'
    ForAllSubplots('caxis([0.02,0.15])');
    ylabel(cax,'Probability');
  case 'count'
    %ForAllSubplots('caxis([0,600])');
    ylabel(cax,'Count');
end



mode = modes{2};
stateInds = {4,[5,6],[7,8]};
stateLabels = {'Rear','high','low'};
[hfig,fig,fax,sax] = set_figure_layout(figure(6660010),'A4','portrait',[],1.5,1.5,0.05,0.05);
for sts = 1:numel(stateInds)
anatLoc = 'CA1 & CA3';
ind = ( spkv.map(:,1)==20 ...
       |spkv.map(:,1)==18 ...
       |spkv.map(:,1)==19 ...        
        |spkv.map(:,1)==21 ...
        |spkv.map(:,1)==22 ...
        |spkv.map(:,1)==23 ...
        |spkv.map(:,1)==29 ...                
        |spkv.map(:,1)==6 ...
        |spkv.map(:,1)==7 ...
        |spkv.map(:,1)==30 ...        
        |spkv.map(:,1)==26 ...
        |spkv.map(:,1)==27 ...
        |spkv.map(:,1)==28 ...
        |spkv.map(:,1)==29 ...        
        |spkv.map(:,1)==1 ...        
        |spkv.map(:,1)==2) ...
      & spkv.stc(:,1)==1 & any(ismember(spkv.stc(:,stateInds{sts}),stateInds{sts}),2) ...
      & sqrt(sum(spkv.ego.^2,2)) < 300 ...
      & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
      & abs(spkv.ego(:,2)) < 150;
subplot2(5,numel(stateInds),2,sts); 
    mp =[];
    for v = 1:numel(vCtr),
        mp(:,v) = histcounts(spkv.rhp(ind&v==vInd),rEds);
    end
    switch mode
      case 'prob'
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(rEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'Nose Speed cm/s'});
    end
    xlabel('\phi_{rhm}');
    title(stateLabels{sts});
    caxis([0.00,0.15]);
subplot2(5,numel(stateInds),3,sts);
    mp =[];
    for v = 1:numel(vCtr)
        mp(:,v) = histcounts(spkv.phz(ind&v==vInd),tEds);
    end
    switch mode
      case 'prob'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   bsxfun(@rdivide,mp,sum(mp))'),...
            'EdgeColor','none');
      case 'count'
        set(pcolor(tEds(1:end-1),                ...
                   vEds(1:end-1),                ...
                   mp'),                         ...
            'EdgeColor','none');
    end
    axis('xy');
    colormap('jet');
    if sts==1
        ylabel({[anatLoc],'Nose Speed cm/s'});
    end
    xlabel('\phi_{theta}');
    if sts==1
        ylabel({[anatLoc],'Nose Speed cm/s'});
    end
    if sts == 3 
        cax = colorbar(gca());
        sax = gca();
        cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        pause(0.1);
        cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        caxis([0.00,0.15]);
    end
        
subplot2(5,numel(stateInds),4,sts);
    set(pcolor(tEds(1:end-1),                ...
               vEds(1:end-1),                ...
               log10(mp)'),                         ...
        'EdgeColor','none');
    if sts == 3 
        cax = colorbar(gca());
        sax = gca();
        cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        pause(0.1);
        cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        ylabel(cax,'Log10 Count');        
    end
    if sts==1
        ylabel({[anatLoc],'Nose Speed cm/s'});
    end
    xlabel('\phi_{theta}');
    if sts==1
        ylabel({[anatLoc],'Nose Speed cm/s'});
    end
    caxis([0,4]);
end




% $$$ 
% $$$     spkv.map   = spkmap;
% $$$     spkv.drz   = spkdrz;
% $$$     spkv.ddz   = spkddz;    
% $$$     spkv.hrz   = spkhrz;
% $$$     spkv.ghz   = spkghz;
% $$$     spkv.pch   = spkpch;
% $$$     spkv.ego   = spkego;
% $$$     spkv.phz   = spkphz;
% $$$     spkv.vxy   = spkvxy;
% $$$     spkv.stc   = spkstc;
% $$$     spkv.avl   = spkavl;
% $$$     spkv.rhp   = spkrhp;    
% $$$     spkv.trans = spktrans;



ind = (spkv.map(:,1)==18 ...
       |spkv.map(:,1)==19 ...        
        |spkv.map(:,1)==21 ...
        |spkv.map(:,1)==22 ...
        |spkv.map(:,1)==23 ...
        |spkv.map(:,1)==24 ...        
        |spkv.map(:,1)==29 ...                
        |spkv.map(:,1)==25) ...
      & spkv.stc(:,1)==1 & (spkv.stc(:,7)==7) ...      
      & sqrt(sum(spkv.ego.^2,2)) < 350 ...
      & spkv.vxy(:,2)<10 ...
      & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
      & abs(spkv.ego(:,2)) < 1000;

figure();
subplot(211);
hist2(spkv.ego(ind&spkv.phz>pi&spkv.rhp<-1&spkv.rhp>-3,:),linspace(-300,400,20),linspace(-300,300,20));
subplot(212);
hist2(spkv.ego(ind&spkv.phz>pi&spkv.rhp>1&spkv.rhp<3,:),linspace(-300,400,20),linspace(-300,300,20));





ind = ( spkv.map(:,1)==3 ...
       |spkv.map(:,1)==4 ...
       |spkv.map(:,1)==5) ...        
      & spkv.stc(:,1)==1 & (spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 300 ...
      & sqrt(sum(spkv.pfc.^2,2)) < 350 ...
      & abs(spkv.ego(:,2)) < 150;

ind = ( spkv.map(:,1)==6 ...
       |spkv.map(:,1)==7 ...
       |spkv.map(:,1)==30 ...        
        |spkv.map(:,1)==26 ...
        |spkv.map(:,1)==27 ...
        |spkv.map(:,1)==28 ...
        |spkv.map(:,1)==1 ...        
        |spkv.map(:,1)==2) ...
       & (spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...      
      & abs(spkv.ego(:,2)) < 150;

% Ed10
ind = ( spkv.map(:,1)==6 ...
       |spkv.map(:,1)==7 ...
        |spkv.map(:,1)==27) ...
       & (spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...      
      & abs(spkv.ego(:,2)) < 150;

% jg05
ind = ( spkv.map(:,1)==30) ...
       & (spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...      
      & abs(spkv.ego(:,2)) < 150;

% ER06
ind = ( spkv.map(:,1)==26) ...
       & (spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...      
      & abs(spkv.ego(:,2)) < 150;
% er01
ind = ( spkv.map(:,1)==28 ...
        |spkv.map(:,1)==1 ...        
        |spkv.map(:,1)==2) ...
       & spkv.stc(:,1)==1 &(spkv.stc(:,7)==7) ...
      & sqrt(sum(spkv.ego.^2,2)) < 250 ...
      & sqrt(sum(spkv.pfc.^2,2))<350 ...      
      & abs(spkv.ego(:,2)) < 150;



cshift = 10000;
cshift = 0;

clear('xcomp','ycomp','zcomp');
xcomp.data = spkv.phz(ind); xcomp.label = 'theta phase';
xcomp.edgs = linspace(0,2*pi,11);
ycomp.data = circshift(spkv.rhp(ind),cshift); ycomp.label = 'rhm phase';
%ycomp.data = circshift(spkv.pcp(ind),cshift); ycomp.label = 'rhm phase';
% $$$ %ycomp.data = spkv.rhp(ind); ycomp.label = 'rhm phase';
ycomp.edgs = linspace(-pi,pi,11);
% $$$ vdata = spkv.pch(ind);
% $$$ vClim = [-1.5,0.5];
% $$$ vlabel = 'head pitch';
vdata = spkv.ego(ind,1);
vClim = [-40,80];
vlabel = 'ego';
% $$$ vdata = abs(spkv.avl(ind,2));
% $$$ vClim = [0,10];
% $$$ vlabel = 'ego';
% $$$ vdata = spkv.ghz(ind);
% $$$ vClim = [-0.2,0.2];
% $$$ vlabel = 'ghz';

[r,p] = circ_corrcc(xcomp.data,ycomp.data)
[r,p] = circ_corrcl(xcomp.data,ycomp.data)


[xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,vdata);
zmean = zcomp.mean;   zmean(zcomp.count<5) = nan;
zmedian = zcomp.median;   zmean(zcomp.count<5) = nan;
zstd  = zcomp.std;    zstd(zcomp.count<5) = nan;
zcount = zcomp.count; zcount(zcomp.count<5) = nan;

figure();
subplot(511); imagesc(xcomp.ctrs,ycomp.ctrs,zmean'); axis('xy');
    cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(512); imagesc(xcomp.ctrs,ycomp.ctrs,zmedian'); axis('xy');
    cax = colorbar(); ylabel(cax,['Median ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(513); imagesc(xcomp.ctrs,ycomp.ctrs,zstd'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',vlabel]); ylabel(ycomp.label);
    colormap('jet');
subplot(514); 
imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
%imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2))'); axis('xy');    
%imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,1))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
    xlabel(xcomp.label);
subplot(515); 
%imagesc(xcomp.ctrs,ycomp.ctrs,zcount'); axis('xy');
imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2))'); axis('xy');    
%imagesc(xcomp.ctrs,ycomp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,1))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
    xlabel(xcomp.label);
    %caxis([0.125,0.165])
% $$$ caxis([0.18,0.22])
% $$$ caxis([0.22,0.27])
caxis([0.08,0.12])

figure()
hold('on');
plot(angle(sum(bsxfun(@times,bsxfun(@rdivide,zcount,sum(zcount,2)),exp(i.*ycomp.ctrs)'))),...
     abs(sum(bsxfun(@times,bsxfun(@rdivide,zcount,sum(zcount,2)),exp(i.*ycomp.ctrs)'))),'-+r');
xlim([-pi,pi])

figure()
polarplot(angle(sum(bsxfun(@times,bsxfun(@rdivide,zcount,sum(zcount,2)),exp(i.*ycomp.ctrs)'))),...
     abs(sum(bsxfun(@times,bsxfun(@rdivide,zcount,sum(zcount,2)),exp(i.*ycomp.ctrs)'))),'-+g');
hold('on');


rthresh = 0.7;
figure,
subplot(311);
rose(circ_dist(spkv.phz(ind&abs(spkv.rhp)>rthresh),spkv.rhp(ind&abs(spkv.rhp)>rthresh)),31)
subplot(312);
rose(circ_dist(spkv.phz(ind&abs(spkv.rhp)<rthresh),spkv.rhp(ind&abs(spkv.rhp)<rthresh)),31)
subplot(313);
rose(circ_dist(circshift(spkv.phz(ind&abs(spkv.rhp)<rthresh),1000),spkv.rhp(ind&abs(spkv.rhp)<rthresh)),31)


configure_default_args();
MjgER2016_load_data();


tind = 20;
Trial = Trials{tind};
sampleRate = 250;
stc = copy(Trial.stc);
xyz = preproc_xyz(Trial,'trb',sampleRate);
vxy = vel(filter(copy(xyz),'ButFilter',4,2.4,'low'),{'spine_lower','nose','hcom'},[1,2]);
rhm = fet_rhm(Trial,sampleRate);
rhmPhz = phase(rhm,[5,12]);
rhmPhz.data = circshift(rhmPhz.data,-250);

% $$$ ang = create(MTADang,Trial,xyz);
% $$$ figure,plot(ang(:,'spine_upper','hcom',3));
% $$$ dang = copy(xyz);
% $$$ dang.data = ang(:,'spine_upper','hcom',3);
% $$$ dang.filter('ButFilter',4,[4,12],'bandpass');
% $$$ figure,
% $$$ plot(nunity(diff(dang.data)))
% $$$ hold('on');
% $$$ plot(nunity(rhm.data)+4,'m')
% $$$ 
% $$$ ddang = copy(dang);
% $$$ ddang.data = circshift(dang.data,-1)-circshift(dang.data,1);
% $$$ ddang.filter('ButFilter',4,[4,12],'bandpass');
% $$$ ddang.data = circshift(ddang.data,-1)-circshift(ddang.data,1);
% $$$ rhmPhz = phase(ddang,[5,12]);

rhmn = fet_rhm_hcom(Trial,sampleRate);
figure,plot(bsxfun(@plus,rhmn.data,[-2,0,2]))

rhmn.data = [circshift(rhmn(:,1),-1)-circshift(rhmn(:,1),1),...
             circshift(rhmn(:,2),-1)-circshift(rhmn(:,2),1),...
             circshift(rhmn(:,3),-1)-circshift(rhmn(:,3),1)];
rhmn.filter('ButFilter',4,[4,12],'bandpass');
rhmn.data = [circshift(rhmn(:,1),-1)-circshift(rhmn(:,1),1),...
             circshift(rhmn(:,2),-1)-circshift(rhmn(:,2),1),...
             circshift(rhmn(:,3),-1)-circshift(rhmn(:,3),1)];


figure,
plot(bsxfun(@plus,nunity(rhmn.data),[-6,0,6]));
hold('on');
plot(nunity(rhm.data)+3,'m');
plot(nunity(pchn.data)+1,'k');

rhmn.data = rhmn(:,1);

rhmnPhz = phase(rhmn,[5,12]);


figure,plot(circ_dist(rhmnPhz.data,rhmPhz.data))



hvxy = fet_href_HXY(Trial);
hvxy.data = circshift(hvxy(:,1),1)-circshift(hvxy(:,1),-1);
hvxy.resample(sampleRate);
fhvxy = filter(copy(hvxy),'ButFilter',4,[5,14],'bandpass');
rhmPhz = phase(fhvxy,[4,12]);

pch = fet_HB_pitchB(Trial,sampleRate);
pch.data = pch.data(:,1);
pchn = copy(pch);
pchn.filter('ButFilter',4,[3,14],'bandpass');
pchn.data = circshift(pchn(:,1),-1)-circshift(pchn(:,1),1);
pchn.filter('ButFilter',4,[3,14],'bandpass');
pchn.data = circshift(pchn(:,1),-1)-circshift(pchn(:,1),1);
%pchPhz = phase(pch,[5,12]);
pchnPhz = phase(pchn,[5,12]);
rhmPhz = phase(pchn,[5,12]);

figure,
plot(nunity(rhm([stc{'lloc'}])));
plot(nunity(pchn([stc{'lloc'}])),'m')

figure,rose(circ_dist(pchnPhz([stc{'lloc'}]),rhmPhz([stc{'lloc'}])));

xcomp.data = spkv.phz(ind); xcomp.label = 'theta phase';
xcomp.edgs = linspace(0,2*pi,6);
ycomp.data = circshift(spkv.rhp(ind),cshift); ycomp.label = 'rhm phase';
ycomp.edgs = linspace(-pi,pi,6);
vdata = spkv.ego(ind,1);
vClim = [-30,80];
vlabel = 'ego';

clear('xcp','ycp','zcp');
xcp.edgs = linspace(-pi,pi,30);
xcp.data = circ_dist(pchnPhz([stc{'lloc'}]),rhmPhz([stc{'lloc'}])-0.4);
xcp.label = 'rhmpch pdiff';
ycp.edgs = linspace(-1.4,-0.6,30);
ycp.data = pch([stc{'lloc'}],1);
ycp.label = 'hpitch';
vdata = vxy([stc{'lloc'}],2);
vClim = [0,30];
vlabel = 'nspeed';

[xcp,ycp,zcp] = compute_2d_discrete_stats(xcp,ycp,vdata);
zmean = zcp.mean;   zmean(zcp.count<5) = nan;
zmedian = zcp.median;   zmean(zcp.count<5) = nan;
zstd  = zcp.std;    zstd(zcp.count<5) = nan;
zcount = zcp.count; zcount(zcp.count<5) = nan;

figure();
subplot(511); imagesc(xcp.ctrs,ycp.ctrs,zmean'); axis('xy');
    cax = colorbar(); ylabel(cax,['Mean ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(512); imagesc(xcp.ctrs,ycp.ctrs,zmedian'); axis('xy');
    cax = colorbar(); ylabel(cax,['Median ',vlabel]); colormap('jet');  caxis([vClim]);
subplot(513); imagesc(xcp.ctrs,ycp.ctrs,zstd'); axis('xy');
    cax = colorbar(); ylabel(cax,['Std ',vlabel]); ylabel(ycp.label);
    colormap('jet');
subplot(514); 
imagesc(xcp.ctrs,ycp.ctrs,zcount'); axis('xy');
    cax = colorbar(); ylabel(cax,'Count'); colormap('jet'); 
    xlabel(xcp.label);
subplot(515); 
imagesc(xcp.ctrs,ycp.ctrs,bsxfun(@rdivide,zcount,sum(zcount,2))'); axis('xy');    
    cax = colorbar(); ylabel(cax,'Probability'); colormap('jet'); 
    xlabel(xcp.label);





phz = load_theta_phase(Trial, ...
                       rhm, ...
                       sessionList(tind).subject.channelGroup.theta,...
                       sessionList(tind).subject.correction.thetaPhase);
lfp = Trial.load('lfp',65:96);

sind = logical(get(resample(cast([stc{'lloc&theta'}],'TimeSeries'),xyz),'data'));



flfp = copy(lfp);
flfp.filter('ButFilter',4,[0.8,20],'bandpass');
% $$$ csdl = copy(lfp);
% $$$ csdl.data = csdl.data(:,1:27);
% $$$ for c = 1:27;
% $$$     csdl.data(:,c) = (flfp(:,c)+flfp(:,c+5)-2.*flfp(:,c+3))./100;
% $$$ end

csdl = copy(lfp);
csdl.data = csdl.data(:,1:29);
for c = 1:29;
    csdl.data(:,c) = (flfp(:,c)+flfp(:,c+3)-2.*flfp(:,c+1))./100;
end

csdl = copy(lfp);
csdl.data = csdl.data(:,1:29);
for c = 1:29;
    csdl.data(:,c) = (flfp(:,c)-flfp(:,c+3))./100;
end

rcsdl = resample(copy(csdl),rhm);



sind = logical(get(resample(cast([stc{states{3}}],'TimeSeries'),rhm),'data'));
sind = logical(get(resample(cast([stc{states{5}}],'TimeSeries'),rhm),'data'));
%sind = sind & abs(circ_dist(pchnPhz(:),rhmPhz(:)-0.4))<1;

tshifts = fliplr(-25:25);
mrcsdl = [];
srcsdl = [];
tts = 26
% $$$ for tts = 1:numel(tshifts)
    rind = sind &vxy(:,2)<12&vxy(:,2)>2 ;    
    rind = circshift(rind,tshifts(tts));
    for pr = 1:nbins,
        for pt = 1:tbins,
        ind = phzInds==pt & rhmPhzInds==pr & rind;
        ind(ind) = 1==(rem(1:sum(ind),2));
        ind(ind) = randn([sum(ind),1])>0;
        mrcsdl(tts,:,pr,pt) = mean(rcsdl(ind,:));
        srcsdl(tts,:,pr,pt) = std(rcsdl(ind,:));
        end
    end

% $$$ for tts = 1:numel(tshifts)
    rind = sind &vxy(:,2)<12&vxy(:,2)>2 ;    
    rind = circshift(rind,tshifts(tts));

phzTstat = [];    
phzDf = [];    
    for pt = 1:tbins,
        indR = phzInds==pt & rhmPhzInds==3 & rind;
        indR(indR) = 1==(rem(1:sum(indR),2));
        %indR(indR) = randn([sum(indR),1])>0;
        indL = phzInds==pt & rhmPhzInds==6 & rind;
        indL(indL) = 1==(rem(1:sum(indL),2));
        %indL(indL) = randn([sum(indL),1])>0;        
        phzTstat(pt,:) = (mean(rcsdl(indR,:))-mean(rcsdl(indL,:)))./(sqrt(((sum(indR)-1).*std(rcsdl(indR,:)).^2+(sum(indL)-1).*std(rcsdl(indL,:)).^2)./(sum(indL)+sum(indR)-2)).*sqrt(1/sum(indL)+1/sum(indR)));
        phzDf(pt) = sum(indL)+sum(indR)-2;
    end
% $$$ figure,imagesc(phzTstat')
% $$$ colormap('jet');
% $$$ caxis([3,3.1])

phzPval = [];
for c = 1:size(phzTstat,2)
    for pt = 1:tbins
        phzPval(pt,c) = 1-tcdf(abs(phzTstat(pt,c)),phzDf(pt));
    end
end
figure,imagesc(phzPval'); caxis([0,0.001]);
colormap('jet');

figure,
hold('on');

c = 13;
mc1 = [];
mc2 = [];
for pt = 1:tbins;
    disp(pt);
for iter = 1:1000
pr = 4;
ind = phzInds==pt & rhmPhzInds==pr & rind;
ind(ind) = 1==circshift(rem(1:sum(ind),2),double(randn([1,1])>0));
ind(ind) = randn([sum(ind),1])>0;
mc1(iter,pt) = mean(rcsdl(ind,c));
%histogram(rcsdl(ind,c),linspace(-100,100,30),'EdgeColor','none','FaceAlpha',0.5)    
pr = 8;
ind = phzInds==pt & rhmPhzInds==pr & rind;
ind(ind) = 1==circshift(rem(1:sum(ind),2),double(randn([1,1])>0));
ind(ind) = randn([sum(ind),1])>0;
mc2(iter,pt) = mean(rcsdl(ind,c));
end
end
%histogram(rcsdl(ind,c),linspace(-100,100,30),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','r')    

figure();
for pt = 1:tbins
    subplot(tbins,1,pt);
    hold('on');    
    histogram(mc1(:,pt),linspace(-40,40,50),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','c')    
    histogram(mc2(:,pt),linspace(-40,40,50),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','r')    
    Lines(prctile(mc1(:,pt),95),[],'c');
    Lines(prctile(mc2(:,pt),5),[],'r');
    Lines(prctile(mc1(:,pt),5),[],'c');
    Lines(prctile(mc2(:,pt),95),[],'r');
end


figure();hold('on');
errorbar([phzCntr,phzCntr+2*pi],[mean(mc1),mean(mc1)],repmat(prctile(mc1(:,:),95)-mean(mc1),[1,2]),'c')
errorbar([phzCntr,phzCntr+2*pi],[mean(mc2),mean(mc2)],repmat(prctile(mc2(:,:),95)-mean(mc2),[1,2]),'r')
errorbar([phzCntr,phzCntr+2*pi]+0.1,[mean(mc1),mean(mc1)],repmat(sq(srcsdl(tts,c,4,:))',[1,2]),'c')
errorbar([phzCntr,phzCntr+2*pi]+0.1,[mean(mc2),mean(mc2)],repmat(sq(srcsdl(tts,c,8,:))',[1,2]),'r')
Lines(2*pi,[],'k')

% $$$ 
% $$$ tshifts = fliplr(-25:25);
% $$$ mvcsdl = [];
% $$$ tts = 26
% $$$ % $$$ for tts = 1:numel(tshifts)
% $$$     rind = sind &vxy(:,2)<12&vxy(:,2)>2 ;    
% $$$     rind = circshift(rind,tshifts(tts));
% $$$     for pr = 1:nbins,
% $$$         for pt = 1:tbins,
% $$$         ind = phzInds==pt & velInds==pr & rind;
% $$$         mvcsdl(tts,:,pr,pt) = mean(rcsdl(ind,:));
% $$$         end
% $$$     end
% $$$ 
% $$$ figure;
% $$$ for p = 1:nbins,
% $$$     subplot(1,nbins,p);
% $$$     imagesc(sq(mvcsdl(tts,:,p,:))-sq(mrcsdl(tts,:,5,:)));
% $$$     caxis([-10,10])
% $$$     %imagesc(sq(mvcsdl(tts,:,p,:))-sq(mrcsdl(tts,:,7,:)));
% $$$     %caxis([-45,45])
% $$$     colormap('jet');
% $$$ end
% $$$ 
% $$$ figure;
% $$$ for p = 1:nbins,
% $$$     subplot(1,nbins,p);
% $$$ % $$$     imagesc(sq(mrcsdl(tts,:,p,:))-sq(mrcsdl(tts,:,10,:)));
% $$$ % $$$     caxis([-10,10])
% $$$     imagesc(sq(mrcsdl(tts,:,p,:)));
% $$$     caxis([-45,45])
% $$$     colormap('jet');
% $$$ end

% $$$ end

tshifts = fliplr(-25:25);
mrcsdlShuff = [];
srcsdlShuff = [];
% $$$ for tts = 1:numel(tshifts)
    disp(tts)
    tic
    rind = sind & vxy(:,2)<12 &vxy(:,2)>2;
    rind = circshift(rind,tshifts(tts));
    rhmPhzIndsRperm = rhmPhzInds(rind);
    thtPhzIndsRperm = phzInds(rind);
    rcsdlRperm      = rcsdl(rind,:);
    for iter = 1:100    
        rhmPhzIndsRperm = rhmPhzIndsRperm(randperm(size(rhmPhzIndsRperm,1)));                
        for pr = 1:nbins,
            for pt = 1:tbins,
                ind = thtPhzIndsRperm==pt & rhmPhzIndsRperm==pr;
                ind(ind) = randn([sum(ind),1])>0;
                mrcsdlShuff(tts,:,pr,pt,iter) = mean(rcsdlRperm(ind,:));
                srcsdlShuff(tts,:,pr,pt,iter) = std(rcsdlRperm(ind,:));
            end
        end
    end
    toc
% $$$ end

    
figure,
for p = 1:nbins
    sax = subplot2(4,nbins,1,p);
    imagesc(sq(mrcsdl(tts,:,p,:)));
    caxis([-40,40]);
    if p==nbins
        cax = colorbar();
        cax.Position(1) = sum(sax.Position([1,3]));
    end
    sax = subplot2(4,nbins,2,p);        
    imagesc(sq(srcsdl(tts,:,p,:)));
    if p==nbins
        cax = colorbar();
        cax.Position(1) = sum(sax.Position([1,3]));
    end
    caxis([0,35]);
    sax = subplot2(4,nbins,3,p);        
    %imagesc(sq(mean(srcsdlShuff(tts,:,p,:,:),5)));
    imagesc(sq(std(srcsdlShuff(tts,:,p,:,:),[],5)));    
    if p==nbins
        cax = colorbar();
        cax.Position(1) = sum(sax.Position([1,3]));
    end
    caxis([0,5]);
    sax = subplot2(4,nbins,4,p);        
    imagesc(sq(srcsdl(tts,:,p,:))-sq(mean(srcsdlShuff(tts,:,p,:,:),5)));
    if p==nbins
        cax = colorbar();
        cax.Position(1) = sum(sax.Position([1,3]));
    end
    caxis([-5,0]);
    
end
colormap('jet');




% $$$ mrcsdlv = [];
% $$$ for pr = 1:16,
% $$$     for pt = 1:16,
% $$$         for tts = 1:numel(tshifts)
% $$$         ind = phzInds==pt & rhmPhzInds==pr & sind &vxy(:,2)>10 ;
% $$$         ind = circshift(ind,tshifts(tts));
% $$$         mrcsdlv(tts,:,pr,pt) = mean(rcsdl(ind,:));
% $$$         end
% $$$     end
% $$$ end

% $$$ figure,imagesc(rcsdl(ind,:)')
% $$$ figure,plot(mean(rcsdl(ind,:)))

figure();
for s = 1:16
    subplot(2,8,s);
    imagesc(sq(mrcsdl(:,3,s,:))');
    caxis([-20,20])
end
colormap('jet')


figure();
for s = 1:25
    subplot(5,5,s);
    imagesc(sq(mrcsdl(20,s,:,:))');
    caxis([-40,40])
end
colormap('jet')

figure;
tpi = 8;
%rpi = 12;
chan = 11;
for rpi = 1:16;
subplot2(2,16,1,rpi);
im = repmat(sq(mrcsdl(:,chan,rpi,:)),[3,3]);
im = imgaussfilt(im,3);
im = im(52:102,17:33);
imagesc(im');
Lines(26,[],'k');
subplot2(2,16,2,rpi);
im = repmat(sq(mrcsdlv(:,chan,rpi,:)),[3,3]);
im = imgaussfilt(im,3);
im = im(52:102,17:33);
imagesc(im');
Lines(26,[],'k');
colormap('jet')
end

figure();
for s = 1:25
    subplot(5,5,s);
    imagesc(sq(mrcsdlv(45,s,:,:))');
    caxis([-30,30])
end
colormap('jet')


figure();
for s = 1:25
    subplot(5,5,s);
    imagesc(sq(mrcsdl(26,s,:,:))');
end

figure();
for t = 1:10;
imo = repmat(sq(mrcsdl(t+20,:,13,:)),[1,3]);
imo = imgaussfilt(imo,1);
imo = imo(:,17:33);
for rpi = 1:16,
subplot2(10,16,t,rpi);
im = repmat(sq(mrcsdl(t+20,:,rpi,:)),[1,3]);
im = imgaussfilt(im,1);
im = im(:,17:33);
imagesc(im-imo);
title(num2str(t+20));
end
end
colormap jet
ForAllSubplots('caxis([-15,15])')
ForAllSubplots('caxis([3,10])')
ForAllSubplots('caxis([-10,-3])')


figure,plot(sq(mrcsdl(26,19,:,9)));


    


figure();
for t = 1:10;
imo = repmat(sq(mrcsdl(t+20,:,13,:)),[1,3]);
imo = imgaussfilt(imo,1);
imo = imo(:,17:32);
for rpi = 1:16,
subplot2(10,16,t,rpi);
im = repmat(sq(mrcsdl(t+20,:,rpi,:)),[1,3]);
im = imgaussfilt(im,1);
im = im(:,17:32);
imagesc(im-imo);
title(num2str(t+20));
end
end
colormap jet
ForAllSubplots('caxis([-15,15])')
ForAllSubplots('caxis([3,10])')
ForAllSubplots('caxis([-10,-3])')



t = 26;
p = 14;
ims = zeros(size(mrcsdl,2),16,16,100);
for iter = 1:100;
imo = repmat(sq(mrcsdlShuff(t,:,p,:,iter)),[1,3]);
imo = imgaussfilt(imo,1);
imo = imo(:,17:32);
for rpi = 1:16
    imt = repmat(sq(mrcsdlShuff(t,:,rpi,:,iter)),[1,3]);
    imt = imgaussfilt(imt,1);
    imt = imt(:,17:32);
    ims(:,rpi,:,iter) = imt-imo;
end
end
im = zeros(size(mrcsdl,2),16,16);
imo = repmat(sq(mrcsdl(t,:,p,:)),[1,3]);
imo = imgaussfilt(imo,1);
imo = imo(:,17:32);
for rpi = 1:16
    imt = repmat(sq(mrcsdl(t,:,rpi,:)),[1,3]);
    imt = imgaussfilt(imt,1);
    imt = imt(:,17:32);
    im(:,rpi,:) = imt-imo;
end
figure,
for p = 1:16
    subplot2(4,16,1,p);
    imagesc(sq((im(:,p,:)-mean(ims(:,p,:,:),4))./std(ims(:,p,:,:),[],4)))
    colormap('jet');
    caxis([-10,10]);
    %caxis([2,5]);    
    %caxis([-5,-2]);    
    subplot2(4,16,2,p);
    imagesc(sq((im(:,p,:))))
    colormap('jet');
    caxis([-15,15]);   
    subplot2(4,16,3,p);
    imagesc(sq(mrcsdl(26,:,p,:)));
    colormap('jet');
    caxis([-35,35]);       
    subplot2(4,16,4,p);
    imagesc(sq(mean(mrcsdlShuff(26,:,p,:,:),5)));
    colormap('jet');
    caxis([-35,35]);           
end


t = 26;
imd = zeros(size(mrcsdl,2),nbins,tbins,100);
for rpi = 1:nbins
    imo = repmat(sq(mrcsdl(t,:,rpi,:)),[1,3]);
    imo = imgaussfilt(imo,1);
    imo = imo(:,tbins+1:tbins*2);
    for iter = 1:100;
        imt = repmat(sq(mrcsdlShuff(t,:,rpi,:,iter)),[1,3]);
        imt = imgaussfilt(imt,1);
        imt = imt(:,tbins+1:tbins*2);
        imd(:,rpi,:,iter) = imo-imt;
    end
end

figure,
for p = 1:nbins
    subplot2(3,nbins,1,p);
    imagesc(sq(mean(imd(:,p,:,:),4)));
    caxis([-5,5]);
    colormap(gca(),'jet');    
    subplot2(3,nbins,2,p);
    imagesc(sq(std(imd(:,p,:,:),[],4)));
    caxis([0,1]);
    colormap(gca(),'jet');
    subplot2(3,nbins,3,p);    
    imagesc(abs(sq(mean(imd(:,p,:,:),4))./sq(std(imd(:,p,:,:),[],4))));    
    caxis([3,5]);    
    colormap(gca(),'parula');    
end


    



t = 26;
p = 13;
imd = zeros(size(mrcsdl,2),nbins,nbins,100);
for rpi = 1:nbins
    for iter = 1:100;
        imo = repmat(sq(mrcsdlShuff(t,:,rpi,:,101-iter)),[1,3]);
        imo = imgaussfilt(imo,1);
        imo = imo(:,17:32);
        imt = repmat(sq(mrcsdlShuff(t,:,rpi,:,iter)),[1,3]);
        imt = imgaussfilt(imt,1);
        imt = imt(:,17:32);
        imd(:,rpi,:,iter) = imo-imt;
    end
end

figure,
for p = 1:nbins
    subplot2(1,nbins,1,p);
    imagesc(sq(mean(imd(:,p,:,:),4)));
end
colormap('jet');
ForAllSubplots('caxis([-5,5])')

ForAllSubplots('caxis([-3,3])')
ForAllSubplots('caxis([3,4])')
ForAllSubplots('caxis([-4,-3])')
ForAllSubplots('caxis([-10,10])')

figure,
hold('on')
lc = jet(27);
plot(sq(mean(mean(imd(10,:,:,:),4),3)))
plot(sq(mean(mean(imd(16,:,:,:),4),3)))
plot(sq(mean(mean(imd(19,:,:,:),4),3)))
plot(sq(mean(mean(imd(25,:,:,:),4),3)))

figure,
hold('on')
lc = jet(27);
for c = 1:27
plot(sq(mean(mean(imd(c,:,:,:),4),3)),'Color',lc(c,:))
end