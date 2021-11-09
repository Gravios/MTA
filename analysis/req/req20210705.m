% req20210705
% 
% XYHB 30sps 300ms

global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

sampleRate = 30;

%%%<<< XYHB Decoding -------------------------------------------------------------------------------

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

tind = [3:5,17:25];
cind = 1;
dc = {};
for t = tind;
    Trial = Trials{t};
    tunits = units{t};
    thetaRefChan = sessionList(t).thetaRefGeneral;

    halfSpkWindow = 0.15;
    ufrWindow = 0.15;

    smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];

    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    pfs = compute_xyhb_ratemaps( Trial, tunits);

    ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
    mask = ds.mask;

    spk = Trial.load('spk', sampleRate, '', tunits, 'deburst');
    ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

    dc{cind} = decode_ufr_boxcar(Trial, ...
                                 tunits, ...
                                 sampleRate, ...
                                 ufr, ...
                                 pfs, ...
                                 mask,...
                                 halfSpkWindow,...
                                 smoothingWeights,...
                                 'overwrite',true);

    dc{cind}.pfs = pfs;

    dc{cind}.xyz  = copy(xyz);
    dc{cind}.xyz.data  = xyz.data(dc{cind}.ind,:,:);

    dc{cind}.hRot = 0.17;
    dc{cind}.vxy  = filter(copy(xyz),'ButFilter',4,1.5,'low');;
    dc{cind}.vxy  = dc{cind}.vxy.vel('hcom',[1,2]);
    dc{cind}.vxy.data  = dc{cind}.vxy(dc{cind}.ind,1);
    dc{cind}.lvxy = copy(dc{cind}.vxy);
    dc{cind}.lvxy.data(dc{cind}.lvxy.data<=0) = 0.001;
    dc{cind}.lvxy.data = log10(dc{cind}.lvxy.data);

    dc{cind}.hvec = dc{cind}.xyz(:,'nose',[1,2])-dc{cind}.xyz(:,'hcom',[1,2]);
    dc{cind}.hvec = sq(bsxfun(@rdivide,dc{cind}.hvec,sqrt(sum(dc{cind}.hvec.^2,3))));
    dc{cind}.hvec = cat(3,dc{cind}.hvec,sq(dc{cind}.hvec)*[0,-1;1,0]);
    dc{cind}.hvec = multiprod(dc{cind}.hvec,                           ...
                              [cos(dc{cind}.hRot),-sin(dc{cind}.hRot); ...
                        sin(dc{cind}.hRot), cos(dc{cind}.hRot)],...
                              [2,3],                                   ...
                              [1,2]);

    dc{cind}.tvec = circshift(xyz(:,'hcom',[1,2]),-round(250*0.1)) - ...
        circshift(xyz(:,'hcom',[1,2]),round(250*0.1));
    dc{cind}.tvec = sq(bsxfun(@rdivide,dc{cind}.tvec,sqrt(sum(dc{cind}.tvec.^2,3))));
    dc{cind}.tvec = cat(3,dc{cind}.tvec,sq(dc{cind}.tvec)*[0,-1;1,0]);
    dc{cind}.tvec = dc{cind}.tvec(dc{cind}.ind,:,:);

    dc{cind}.fet  = fet_HB_pitchB(Trial,dc{cind}.sampleRate);
    dc{cind}.fet.data  = dc{cind}.fet.data(dc{cind}.ind,:);

    dc{cind}.hdist = filter(copy(xyz),'ButFilter',3,5,'low');
    xycoor =    xyz(:,'hcom',[1,2]);
    [~,dc{cind}.hdist.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hdist.data = dc{cind}.hdist.data(dc{cind}.ind);

    % COM 
    dc{cind}.ecom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.com},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.ecom = dc{cind}.ecom{1};
    % SAX 
    dc{cind}.esax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.sax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.esax = dc{cind}.esax{1};
    % MAX 
    dc{cind}.emax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.max},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.emax = dc{cind}.emax{1};

    % HVANG
    dc{cind}.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
    xycoor = cat(2,...
                 dc{cind}.hvang(:,'spine_upper',[1,2])-dc{cind}.hvang(:,'bcom',[1,2]),...
                 dc{cind}.hvang(:,'nose',[1,2])-dc{cind}.hvang(:,'hcom',[1,2]));
    dc{cind}.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    % Positive: CCW (Left)     Negative: CW (Right)
    dc{cind}.hvang.data = circ_dist(circshift(dc{cind}.hvang.data(:,2),-10),...
                                    circshift(dc{cind}.hvang.data(:,2),+10));
    dc{cind}.hvang = dc{cind}.hvang.data(dc{cind}.ind);

    
    
% COMPUTE corrected hbang
    dc{cind}.hbang = filter(copy(xyz),'ButFilter',3,9,'low');    
    xycoor = cat(2,...
                 dc{cind}.hbang(:,'spine_upper',[1,2])-dc{cind}.hbang(:,'bcom',[1,2]),...
                 dc{cind}.hbang(:,'nose',[1,2])-dc{cind}.hbang(:,'hcom',[1,2]));
    dc{cind}.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hbang.data = circ_dist(dc{cind}.hbang.data(:,2),dc{cind}.hbang.data(:,1));
    dc{cind}.hbang = dc{cind}.hbang.data(dc{cind}.ind)+0.2+-0.4*double(cind>4);

% GENERATE state collection matrix    
    dc{cind}.stcm = stc2mat(Trial.stc,xyz,states);
    dc{cind}.stcm = dc{cind}.stcm(dc{cind}.ind,:);

% COMPUTE head and body relative velocity vectors
    hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    dc{cind}.hvfl = hvfl(dc{cind}.ind,:);
    dc{cind}.bvfl = bvfl(dc{cind}.ind,:);    
    
    cind  = cind +1;
end                     

%%%>>>


%%%<<< LOCATE optimal ego centric position for decoding XY -----------------------------------------

xshifts = -80:5:80;
yshifts = -80:5:80;
medianPosError = {};
mmpe = {};
for cind = 1:numel(dc),
    ind = dc{cind}.stcm(:,1)==1 & dc{cind}.stcm(:,2)~=2  & dc{cind}.ucnt>=2;
    medianPosError{cind} = nan([numel(xshifts),numel(yshifts)]);
    for xshift = 1:numel(xshifts),
        for yshift = 1:numel(yshifts),
            medianPosError{cind}(xshift,yshift) = median(sqrt(sum([dc{cind}.ecom(ind,1)+xshifts(xshift),...
                                                                   dc{cind}.ecom(ind,2)+yshifts(yshift)].^2,...
                                                                   2)...
                                                              )...
                                                         );
        end
    end
    mmpe{cind} = LocalMinimaN(medianPosError{cind},1000,10);    
end


figure,
for cind = 1:numel(dc)
    subplot(3,4,cind);
    set(pcolor(yshifts,xshifts,medianPosError{cind}),'EdgeColor','none');
    hold('on');
    scatter(yshifts(mmpe{cind}(2)),xshifts(mmpe{cind}(1)),30,'r','Filled');
end

%%%>>>


%%%<<< PLOT timeseries XYHB ------------------------------------------------------------------------

cind = 7;
figure,
ind = dc{cind}.stcm(:,1)==1  & dc{cind}.ucnt>=3;
subplot(411);
    hold('on');
    plot(dc{cind}.ind(ind),dc{cind}.com(ind,1),'.r');
    plot(dc{cind}.ind(ind),dc{cind}.sax(ind,1),'.g');    
    plot(dc{cind}.ind(ind),dc{cind}.xyz(ind,'hcom',1),'.b');
subplot(412);
    hold('on');
    plot(dc{cind}.ind(ind),dc{cind}.com(ind,2),'.r');
    plot(dc{cind}.ind(ind),dc{cind}.sax(ind,2),'.g');        
    plot(dc{cind}.ind(ind),dc{cind}.xyz(ind,'hcom',2),'.b');
subplot(413);
    hold('on');
    plot(dc{cind}.ind(ind),dc{cind}.com(ind,3),'.r');
    plot(dc{cind}.ind(ind),dc{cind}.sax(ind,3),'.g');            
    plot(dc{cind}.ind(ind),dc{cind}.fet(ind,1),'.b');
subplot(414);
    hold('on');
    plot(dc{cind}.ind(ind),dc{cind}.com(ind,4),'.r');
    plot(dc{cind}.ind(ind),dc{cind}.sax(ind,4),'.g');                
    plot(dc{cind}.ind(ind),dc{cind}.fet(ind,2),'.b');
    %linkaxes(findobj([figure(22),figure(23)],'Type','Axes'),'x');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

xlim([27770,28474]);

figure,
hist(abs(dc{cind}.ecom(ind,4)),100)
Lines(prctile(abs(dc{cind}.ecom(ind,4)),95),[],'r');

%%%>>>


%%%<<< CONVERT cell to array -----------------------------------------------------------------------
cind = 1;
dca = dc{cind};
dca.sessionId = tind(cind).*ones(size(dc{cind}.ind));

%dca.xyz = dca.xyz.data;
dca.pos = sq(dca.xyz(:,'hcom',:));
dca.fet = dca.fet.data;
dca.vxy = dca.vxy.data;
dca.lvxy = dca.lvxy.data;
dca.hdist = dca.hdist.data;

dca.uinc = dc{1}.uinc;
maxUinc = max(cell2array(cf(@(d) size(d.uinc,2),dc)))
dca.uinc = cat(2,dc{1}.uinc,nan(abs([[0,maxUinc]-size(dc{1}.uinc)])));

for cind = 2:numel(dc),
    dca.sessionId = cat(1,dca.sessionId,tind(cind).*ones(size(dc{cind}.ind)));
    dca.ind = cat(1,dca.ind,dc{cind}.ind);
    dca.max = cat(1,dca.max,dc{cind}.max);
    dca.com = cat(1,dca.com,dc{cind}.com);    
    dca.sax = cat(1,dca.sax,dc{cind}.sax);        
    dca.post = cat(1,dca.post,dc{cind}.post);
    dca.ucnt = cat(1,dca.ucnt,dc{cind}.ucnt);    
    dca.uinc = cat(1,dca.uinc,cat(2,dc{cind}.uinc,nan(abs([[0,maxUinc]-size(dc{cind}.uinc)]))));    
    %dca.xyz = cat(1,dca.xyz,dc{cind}.xyz.data);    
    dca.pos = cat(1,dca.pos,sq(dc{cind}.xyz(:,'hcom',:)));
    dca.vxy = cat(1,dca.vxy,dc{cind}.vxy.data);        
    dca.lvxy = cat(1,dca.lvxy,dc{cind}.lvxy.data);            
    dca.hvec = cat(1,dca.hvec,dc{cind}.hvec);    
    dca.tvec = cat(1,dca.tvec,dc{cind}.tvec);        
    dca.fet = cat(1,dca.fet,dc{cind}.fet.data);
    dca.hvang = cat(1,dca.hvang,dc{cind}.hvang);    
    dca.hbang = cat(1,dca.hbang,dc{cind}.hbang);        
    dca.hdist = cat(1,dca.hdist,dc{cind}.hdist.data);        
    
    dca.ecom  = cat(1,dca.ecom,dc{cind}.ecom);        
    dca.esax  = cat(1,dca.esax,dc{cind}.esax);        
    dca.emax  = cat(1,dca.emax,dc{cind}.emax);
    
    dca.stcm = cat(1,dca.stcm,dc{cind}.stcm);    
    dca.hvfl = cat(1,dca.hvfl,dc{cind}.hvfl);        
    dca.bvfl = cat(1,dca.bvfl,dc{cind}.bvfl);            
    
    dca.hRot  = cat(1,dca.hRot,dc{cind}.hRot);
    dca.pfs   = cat(1,dca.pfs,dc{cind}.pfs);
end


dca.posi = discretize(dca.pos,linspace(-500,500,21));
dca.feti(:,1) = discretize(dca.fet(:,1),linspace(-2,0.8,29));
dca.feti(:,2) = discretize(dca.fet(:,2),linspace(-0.8,2,29));


%%%>>>


%%%<<< COMPUTE bhv error as a function of spatial position -----------------------------------------

posBins = [-500:50:500];
posBinc = [posBins(2:end)+posBins(1:end-1)]/2;

dca.posi = discretize(dca.pos,posBins);

ind = dca.stcm(:,1)==1  &  dca.ucnt>=2 & nniz(dca.posi) & ismember(dca.sessionId,[3:5,17:23]);

bhvHpErrStdXY = accumarray(dca.posi(ind,1:2),dca.esax(ind,3),numel(posBins)-[1,1],@std);
bhvHpErrStdXY = bhvHpErrStdXY.*mask(:,:,20,20);
bhvHpErrStdXY(bhvHpErrStdXY==0) = nan;

bhvBpErrStdXY = accumarray(dca.posi(ind,1:2),dca.esax(ind,4),numel(posBins)-[1,1],@std);
bhvBpErrStdXY = bhvBpErrStdXY.*mask(:,:,20,20);
bhvBpErrStdXY(bhvBpErrStdXY==0) = nan;

figure();
subplot(121);
set(pcolor(posBinc,posBinc,bhvHpErrStdXY'),'EdgeColor','none');
subplot(122);
set(pcolor(posBinc,posBinc,bhvBpErrStdXY'),'EdgeColor','none');

%%%>>>


%%%<<< COMPUTE the egocentric shift to match head position -----------------------------------------
ind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2  & dca.ucnt>=2;
xshifts = -80:5:80;
yshifts = -80:5:80;
medianPosErrorAll = nan([numel(xshifts),numel(yshifts)]);
for xshift = 1:numel(xshifts),
    for yshift = 1:numel(yshifts),
        medianPosErrorAll(xshift,yshift) = median(sqrt(sum([dca.ecom(ind,1)+xshifts(xshift),...
                                                         dca.ecom(ind,2)+yshifts(yshift)].^2,2)));
    end
end
mmpeAll = LocalMinimaN(medianPosErrorAll,300,10);
figure,
set(pcolor(yshifts,xshifts,medianPosErrorAll),'EdgeColor','none');
hold('on');
scatter(yshifts(mmpeAll(2)),xshifts(mmpeAll(1)),30,'r','Filled');
colorbar();

%%%>>>


ind = dca.stcm(:,1)==1 & dca.ucnt>=2;
ind = dca.stcm(:,1)==1 & dca.stcm(:,4)==4 & dca.ucnt>=2;
ind = dca.stcm(:,1)==1 & dca.stcm(:,3)==3 | dca.stcm(:,5)==5 & dca.ucnt>=2;
figure();
bhvError = sqrt(sum([dca.com(ind,3)-dca.fet(ind,1),dca.com(ind,4)-dca.fet(ind,2)].^2,2));
histogram(bhvError,linspace(0,2,50));
Lines(median(bhvError),[],'r');
hold('on')
bhvErrorShuff = sqrt(sum([dca.com(ind,3)-circshift(dca.fet(ind,1),-15000),...
                          dca.com(ind,4)-circshift(dca.fet(ind,2),-15000)].^2,2));
histogram(bhvErrorShuff,linspace(0,2,50));
Lines(median(bhvErrorShuff),[],'r');


% STD of error head pitch
figure();
ind =     dca.stcm(:,1)==1 ...
       & (dca.stcm(:,5)==5 | dca.stcm(:,4)==4 | dca.stcm(:,3)==3 | dca.stcm(:,6)==6 |dca.stcm(:,2)==2) ...
       &  dca.ucnt>=2;
shifts = -120:120;
pstd = [];
for shft = 1:numel(shifts)
    pstd(end+1) = std([dca.com(ind,3)-circshift(dca.fet(ind,1),shifts(shft))]);
    %pstd(end+1) = std([dca.fet(ind,1)-circshift(dca.fet(ind,1),shifts(shft))]);
end
subplot(211);
plot(shifts,pstd);
title('STD of shifted Head pitch error');
grid('on');
% STD of error Body pitch
pstd = [];
for shft = 1:numel(shifts)
    %pstd(end+1) = std([dca.fet(ind,2)-circshift(dca.fet(ind,2),shifts(shft))]);
    pstd(end+1) = std([dca.com(ind,4)-circshift(dca.fet(ind,2),shifts(shft))]);
end
subplot(212);
plot(shifts,pstd);
title('STD of shifted Body pitch error');
grid('on');




ind = dca.stcm(:,1)==1 & dca.stcm(:,3)==3 | dca.stcm(:,5)==5 & dca.ucnt>=2;
shft = 5;
figure();
subplot(211);
bhvError = [dca.com(ind,3)-circshift(dca.fet(ind,1),shft),dca.com(ind,4)-circshift(dca.fet(ind,2),shft)];
histogram2(bhvError(:,1),bhvError(:,2),linspace(-2,2,30),linspace(-2,2,30),'DisplayStyle','tile');



ind = dca.stcm(:,1)==1 & dca.ucnt>=2;
ind = dca.stcm(:,1)==1 & dca.stcm(:,3)==3 | dca.stcm(:,5)==5 & dca.ucnt>=2;
figure();
subplot(211);
bhvError = [dca.com(ind,3)-dca.fet(ind,1),dca.com(ind,4)-dca.fet(ind,2)];
histogram2(bhvError(:,1),bhvError(:,2),linspace(-2,2,30),linspace(-2,2,30),'DisplayStyle','tile');

subplot(212);
bhvErrorShuff = [dca.com(ind,3)-circshift(dca.fet(ind,1),-5000),...
                 dca.com(ind,4)-circshift(dca.fet(ind,2),-5000)];
histogram2(bhvErrorShuff(:,1),...
           bhvErrorShuff(:,2),...
           linspace(-2,2,30),...
           linspace(-2,2,30),...
           'DisplayStyle','tile');




ind = dca.stcm(:,1)==1 & dca.ucnt>=2;
ind = dca.stcm(:,1)==1 & dca.stcm(:,3)==3 | dca.stcm(:,5)==5 & dca.ucnt>=2;
figure();
hold('on');
bhvError = sqrt(sum([dca.com(ind,3)-dca.fet(ind,1)].^2,2));
histogram(bhvError,linspace(0,2,60));
Lines(median(bhvError),[],'r');
bhvErrorShuff = sqrt(sum([dca.com(ind,3)-circshift(dca.fet(ind,1),-5000)].^2,2));
histogram(bhvErrorShuff,linspace(0,2,60));
Lines(median(bhvErrorShuff),[],'r');



ind = dca.stcm(:,1)==1 & dca.ucnt>=2;
ind = dca.stcm(:,1)==1 & (dca.stcm(:,3)==3 | dca.stcm(:,5)==5) & dca.ucnt>=2;
figure();
hold('on');
bhvError = [dca.com(ind,3)-dca.fet(ind,1)];
histogram(bhvError,linspace(-2,2,30));
Lines(median(bhvError),[],'r');
bhvErrorShuff = [dca.com(ind,3)-circshift(dca.fet(ind,1),-5000)];
histogram(bhvErrorShuff,linspace(-2,2,30));
Lines(median(bhvErrorShuff),[],'r');


%%%<<< Behavior decoded vs obsevered ---------------------------------------------------------------

hfig = set_figure_layout(figure(),'A4','portrait');
corrHP = [];
corrBP = [];
nSamples = 1e3;
subplotCount = 1;
for uc = 1:15;
sind = dca.stcm(:,1)==1 & dca.ucnt==uc;
hSamples = randsample(find(sind & (dca.stcm(:,3)==3|dca.stcm(:,4)==4)),nSamples);
vecHSamples = false(size(dca.ind));
vecHSamples(hSamples) = true;

lSamples = randsample(find(sind & (dca.stcm(:,5)==5|dca.stcm(:,6)==6)),nSamples);
vecLSamples = false(size(dca.ind));
vecLSamples(lSamples) = true;

rSamples = randsample(find(sind & dca.stcm(:,2)==2),nSamples);
vecRSamples = false(size(dca.ind));
vecRSamples(rSamples) = true;

if mod(uc,5)==0,
    subplot2(3,3,1,subplotCount);
    plot(dca.fet(ind,1),dca.com(ind,3),'.');
    title({'Head Pitch',['Unit Count: ',num2str(uc)]});
    xlabel('Observed (rad)');
    ylabel('Decoded (rad)');
    xlim([-2,1]);
    ylim([-2,1]);
    daspect(gca(),[1,1,1]);    
    
    subplot2(3,3,2,subplotCount);
    plot(dca.fet(ind,2),dca.com(ind,4),'.');
    title({'Body Pitch Pitch',['Unit Count: ',num2str(uc)]});
    xlabel('Observed (rad)');
    ylabel('Decoded (rad)');
    xlim([-0.25,1.5]);
    ylim([-0.25,1.5]);
    daspect(gca(),[1,1,1]);
    
    subplotCount = subplotCount + 1;
end    

ind =  (vecHSamples | vecLSamples | vecRSamples) ;
corrHP(uc) = corr(dca.com(ind,3),dca.fet(ind,1));
corrBP(uc) = corr(dca.com(ind,4),dca.fet(ind,2));
end
subplot2(3,3,3,[1:3]);
hold('on');
plot(1:15,corrHP,'-+');
plot(1:15,corrBP,'-+');
ylabel('Correlation');
xlabel('Unit Count');
ylim(gca(),[0.25,0.8]);
xlim(gca(),[0,16]);
box(gca(),'on');
legend(gca(),{'Head','Body'},'Location','NorthWest');
title('Correlation of decoded vs observed pitches varying with coactive Units');

%%%>>>


%%%<<< STUFF ---------------------------------------------------------------------------------------

figure();
subplot(311);
bhvError = [dca.com(ind,4)-dca.fet(ind,2)];
histogram(bhvError,linspace(-2,2,30));
Lines(std(bhvError),[],'r');
Lines(-std(bhvError),[],'r');
subplot(312);
bhvErrorShuff = [dca.com(ind,4)-circshift(dca.fet(ind,2),-15000)];
histogram(bhvErrorShuff,linspace(-2,2,30));
Lines(std(bhvErrorShuff),[],'r');
Lines(-std(bhvErrorShuff),[],'r');
subplot(313);
bhvErrorShuffAuto = [dca.fet(ind,2)-circshift(dca.fet(ind,2),-15000)];
histogram(bhvErrorShuffAuto,linspace(-2,2,30));
Lines(std(bhvErrorShuffAuto),[],'r');
Lines(-std(bhvErrorShuffAuto),[],'r');
linkaxes(findobj(gcf(),'Type','Axes'),'y');


ind = dca.stcm(:,1)==1 & (dca.stcm(:,5)==5) & dca.ucnt>=2;
figure();
hold('on');
bhvError = [dca.com(ind,4)-dca.fet(ind,2)];
set(histogram(bhvError,linspace(-2,2,30)),'FaceColor','g');
Lines(median(bhvError),[],'r');
bhvErrorShuff = [dca.com(ind,4)-circshift(dca.fet(ind,2),-15000)];
set(histogram(bhvErrorShuff,linspace(-2,2,30)),'FaceColor','c');
Lines(median(bhvErrorShuff),[],'r');
bhvErrorShuffAuto = [dca.fet(ind,2)-circshift(dca.fet(ind,2),-15000)];
set(histogram(bhvErrorShuffAuto,linspace(-2,2,30)),'FaceColor','r');
Lines(median(bhvErrorShuffAuto),[],'r');




ind = dca.stcm(:,1)==1 & (dca.stcm(:,3)==3 | dca.stcm(:,5)==5) & dca.ucnt>=2;
figure();
hold('on');
bhvError = [dca.com(ind,3)-dca.fet(ind,1)];
set(histogram(bhvError,linspace(-2,2,30)),'FaceColor','g');
Lines(median(bhvError),[],'r');
std(bhvError)
bhvErrorShuff = [dca.com(ind,3)-circshift(dca.fet(ind,1),-15000)];
set(histogram(bhvErrorShuff,linspace(-2,2,30)),'FaceColor','c');
Lines(median(bhvErrorShuff),[],'r');
std(bhvErrorShuff)
bhvErrorShuffAuto = [dca.fet(ind,1)-circshift(dca.fet(ind,1),-15000)];
set(histogram(bhvErrorShuffAuto,linspace(-2,2,30)),'FaceColor','r');
Lines(median(bhvErrorShuffAuto),[],'r');
std(bhvErrorShuffAuto)



ind = dca.stcm(:,1)==1 & (dca.stcm(:,3)==3 | dca.stcm(:,5)==5 | dca.stcm(:,2)==2 ) & dca.ucnt>=2;
figure();
hold('on');
bhvError = [dca.com(ind,4)-dca.fet(ind,2)];
histogram(bhvError,linspace(-2,2,30));
Lines(median(bhvError),[],'r');
bhvErrorShuff = [dca.com(ind,4)-circshift(dca.fet(ind,2),-5000)];
histogram(bhvErrorShuff,linspace(-2,2,30));
Lines(median(bhvErrorShuff),[],'r');





cind = 7;
stc = copy(Trials{tind(cind)}.stc);
hfig = set_figure_layout(figure(),'A4','landscape');
ind = dc{cind}.stcm(:,1)==1;
subplot2(5,1,[1,2],1);
    hold('on');
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.ecom(ind,3),'.r');
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.com(ind,3),'.c');            
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.fet(ind,1),'.b');
    grid('on');
    ylim([-2,2]);
    Lines([],0,'k');
    ylabel('Head Pitch');    
    legend(gca(),{'Error','Decoded','Observed'});
subplot2(5,1,[3,4],1);
    hold('on');
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.ecom(ind,4),'.r');
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.com(ind,4),'.c');                
    plot(dc{cind}.ind(ind)./sampleRate,dc{cind}.fet(ind,2),'.b');
    grid('on');
    ylim([-2,2]);
    Lines([],0,'k');    
    ylabel('Body Pitch');
    legend(gca(),{'Error','Decoded','Observed'});    
subplot2(5,1,[5],1);
    plotSTC(stc,1,'text');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');
xlim([310,340]);
suptitle('Decoded Behavior: Timeseries')
    
%%%>>>


%%%<<<  Theta(dsc,asc) : hvfl vs ego forward, hbang vs ego lateral ---------------------------------

figure();
norm = 'xprob';
%norm = '';
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3;
sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & dca.ucnt>=2;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & dca.sessionId==20;
subplot(221);
    ind = ismember(dca.iphz,[3:10])  & sind;
    hist2([dca.esax(ind,1),...
           dca.bvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(222);
    ind = ismember(dca.iphz,[16:23]) & sind;
    hist2([dca.esax(ind,1),...
           dca.bvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(223);
    ind = ismember(dca.iphz,[3:10]) & sind;
    hist2([dca.esax(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('ego lateral');
    ylabel('hba rad');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(224);
    ind = ismember(dca.iphz,[16:23]) & sind;
    hist2([dca.esax(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('ego lateral');
    ylabel('hba rad');
    Lines([],0,'k');
    Lines(0,[],'k');

%%%>>>




%%%<<< MAIN FIGURE decoding ------------------------------------------------------------------------
t = 20; 
exyz = preproc_xyz(Trials{t},'trb');
exyz.resample(dc{t==tind}.sampleRate);
efet =  fet_HB_pitchB(Trials{t},dc{t==tind}.sampleRate);
ets = [1:size(exyz,1)]./sampleRate;
maskcirc = create_tensor_mask(dc{1}.pfs.adata.bins(1:2));
maskcirc = sq(mask(:,:,14,14));
maskbhv = sq(mask(11,11,:,:));

eind =   dc{t==tind}.ucnt >= 2;

%exampleRange = [3332,3417];
exampleRange = [8950,10100];
%exampleRange = [27770,28474];
exampleTimestamps = (dc{t==tind}.ind-1)/dc{t==tind}.sampleRate;
%eTS = [290,345];
%eTS = [530,600];
%eTS = [1000,1045];
eInds = find(WithinRanges(dc{t==tind}.ind,exampleRange));
exampleSubsetRange = [9101,9200];
exampleSubsetRange = [9650,9950];

eIndsSub = find( WithinRanges( dc{t==tind}.ind, exampleSubsetRange));

exampleTimestampsSubset = (dc{t==tind}.ind(eIndsSub)-1)./dc{t==tind}.sampleRate;

enanmaskSubset = double( dc{t==tind}.ucnt(eIndsSub) >= 3     ...
                 & dc{t==tind}.stcm(eIndsSub,1)==1);
enanmaskSubset(~enanmaskSubset)=nan;

enanmask = double( dc{t==tind}.ucnt(eInds) >= 3     ...
                 & dc{t==tind}.stcm(eInds,1)==1);
enanmask(~enanmask)=nan;



[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],1.5,1.5,0,0.2);

%%%<<< PLOT example trajectory within physical space

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height-fig.subplot.verticalPadding, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.verticalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc{t==tind}.pfs.adata.bins{1:2},~maskcirc');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

circle(0,0,480,'-k');
plot(exyz(exampleSubsetRange,'hcom',1),...
     exyz(exampleSubsetRange,'hcom',2),...
     '-k','LineWidth',1);
plot(sq(dc{t==tind}.com(eIndsSub,1)).*enanmaskSubset,...
     sq(dc{t==tind}.com(eIndsSub,2)).*enanmaskSubset,...
     '-c','LineWidth',1);
xlim([-500,500]);
ylim([-500,500]);
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');


%%%>>>

%%%<<< PLOT example trajectory within behavior space
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, -fig.subplot.height-fig.subplot.verticalPadding, 1, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+fig.subplot.verticalPadding,  ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');
imagesc(dc{1}.pfs.adata.bins{3:4},~maskbhv');
colormap(sax(end),'gray')
caxis(sax(end),[-2,1])

axis('xy');
plot(efet(exampleSubsetRange,1),...
     efet(exampleSubsetRange,2),...
     '-k','LineWidth',1);
plot(sq(dc{t==tind}.com(eIndsSub,3)).*enanmaskSubset,...
     sq(dc{t==tind}.com(eIndsSub,4)).*enanmaskSubset,...
     '-c','LineWidth',1);
xlim(sax(end),[-1.8,0.6]);
ylim(sax(end),[-0.6,1.8]);
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');

%%%>>>

%%%<<< PLOT timeseries of xcoord vs decoded xcoord

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     exyz(exampleRange,'hcom',1),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.com(eInds,1)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.com(eIndsSub,1)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,400,'X');
%%%>>>

%%%<<< PLOT timeseries of ycoord vs decoded ycoord
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(2, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     exyz(exampleRange,'hcom',2),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.com(eInds,2)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.com(eIndsSub,2)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,400,'Y');
%%%>>>

%%%<<< PLOT timeseries of head pitch vs decoded head pitch
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     efet(exampleRange,1),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.com(eInds,3)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.com(eIndsSub,3)).*enanmaskSubset,...
     '-c','LineWidth',1);
sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
ylim([-1.6,0.5]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,0.2,'HP');
%%%>>>

%%%<<< PLOT timeseries of body pitch vs decoded body pitch

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(4, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plot(ets(exampleRange(1):exampleRange(2)),...
     efet(exampleRange,2),...
     '-k','LineWidth',1);
plot(exampleTimestamps(eInds),...
     sq(dc{t==tind}.com(eInds,4)).*enanmask,...
     '-r','LineWidth',1);
plot(exampleTimestampsSubset,...
     sq(dc{t==tind}.com(eIndsSub,4)).*enanmaskSubset,...
     '-c','LineWidth',1);

sax(end).XTickLabels = [];
sax(end).YTick = [];
xlim(ets([exampleRange(1),exampleRange(2)]));
ylim([-0.5,1.6]);
box(sax(end),'on');
text(ets(exampleRange(1))+0.25,1.3,'BP');

%%%>>>

%%%<<< PLOT timeseries of states

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(5, 0, 3, 0.4);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                ...
                  'Position',[fig.page.xpos(xind)+xOffSet,              ...
                              fig.page.ypos(yind)+yOffSet,              ...
                              fig.subplot.width*7,                      ...
                              fig.subplot.height],                      ...
                  'FontSize', 8,                                        ...
                  'LineWidth',1);
hold('on');
plotSTC(Trials{t}.stc,1,'text',states([1,fliplr(2:6)]),'kbbggr');
xlim(ets([exampleRange(1),exampleRange(2)]));
box(sax(end),'on');


%%%>>>

%%%<<< PLOT HB error conditioned on XY

% ADJUST subplot coordinates

[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height-fig.subplot.verticalPadding, 10, 0.6);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.horizontalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.posi(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
        & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.posi(ind,1:2),sqrt(sum(dca.ecom(ind,3:4).^2,2)),[20,20],@mean);
eCount = accumarray(dca.posi(ind,1:2),ones([sum(ind),1]),[20,20],@sum);
out(eCount<100) = nan;
circle(0,0,480,'-k');

set(pcolor(dca.pfs(1).adata.bins{1}-25,...
           dca.pfs(1).adata.bins{2}-25,out'),'EdgeColor','none');
axis('xy');
xlim([-500,500]);
ylim([-500,500]);
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
cax = colorbar(sax(end),'eastoutside');
cax.Units = 'Centimeters';
cax.Position(1) = fig.page.xpos(xind)+xOffSet+0.1 + fig.subplot.width*2+ fig.subplot.verticalPadding;
cax.Label.String = 'rad';
%%%>>>

%%%<<< PLOT XY error conditioned on HB

% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3,-fig.subplot.height-fig.subplot.verticalPadding, 10, 0.6);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.horizontalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.feti(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
      & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.feti(ind,1:2),sqrt(sum(bsxfun(@plus,dca.ecom(ind,1:2),[-35,0]).^2,2)),[28,28],@mean);
eCount = accumarray(dca.feti(ind,1:2),ones([sum(ind),1]),[28,28],@sum);
out(eCount<100) = nan;

set(pcolor(dca.pfs(1).adata.bins{3},...
           dca.pfs(1).adata.bins{4},out'/10),'EdgeColor','none');
axis('xy');
xlim([-1.6,0.5]);
ylim([-0.5,1.6]);
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
cax = colorbar(sax(end),'eastoutside');
cax.Units = 'Centimeters';
cax.Position(1) = fig.page.xpos(xind)+xOffSet+0.1 + fig.subplot.width*2+ fig.subplot.verticalPadding;
cax.Label.String = 'cm';
%%%>>>

%%%<<< PLOT XY error as function of included units
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(1,-fig.subplot.height-fig.subplot.verticalPadding, 14, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.horizontalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);
hist2([dca.ucnt(ind),sqrt(sum(bsxfun(@plus,dca.esax(ind,1:2), [-35,0]).^2,2))/10],...
      linspace(0.5,20.5,21),...
      linspace(0,80,41),'');

%%%>>>

%%%<<< PLOT HB error as function of included units
% ADJUST subplot coordinates
[yind, yOffSet, xind, xOffSet] = deal(3,-fig.subplot.height-fig.subplot.verticalPadding, 14, 0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.horizontalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);
hist2([dca.ucnt(ind),sqrt(sum(dca.esax(ind,3:4).^2,2))],...
      linspace(0.5,20.5,21),...
      linspace(0,1.5,41),'');

%%%>>>

% ENDFIG

%%%<<< PLOT MUD XY vs HB error

bfsn = cf(@(t,u)  compute_bhv_ratemaps(t,u),          Trials, units);
bfss = cf(@(t,u)  compute_bhv_ratemaps_shuffled(t,u), Trials, units);

[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfsn, units, [], [], false);
numComp = size(eigVecs,2);
bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfsn)';


ngBins = sum(validDims);

bhvInfoN = {};
for t = 1:numel(bfsn)
for u = 1:numel(units{t})
rmap = bfsn{t}.data.rateMap(validDims,bfsn{t}.data.clu==units{t}(u));
MRate = mean(rmap,'omitnan');
bhvInfoN{t}(u) = sum(1/ngBins.*(rmap./MRate).*log2(rmap./MRate),'omitnan');
end
end

bhvInfoS = {};
for t = 1:numel(bfsn)
for u = 1:numel(units{t})
rmap = sq(bfss{t}.data.rateMap(validDims,bfsn{t}.data.clu==units{t}(u),:));
MRate = mean(rmap,'omitnan');
bhvInfoS{t}(u,:) = sum(1/ngBins.*(rmap./MRate).*log2(rmap./MRate),'omitnan');
end
end


% $$$ u = 60;
% $$$ (bhvInfoN{20}(66)-mean(bhvInfoS{20}(66,:),2))./std(bhvInfoS{20}(66,:),[],2)


for ind = 1:size(dca.uinc,1)
    nzd = nonzeros(bhvInfoN{dca.sessionId(ind)}.*dca.uinc(ind,~isnan(dca.uinc(ind,:))));
    dca.BhvInfoM(ind,1) = mean(nzd);
    dca.BhvInfoS(ind,1) = std(nzd);    
end

ind = dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2 | dca.stcm(:,3)==3 | dca.stcm(:,4)==4 | dca.stcm(:,5)==5 | dca.stcm(:,6)==6);

[uPosErr,xPosErr] = MakeUniformDistr(sqrt(sum(bsxfun(@plus,dca.esax(ind,1:2), [-35,0]).^2,2))/10);
[uBhvErr,xBhvErr] = MakeUniformDistr(sqrt(sum(dca.esax(ind,3:4).^2,2)));

uPosErrInd = discretize(uPosErr,linspace([xPosErr([1,end])',31]));
uBhvErrInd = discretize(uBhvErr,linspace([xBhvErr([1,end])',31]));

figure();
hist2([uPosErr,uBhvErr],linspace(0,89,31),linspace(0,3.66,31));

ucntErr = dca.ucnt(ind);


nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ucntErr);
out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ucntErr(nind),[30,30],@mean);

ubhvInfoM = dca.BhvInfoM(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoM);
out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoM(nind),[30,30],@mean);

ubhvInfoS = dca.BhvInfoS(ind);
nind = nniz(uPosErrInd) & nniz(uBhvErrInd) & nniz(ubhvInfoS);
out = accumarray([uPosErrInd(nind),uBhvErrInd(nind)],ubhvInfoS(nind),[30,30],@mean);

xticks = xPosErr(round([1:size(xPosErr,1)/5:size(xPosErr,1),size(xPosErr,1)]));
yticks = xBhvErr(round([1:size(xBhvErr,1)/5:size(xBhvErr,1),size(xBhvErr,1)]));

figure,
imagesc(out')
set(gca(),'XTick',0:6:30);
set(gca(),'YTick',0:6:30);
set(gca(),'XTickLabel',xticks);
set(gca(),'YTickLabel',yticks);
axis('xy');
colormap('jet');

%%%>>>
%%%<<< PLOT unit count Condtiond on XY,HB error mud
%%%>>>




%%%>>>

%%%<<< SUPFIG conditional error for each session ---------------------------------------------------

[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'A4','landscape',[],3,3,0.5,0.5);


%%%<<< PLOT HB error conditioned on XY

for cind = 1:numel(tind)
% ADJUST subplot coordinates


[yind, yOffSet, xind, xOffSet] = ...
    deal(floor(cind/4)+1,-fig.subplot.height-fig.subplot.verticalPadding, mod(cind,4)+double(~mod(cind,4)), 0.6);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width, ...
                              fig.subplot.height],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.posi(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
        & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.posi(ind,1:2),sqrt(sum(dca.ecom(ind,3:4).^2,2)),[20,20],@mean);
eCount = accumarray(dca.posi(ind,1:2),ones([sum(ind),1]),[20,20],@sum);
out(eCount<100) = nan;
circle(0,0,480,'-k');

set(pcolor(dca.pfs(1).adata.bins{1}-25,...
           dca.pfs(1).adata.bins{2}-25,out'),'EdgeColor','none');
axis('xy');
xlim([-500,500]);
ylim([-500,500]);
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
% $$$ cax = colorbar(sax(end),'eastoutside');
% $$$ cax.Units = 'Centimeters';
% $$$ cax.Position(1) = fig.page.xpos(xind)+xOffSet+0.1 + fig.subplot.width*2+ fig.subplot.verticalPadding;
% $$$ cax.Label.String = 'rad';
%%%>>>
end

%%%<<< PLOT XY error conditioned on HB

% ADJUST subplot coordinates

[yind, yOffSet, xind, xOffSet] = deal(3,-fig.subplot.height-fig.subplot.verticalPadding, 10, 0.6);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                              fig.page.ypos(yind)+yOffSet,                      ...
                              fig.subplot.width*2+ fig.subplot.horizontalPadding, ...
                              fig.subplot.height*2+fig.subplot.verticalPadding],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

ind =   nniz(dca.feti(:,1:2)) ...
      & dca.stcm(:,1)==1 ...
      & (dca.stcm(:,2)==2|dca.stcm(:,3)==3|dca.stcm(:,4)==4|dca.stcm(:,5)==5|dca.stcm(:,6)==6) ...
      & dca.ucnt >= 2;
ind(ind) = logical(randn([sum(ind),1])>0.65);
out = accumarray(dca.feti(ind,1:2),sqrt(sum(bsxfun(@plus,dca.ecom(ind,1:2),[-35,0]).^2,2)),[28,28],@mean);
eCount = accumarray(dca.feti(ind,1:2),ones([sum(ind),1]),[28,28],@sum);
out(eCount<100) = nan;

set(pcolor(dca.pfs(1).adata.bins{3},...
           dca.pfs(1).adata.bins{4},out'/10),'EdgeColor','none');
axis('xy');
xlim([-1.6,0.5]);
ylim([-0.5,1.6]);
colormap(sax(end),'jet');
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
cax = colorbar(sax(end),'eastoutside');
cax.Units = 'Centimeters';
cax.Position(1) = fig.page.xpos(xind)+xOffSet+0.1 + fig.subplot.width*2+ fig.subplot.verticalPadding;
cax.Label.String = 'cm';
%%%>>>

%%%>>>
