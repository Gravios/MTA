
MTA_PROJECT_REPORT_PATH = '/storage/share/Projects/BehaviorPlaceCode/';
figurePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP','temporalShift');
create_directory(figurePath);

vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below
el = 'FL';

%%%<<< HBA x HVF
% $$$ % 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% $$$ vlb{end+1} = 'hba';
% $$$ vdc{end+1} = 'head-body angle';
% $$$ vnm{end+1} = 'cha';
% $$$ vbe{end+1} = [-1.2,-0.6,-0.2,0.2,0.6,1.2];
% $$$ vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
% $$$ vun{end+1} = 'rad';

% $$$ vlb{end+1} = 'hva';
% $$$ vdc{end+1} = 'head-body anglular veloctiy';
% $$$ vnm{end+1} = 'chvang';
% $$$ vbe{end+1} = [ -0.3,-0.17,-0.05,0.05,0.17,0.3];
% $$$ vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
% $$$ vun{end+1} = 'rad';

% $$$ vlb{end+1} = 'hvl';
% $$$ vdc{end+1} = 'lateral head speed';
% $$$ %vnm{end+1} = 'cfhrvl';
% $$$ vnm{end+1} = 'dfhrvl';
% $$$ vbe{end+1} = [-40,-24,-16,-8,-3,3,8,16,24,40];
% $$$ %nvbe{end+1} = [-40,-4,4,40];
% $$$ vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
% $$$ vun{end+1} = 'cm/s';

%2
vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
%vbe{end+1} = [-4,4,8,16,24,32,40,48,56,64,80];
vbe{end+1} = [-4,4,10,20,30,40,50,60,85];
%vbe{end+1} = [-4,4,14,28,42,56,70];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>


mirror = false;
stgrps = {[10]};
stlbls = {'Walk'};
% $$$ stgrps = {[9]};
% $$$ stlbls = {'Turn'};
% $$$ stgrps = {11};
% $$$ stlbls = {'Pause'};
nIter = 20;
shifts = 0:8:2^8; % sampling shifts
tShifts =  [-250:5:250];% temporal shifts
tShifts =  [-100:2:50];% temporal shifts
tShifts =  [-100:50];% temporal shifts
tag = DataHash({vlb,stgrps,stlbls,tShifts});
%filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_2d_TS_',tag,'.mat']);




cha = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
% CORRECT variables, by subject 
chzdza = -circ_dist(dhzdza,pi/2);
cblen = dblen+20*double(ismember(dtind,[3,4,5]));
chroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;

% FLIP variable sign where variable is oposite to head-body angle
ac = chbang;
chbang  ( ac<0 ) = -chbang   ( ac<0 );
cfhrvl  ( ac<0 ) = -cfhrvl   ( ac<0 );
chroll  ( ac<0 ) = -chroll   ( ac<0 );
chvang  ( ac<0 ) = -chvang   ( ac<0 );
chzdza  ( ac<0 ) = -chzdza   ( ac<0 );


if exist(filepath,'file'),
    load(filepath);
else

    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),...
                    numel(phzBinCenters),...
                    numel(vbc{1}),...
                    numel(vbc{2}),...
                    2, ...
                    numel(tShifts),...
                    nIter]);
    smJpdf = zeros([3,numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2,numel(tShifts),nIter]);
    sCnt =  zeros([numel(vbc{1}),numel(vbc{2}), 2,numel(tShifts),nIter]);
    % SET ind to walk 
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{1})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<320;%& chbang<0.4;
    for v = 1:numel(vbe),
        vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
    end
    
    
    tdphz = dphz(ind);
    tsmMask = smMask(ind);

    g = 1;
    h = 2;
    xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
    yi = vbi{h};        yb = vbc{h};        yc = numel(yb);

% ITERATE over time-shifted error projected onto current head position    
%for s = 1:numel(tShifts),
for s = [35,51,67]
% COMPUTE shifted position prection projection (PPP)
        dct.error = dct.com;
        dct.errorShift = cf(@(e,x,h,t)  sq(multiprod(permute(bsxfun(@minus,                                   ...
                                           e(:,[1,2],:),                                                      ...
                                           sq(circshift(x(:,'hcom',[1,2]),-tShifts(s)))),                     ...
                                           [1,2,4,3]),                                                        ...
                                           circshift(h(:,:,:),-tShifts(s)),2,[2,3])),                         ...
                            dct.error,dct.xyz,dct.hvec,dct.tvec);
        derrlns = cf(@(e) reshape( sq(e(:,1,:))', [], 1), dct.errorShift);
        derrlns = cat( 1, derrlns{:} );
        derrlts = cf(@(e) reshape( sq(e(:,2,:))', [], 1), dct.errorShift);                      
        derrlts = cat( 1, derrlts{:} );
        lfErr = {derrlns,derrlts+8};
        if mirror,
            lfErr{2}( ac<0 ) = -lfErr{2} ( ac<0 );
        end
        

% COMPUTE JPDF(FERROR, PHZ | HBA, HRVL)        
        for e = 1:2,
            tferr = lfErr{e}(ind);
            for x = 1:xc,
                for y = 1:yc,
                    disp(['s: ',num2str(s),'   x: ',num2str(x),'  y: ',num2str(y),'  e: ',num2str(e)]);
                    for i = 1:nIter,
                        indb = xi==x ...
                               & yi==y ...
                               & circshift(tsmMask, randsample(shifts,1));
                        indb(randn([sum(indb),1]) > 0) = false;
                        sCnt(x,y,e,s,i) = sum(indb);
                        sJpdf(:,:,x,y,e,s,i) = hist2([tferr(indb),tdphz(indb)],...
                                                     ferrorBinEdges{e},...
                                                     phzBins);
                    end
                end
            end
        end

    end
% COMPUTE median JPDF(FERROR, PHZ | HBA, HRVL) accross PHZ
%for s = 1:numel(tShifts)
        for s = [35,51,67]
        for x = 1:xc,
            tic
            for y = 1:yc
                for p = 1:numel(phzBinCenters),
                    for e = 1:2,
                        for i = 1:nIter
                        try
                            [c, index] = unique(cumsum(sJpdf(:,p,x,y,e,s,i)./sum(sJpdf(:,p,x,y,e,s,i))),'last');
                            index(c==0|c==1) = [];
                            c(c==0|c==1) = [];
                            smJpdf(:,p,x,y,e,s,i) = interp1(c,ferrorBinCenters{e}(index),[0.25,0.5,0.75]);
                        end
                        end
                    end
                end
            end
            toc
        end

    end

    save(filepath,...
         'smJpdf',                                                                           ...
         'sJpdf',                                                                            ...
         'sCnt',                                                                             ...
         'vlb',                                                                              ...
         'vdc',                                                                              ...
         'vbe',                                                                              ...
         'vbc',                                                                              ...
         'vnm',                                                                              ...     
         'vun',                                                                              ...          
         'phzBins',                                                                          ...
         'stgrps',                                                                           ...
         'stlbls',                                                                           ...
         'tShifts',                                                                          ...
         '-v7.3'                                                                             ...
         );
end





hfig = figure();

clf(hfig);
hba = 1;
tpa = 2;
tpd = 6;
e = 1;


subplot(221);
hold(gca(),'on');
cmap = cool(size(smJpdf,4));
for x = 1:size(smJpdf,4);
    plot(tShifts/250*1000,                                              ... % X: time shift
         sq(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,e,:,:),ndims(smJpdf))))',  ... % Y: Forward Egocentric Decoded Position Projection on Head FOR
         'Color',cmap(x,:));
end
legend(cf(@(v) ['hvf = ',num2str(v), ' cm/s'], num2cell(vbc{2})));
xlabel('Time Shift (ms)');
ylabel({'Forward Egocentric Decoded Position','Projection on Head FOR (mm)'});
title({'Descending Theta Phase (\phi_{\theta} = 113)','over range of forward head speed (hvf)'});

subplot(222);
hold(gca(),'on');
cmap = cool(size(smJpdf,4));
for x = 1:size(smJpdf,4);
    plot(tShifts/250*1000,                                              ... % X: time shift
         sq(sq(mean(smJpdf(2,phzOrder(tpa),hba,x,e,:,:),ndims(smJpdf))))',  ... % Y: Forward Egocentric Decoded Position Projection on Head FOR
         'Color',cmap(x,:));
end
legend(cf(@(v) ['hvf = ',num2str(v), ' cm/s'], num2cell(vbc{2})));
xlabel('Time Shift (ms)');
ylabel({'Forward Egocentric Decoded Position','Projection on Head FOR (mm)'});
title({'Ascending Theta Phase (\phi_{\theta} = 247)','over range of forward head speed (hvf)'});


subplot(223);
hold('on');
cmap = cool(size(smJpdf,4));
% $$$ ttd = 44;
% $$$ tta = 101;
tta = 96;
ttd = 66;
%ttd = 35;%?
%tta = 50;%?
for x = 1:size(smJpdf,4);
    plot(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,e,:,:),ndims(smJpdf))),...
         sq(mean(smJpdf(2,phzOrder(tpa),hba,x,e,:,:),ndims(smJpdf))),...
         '-+','Color',cmap(x,:))
    plot(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,e,tta,:),ndims(smJpdf))),...
         sq(mean(smJpdf(2,phzOrder(tpa),hba,x,e,tta,:),ndims(smJpdf))),...
         'o','Color',cmap(x,:),'MarkerSize',10,'MarkerFaceColor',cmap(x,:),'MarkerEdgeColor','k')
    plot(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,e,1,:),ndims(smJpdf))),...
         sq(mean(smJpdf(2,phzOrder(tpa),hba,x,e,1,:),ndims(smJpdf))),...
         's','Color',cmap(x,:),'MarkerSize',10,'MarkerFaceColor',cmap(x,:),'MarkerEdgeColor','k')         
    plot(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,e,ttd,:),ndims(smJpdf))),...
         sq(mean(smJpdf(2,phzOrder(tpa),hba,x,e,ttd,:),ndims(smJpdf))),...
         '<','Color',cmap(x,:),'MarkerSize',10,'MarkerFaceColor',cmap(x,:),'MarkerEdgeColor','k')
    
    %plot(sq(mean(smJpdf(2,phzOrder(tpd),hba,x,2,:,:),ndims(smJpdf)))','Color',cmap(x,:))
    %plot(sq(sqrt(sum(mean(smJpdf(2,phzOrder(tpd),hba,x,:,:,:),ndims(smJpdf)).^2,5)))','Color',cmap(x,:))
end
Lines(mean(sq(mean(smJpdf(2,phzOrder(tpd),hba,:,e,ttd,:),ndims(smJpdf)))),[],'k');
Lines([],mean(sq(mean(smJpdf(2,phzOrder(tpa),hba,:,e,tta,:),ndims(smJpdf)))),'k');
ylabel({['Ascending theta phase (\phi_{\theta} = 247)' ],'Forward Decoded Position projected onto','head frame of reference (mm)'})
xlabel({'Descending theta phase (\phi_{\theta} = 113)','Forward Decoded Position projected onto','head frame of reference (mm)'})
title({'Time Shifted projected decoded postition','Squares: -400ms','Triangles: -140ms','Circles: -20ms'});

e = 1;
etheu = [];
for tta = 1:numel(tShifts),
    for tpa = 2:7,
    etheu(tta,tpa-1) = sum((sq(mean(smJpdf(2,phzOrder(tpa),hba,:,e,tta,:),ndims(smJpdf)))-mean(sq(mean( ...
                         smJpdf(2,phzOrder(tpa),hba,:,e,tta,:),ndims(smJpdf))))).^2);
    end
end


subplot(224);
hold('on');
cmap = winter(6);
for p = 1:6
    plot(tShifts/250*1000,imgaussfilt(etheu(:,p),3),'Color',cmap(p,:));
end
ylabel('SSE of mean centered forward projection for each phase');
xlabel('Time Shift (ms)');
legend(cf(@(p)  [ '\phi_{\theta} = ',num2str(p)],num2cell(round(circ_rad2ang(phzBinCenters(phzOrder(2:7))+2*pi*double(phzBinCenters(phzOrder(2:7))<0))))));
title({'Sum of Squared Errors of Forward Egocentric Decoded Position Projection','for each phase (\phi_{\theta}) of theta'});

figureName = ['TimeShifted_',el(e),'EPP_',vlb{1},'_mirrored_',vlb{2},'_',tag];
print(hfig,'-dpng', fullfile(figurePath,[figureName,'.png']));



k = 0;
%for s = [1:5];
%for s = [51];
for s = [35,51,67];
    %for s = [60,96,130];    
    %for s = [1,6,7,8];
    k = k+1;
    hfig = figure();
    hfig.Units = 'centimeters';
    hfig.Position = [(k-1).*6,5,6,24];
    for p = 1:8,
        for e = 1:2
            subplot2(8,2,p,e);
            nmask = double(sq(mean(sCnt(:,:,e,s,1:nIter),ndims(sCnt)))>300);
            nmask(~nmask) = nan;
            %nmask = 1;
            if e==1,
                ca = [-120,200];
            else
                ca = [-60, 60];
            end
 % $$$         imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(smJpdf(2,phzOrder(p),:,:,e,g,s,:))'.*nmask'},...
 % $$$                    ca,'linear',false,'colorMap',@jet);**
 % $$$         axis('xy');
            imagesc(vbc{1},vbc{2},sq(median(smJpdf(2,phzOrder(p),:,:,e,s,:),ndims(smJpdf)))'.*nmask');
            colormap('jet');
            caxis(ca)
            axis('xy');
            
            if e==1&&p==1,
                %title(stlbls{s});
            end
        end
    end
end



j = 1;
for s = [60,96];
%    for s = [35,51];    
%    for s = [51];            
j = j+1;
k = 1;
for x = [1:size(smJpdf,3)];
    k = k+1;
    hfig = figure();
    hfig.Units = 'centimeters';
    hfig.Position = [(k-1).*10,j*24-48,10,24];

for y = 1:numel(vbc{2});
    subplot2(numel(vbc{2}),2,y,1);
    imagesc(ferrorBinCenters{1},...
            phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0), ...
            imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,y,1,s,:),ndims(sJpdf))),[2,0.1])');
    axis('xy');
    Lines(0,[],'k');
    Lines([],pi,'k');
    xlim([-250,350]);
    subplot2(numel(vbc{2}),2,y,2);
    imagesc(ferrorBinCenters{1},...
            phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0), ...
            imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,y,2,s,:),ndims(sJpdf))),[2,0.1])');
    axis('xy');
    Lines(0,[],'k');
    Lines([],pi,'k');
    xlim([-250,250]);
end
colormap('jet');
end
end

for s = 1:numel(tShifts)
    for x = 1:numel(vbc{1})
        for y = 1:numel(vbc{2})
            for e = 1:2,
                for p = 1:numel(phzOrder),
                    dpStd(1,p,x,y,e,s) = std(ferrorBinCenters{1}',                                   ... bins
                                             sq(mean(sJpdf(:,phzOrder(p),x,y,1,s,:),ndims(sJpdf)))); %   weights
                end
            end
        end
    end
end


cmap = cool(numel(vbc{2}));
hfig  = figure();

clf();
for p = 1:numel(phzOrder),
    subplot(numel(phzOrder),1,p);
    hold('on');    
    for x = 1:numel(vbc{2}),
        plot(sq(dpStd(1,p,1,x,1,:))','Color',cmap(x,:))
        xlim([1,101])
        ylim([100,200]);
    end
end

e = 1;
x = 2;
s = 3;
switch s,
  case 1,
    span = 1:numel(vbc{2});
  case 2,
    span = 2:numel(vbc{2})-3;
  case 3,
    span = 3:numel(vbc{2});
  case 4,
    span = 1:numel(vbc{2})-1;
end

figure();
hold('on');
cmap = cool(numel(vbc{2}));
for y = span,
    for i = 1:nIter,
    plot(sq(smJpdf(2,phzOrder(1:end),x,y,e,s,i)),'-+','Color',cmap(y,:));
    end
end
for y = span,
    plot(sq(mean(smJpdf(2,phzOrder(1:end),x,y,e,s,:),ndims(smJpdf))),'-+','Color','k','LineWidth',2);
end
ylim([-25,135]);








g = 1;
h = 2;
figure();
k = 0;
for s = [1:5]
    %for s = [1,6,7,8];
    k = k +1;
    subplot(1,4,k);
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{s})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;% & chbang <0.25;
    hist2([eval([vnm{g},'(ind)']),eval([vnm{h},'(ind)'])],...
          linspace([vbe{g}([1,end]),50]),...
          linspace([vbe{h}([1,end]),50]));
    Lines(vbe{g},[],'k');
    Lines([],vbe{h},'k');
    title(stlbls{s});
end




%%%>>>
