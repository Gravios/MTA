
MTA_PROJECT_REPORT_PATH = '/storage/share/Projects/BehaviorPlaceCode/';
figurePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP','bySession');
create_directory(figurePath);
phzOrder = [5:8,1:4];

figure();
    subplot(121);hist(dhbang((logical(dstcm(:,5))|logical(dstcm(:,3)))&ismember(dtind,[1,2])),100)
    subplot(122);hist(dhbang((logical(dstcm(:,5))|logical(dstcm(:,3)))&~ismember(dtind,[1,2])),100)

figure()
    subplot(121);hist(dblen((logical(dstcm(:,5))|logical(dstcm(:,3)))&ismember(dtind,[1,2])),100)
    subplot(122);hist(dblen((logical(dstcm(:,5))|logical(dstcm(:,3)))&~ismember(dtind,[1,2])),100)

chbang = -dhbang;
cblen = dblen+20*double(ismember(dtind,[3,4,5]));


vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below


%%%<<< Variable Setup
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
%vbe{end+1} = [-1.2,-0.8,-0.4,0.4,0.8,1.2];
%vbe{end+1} = [-1.2,-0.8,-0.4,0.4,0.8,1.2];
vbe{end+1} = [-1.2,-0.3,0.3,1.2];
%vbe{end+1} = [-1.2,-0.6,-0.2,0.2,0.6,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
% $$$ vlb{end+1} = 'hvl';
% $$$ vdc{end+1} = 'lateral head speed';
% $$$ vnm{end+1} = 'dfhrvl';
% $$$ %vbe{end+1} = [-10,10,30,50,80];
% $$$ %vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
% $$$ %vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
% $$$ %vbe{end+1} = [-4,4,10,25,40,55,70,85];
% $$$ vbe{end+1} = [-80,-5,5,80];
% $$$ vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
% $$$ vun{end+1} = 'cm/s';

vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
%vbe{end+1} = [-10,10,30,50,80];
%vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
%vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
%vbe{end+1} = [-4,4,10,25,40,55,70,85];
vbe{end+1} = [-4,4,20,80];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>

%%%<<< computations
stgrps = {9:11};
stlbls = {'ALL'};
nIter = 100;
shifts = 0:8:2^8;
tag = DataHash({vlb,vbe,stgrps,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_2d_PS_ca3_',tag,'.mat']);

lfErr = ferr([1:2]);

if exist(filepath,'file'),
    load(filepath);
else

    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2, numel(tind),nIter]);
    smJpdf = zeros([3,numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2,numel(tind),nIter]);
    sCnt =  zeros([numel(vbc{1}),numel(vbc{2}), 2,numel(tind),nIter]);
    for s = 1:numel(tind),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,[9:11])),2)               ...
              & dtind==tind(s)                                ...
              & duincI                                        ...
              & dpostI                                        ...        
              & dhdist<320;% & chbang<0.4;
        for v = 1:numel(vbe),
            vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
        end
        tdphz = dphz(ind);
        tsmMask = smMask(ind);
        for g = 1,
            for h = 2,
                xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
                yi = vbi{h};        yb = vbc{h};        yc = numel(yb);
                for e = 1:2,
                    tferr = lfErr{e}(ind);
                    % compute JPDF(FERROR, PHZ | HBA, HRVL)
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
        end
    end

    for s = 1:numel(tind)
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
         '-v7.3'                                                                             ...
         );
end
%%%>>>




%%%<<< Summary Figures
xc = numel(vbc{1});
yc = numel(vbc{2});
el = 'FL';
elbl = {'Forward','Lateral'};
hfig = figure();
for s = 1:4;
    for e = 1:2;

        clf(hfig)
        hfig.Units = 'Centimeters';
        hfig.Position(3:4) = [20,10];
        sax = reshape(tight_subplot(yc,xc,[0.01,0.01],0.25,0.1),[xc,yc])';
        for x = 1:xc;
            for y = 1:yc;
                axes(sax(y,x));
                imagesc(ferrorBinCenters{1},...
                        phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0), ...
                        imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,y,e,s,:),ndims(sJpdf))),[2,0.1])');
                axis('xy');
                Lines(0,[],'k');
                Lines([],pi,'k');
                xlim([-250,250+100*double(e==1)]);
                if round(xc/2)==x & yc==y;
                    xlabel(vlb{1});
                end        
                if round(xc/2)==x & y==1;
                    title({Trials{s}.filebase,[vdc{1},' vs ',vdc{2}],[elbl{e} ' Egocentric Phase Preccession']});
                end        
                if x==1 & round(yc/2)==y;
                    ylabel(vlb{2});
                end        
                if x==1&y==yc,
                    xlabel('mm');
                    ylabel({'Theta','phase'});
                end
                if x~=1 | y~=yc,
                    sax(y,x).XTick = [];
                end
                if y~=yc | x>1,
                    sax(y,x).YTick = [];
                end
            end
        end
        colormap('jet');
        drawnow();
        fax = axes('Position',[0,0,1,1],'Visible','off');
        xlim(fax,[0,1]);
        ylim(fax,[0,1]); 

        % LABEL x bins
        line([sax(yc,1).Position(1),sax(yc,xc).Position(1)+sax(yc,xc).Position(3)],...
             (sax(yc,xc).Position(2)-0.1).*[1,1],'Color','k');
        for x = 1:xc,
            text(sax(yc,x).Position(1)+sax(yc,x).Position(3)/2,...
                 sax(yc,xc).Position(2)-0.13,...
                 ['( ',num2str(vbe{1}(x)),' to ',num2str(vbe{1}(x+1)),' )'],...
                 'HorizontalAlignment','center');
        end
        text(mean([sax(yc,1).Position(1),sax(yc,xc).Position(1)+sax(yc,xc).Position(3)]),...
             sax(yc,xc).Position(2)-0.19,...
             {vdc{1},vun{1}},...
             'HorizontalAlignment','center');

        % LABEL y bins
        line((sum(sax(yc,xc).Position([1,3]))+0.02).*[1,1],...
             [sax(yc,xc).Position(2),sum(sax(1,xc).Position([2,4]))],...
             'Color','k');
        for y = 1:yc,
            text(sum(sax(yc,xc).Position([1,3]))+0.03,...
                 sax(y,xc).Position(2)+sax(y,xc).Position(4)/2,...
                 ['( ',num2str(vbe{2}(y)),' to ',num2str(vbe{2}(y+1)),' )'],...
                 'HorizontalAlignment','center',...
                 'Rotation',90);
        end
        text(sum(sax(yc,xc).Position([1,3]))+0.06,...
             sum(sax(round(yc/2),xc).Position([2,4]).*[1,0.5]),...
             {vdc{2},vun{2}},...
             'HorizontalAlignment','center',...
             'Rotation',90);

        figureName = [el(e),'EPP_',vlb{1},'_',vlb{2},'_',Trials{s}.filebase];
        print(hfig,'-dpng', fullfile(figurePath,[figureName,'.png']));
    end
end

%print(hfig,'-depsc',[filePath,'.eps']);
%%%>>>

figure();
x = 3;y = 2;
plot(sq(mean(smJpdf(2,phzOrder(2),x,y,1,:,:),ndims(smJpdf))),...
     sq(mean(smJpdf(2,phzOrder(2),x,y,2,:,:),ndims(smJpdf))),...
     '.');
xlim([-200,200])
ylim([-200,200])

figure();
hold('on');
x = 1;y = 2;
plot(sq(mean(smJpdf(2,phzOrder(5),x,y,1,:,:),ndims(smJpdf))),...
     sq(mean(smJpdf(2,phzOrder(5),x,y,2,:,:),ndims(smJpdf))),...
     '.k','MarkerSize',20);
plot(sq(mean(smJpdf(2,phzOrder(2),x,y,1,:,:),ndims(smJpdf))),...
     sq(mean(smJpdf(2,phzOrder(2),x,y,2,:,:),ndims(smJpdf))),...
     '.r','MarkerSize',20);

xlim([-50,200])
ylim([-100,100])



hfig = figure();
hfig.PaperOrientation = 'landscape';
cmap = jet(size(smJpdf,3));
cmap = ['bgr']';
for y = 1:size(smJpdf,4),
    subplot(1,3,y);
    hold('on');
    for x = 1:size(smJpdf,3),
        for s = 1:7,
            %if mean(sCnt(x,y,e,s,:),ndims(sCnt))<250,
            %    continue;
            %end
            
            plot(sq(mean(smJpdf(2,phzOrder([2,6]),x,y,2,s,:),ndims(smJpdf))),...
                 sq(mean(smJpdf(2,phzOrder([2,6]),x,y,1,s,:),ndims(smJpdf)))+8,...
                 'Color',cmap(x,:),'LineWidth',2);
            plot(sq(mean(smJpdf(2,phzOrder([6]),x,y,2,s,:),ndims(smJpdf))),...
                 sq(mean(smJpdf(2,phzOrder([6]),x,y,1,s,:),ndims(smJpdf)))+8,...
                 '.','Color','k','MarkerSize',10);
            
        end
    end
    ylim([-100,150]);
    xlim([-100,100]);
    daspect([1,1,1]);
    grid('on');
end





figure();
cmap = cool(10);
%cmap = ['rbg']';
for y = 1:3,
    subplot(3,1,y);
    hold('on');
    for x = 1:5;    
        for s = 1:10; 
            plot(sq(mean(smJpdf(2,phzOrder([2]),:,y,1,s,:),ndims(smJpdf))),...
                 sq(mean(smJpdf(2,phzOrder([2]),:,y,2,s,:),ndims(smJpdf)))+8,...
                 'Color',cmap(s,:),'LineWidth',2);
        end
    end
    xlim([-100,200])
ylim([-150,150])

end

plot(sq(mean(smJpdf(2,phzOrder(2),x,y,1,s,:),ndims(smJpdf))),...
     sq(mean(smJpdf(2,phzOrder(2),x,y,2,s,:),ndims(smJpdf))),...
     '.r','MarkerSize',20);

xlim([-50,200])
ylim([-100,100])



% DEMONSTRATED the EPP dependence upon HBA and HVL

% MIRROR LEPP and get stats




chzdza = -circ_dist(dhzdza,pi/2);
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;
lfErr = ferr([1:2]);
a = chbang;
chbang  ( a<0 ) = -chbang   ( a<0 );
lfErr{2}( a<0 ) = -lfErr{2} ( a<0 );
cfhrvl  ( a<0 ) = -cfhrvl   ( a<0 );
chvang  ( a<0 ) = -chvang   ( a<0 );
clear('a');



vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below


%%%<<< HBA x HVF
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.3,0.6,0.9,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
vlb{end+1} = 'hvl';
vdc{end+1} = 'lateral head speed';
vnm{end+1} = 'dfhrvl';
vbe{end+1} = [-60,-4,4,60];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>

stgrps = {9:11}
stlbls = {'ALL'};
nIter = 20;
shifts = 0:8:2^8;
tag = DataHash({vlb,stgrps,nbins,nvars,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_2d_PS_Mirrored_',tag,'.mat']);

if exist(filepath,'file'),
    load(filepath);
else
    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2, numel(tind),nIter]);
    smJpdf = zeros([3,numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2,numel(tind),nIter]);
    sCnt =  zeros([numel(vbc{1}),numel(vbc{2}), 2,numel(tind),nIter]);
    for s = 1:numel(tind),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,[9:11])),2)               ...
              & dtind==tind(s)                                ...
              & duincI                                        ...
              & dpostI                                        ...        
              & dhdist<320;% & chbang<0.4;
        for v = 1:numel(vbe),
            vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
        end
        tdphz = dphz(ind);
        tsmMask = smMask(ind);
        for g = 1,
            for h = 2,
                xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
                yi = vbi{h};        yb = vbc{h};        yc = numel(yb);
                for e = 1:2,
                    tferr = lfErr{e}(ind);
                    % compute JPDF(FERROR, PHZ | HBA, HRVL)
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
        end
    end

    for s = 1:numel(tind)
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
         '-v7.3'                                                                             ...
         );
end


xc = size(sJpdf,3);
yc = size(sJpdf,4);
el = 'FL';

for s = 1:numel(tind);
    for e = 1:2;

    hfig = figure();
        clf(hfig)
        hfig.Units = 'Centimeters';
        hfig.Position(3:4) = [20,10];
        sax = reshape(tight_subplot(yc,xc,[0.01,0.01],0.25,0.1),[xc,yc])';        
        for x = 1:xc
            for y = 1:yc
                axes(sax(y,x));
                imagesc(ferrorBinCenters{1},pbc(phzOrder),imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,yc+1-y,e,s,:),ndims(sJpdf)))',[0.1,3]));
                axis('xy');

                Lines(0,[],'k');
                Lines([],pi,'k');
                xlim([-250,250+100*double(e==1)]);
                if round(xc/2)==x & y==1;
                    title({Trials{s}.filebase,[vdc{1},' vs ',vdc{2}],'Mirrored Lateral Egocentric Phase Preccession'});
                end        
                if x==1&y==yc,
                    xlabel('mm');
                    ylabel({'Theta','phase'});
                end
                if x~=1 | y~=yc,
                    sax(y,x).XTick = [];
                end
                if y~=yc | x>1,
                    sax(y,x).YTick = [];
                end
            end
        end
        colormap('jet');
        drawnow();
        fax = axes('Position',[0,0,1,1],'Visible','off');
        xlim(fax,[0,1]);
        ylim(fax,[0,1]); 

        % LABEL x bins
        line([sax(yc,1).Position(1),sax(yc,xc).Position(1)+sax(yc,xc).Position(3)],...
             (sax(yc,xc).Position(2)-0.1).*[1,1],'Color','k');
        for x = 1:xc,
            text(sax(yc,x).Position(1)+sax(yc,x).Position(3)/2,...
                 sax(yc,xc).Position(2)-0.13,...
                 ['( ',num2str(vbe{1}(x)),' to ',num2str(vbe{1}(x+1)),' )'],...
                 'HorizontalAlignment','center');
        end
        text(mean([sax(yc,1).Position(1),sax(yc,xc).Position(1)+sax(yc,xc).Position(3)]),...
             sax(yc,xc).Position(2)-0.19,...
             {vdc{1},vun{1}},...
             'HorizontalAlignment','center');

        % LABEL y bins
        line((sum(sax(yc,xc).Position([1,3]))+0.02).*[1,1],...
             [sax(yc,xc).Position(2),sum(sax(1,xc).Position([2,4]))],...
             'Color','k');
        for y = 1:yc,
            text(sum(sax(yc,xc).Position([1,3]))+0.03,...
                 sax(y,xc).Position(2)+sax(y,xc).Position(4)/2,...
                 ['( ',num2str(vbe{2}(yc+1-y)),' to ',num2str(vbe{2}(yc+1-y+1)),' )'],...
                 'HorizontalAlignment','center',...
                 'Rotation',90);
        end
        text(sum(sax(yc,xc).Position([1,3]))+0.06,...
             sum(sax(round(yc/2),xc).Position([2,4]).*[1,0.5]),...
             {vdc{2},vun{2}},...
             'HorizontalAlignment','center',...
             'Rotation',90);

% $$$         figureName = [el(e),'EPP_',vlb{1},'_mirrored_',vlb{2},'_',Trials{s}.filebase];
% $$$         print(hfig,'-dpng', fullfile(figurePath,[figureName,'.png']));
    end
end



figure,
for s = 1:10,
    for p = 1:8,
    subplot2(8,10,p,s);
    imagesc(sq(mean(smJpdf(2,phzOrder(p),:,:,2,s,:),ndims(smJpdf)))');axis('xy');
    caxis([-100,100]);
    colormap('jet');
    end
end

