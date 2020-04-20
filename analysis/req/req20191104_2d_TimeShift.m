





dct.error = dct.com;
dct.errorShift = cf(@(e,x,h,t)  cat(2,                                                                ...
                           sq(multiprod(permute(bsxfun(@minus,                                   ...
                                                       e(:,[1,2],:),                             ...
                                                       sq(circshift(x(:,'hcom',[1,2]),25))),                   ...
                                                [1,2,4,3]),...
                                        circshift(h(:,:,:),25),2,[2,3])),...
                           sq(multiprod(permute(bsxfun(@minus,                                   ...
                                                       e(:,[1,2],:),                             ...
                                                       sq(x(:,'hcom',[1,2]))),                   ...
                                                [1,2,4,3]),...
                                        t(:,:,:),2,[2,3]))),...                               
               dct.error,dct.xyz,dct.hvec,dct.tvec);
derrlns = cf(@(e) reshape( sq(e(:,1,:))', [], 1), dct.errorShift);                      derrlns = cat( 1, derrlns{:} );
derrlts = cf(@(e) reshape( sq(e(:,2,:))', [], 1), dct.errorShift);                      derrlts = cat( 1, derrlts{:} );


chzdza = -circ_dist(dhzdza,pi/2);
cblen = dblen+20*double(ismember(dtind,[3,4,5]));
chroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;
lfErr = ferr([1:2]);
lfErr = {derrlns,derrlts};
v = cfhrvl;
a = chbang;
chbang  ( a<0 ) = -chbang   ( a<0 );
lfErr{2}( a<0 ) = -lfErr{2} ( a<0 );
cfhrvl  ( a<0 ) = -cfhrvl   ( a<0 );
chroll ( a<0 ) = -chroll  ( a<0 );
chvang  ( a<0 ) = -chvang   ( a<0 );
chzdza  ( a<0 ) = -chzdza   ( a<0 );
% $$$ cfhrvl( a<0 & v>0) = -cfhrvl( a<0 & v>0);
% $$$ cfhrvl( a<0 & v<0) = -cfhrvl( a<0 & v<0);
clear('a');



%%%<<< Indepth 2d 


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
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
%vbe{end+1} = [-10,10,30,50,80];
%vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
%vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
vbe{end+1} = [-4,4,10,25,40,55,70,85];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>

%%%<<< HBA x HVF
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
vlb{end+1} = 'hvl';
vdc{end+1} = 'lateral head speed';
vnm{end+1} = 'dfhrvl';
%vbe{end+1} = [-10,10,30,50,80];
%vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
%vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
vbe{end+1} = [-60,-40,-25,-10,-4,4,10,25,40,60];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>


%%%<<< HBA x HBD
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
% 2
vlb{end+1} = 'hbd';
vdc{end+1} = 'head-body distance';
vnm{end+1} = 'cblen';
vbe{end+1} = linspace(110,160,9);
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'mm';
%%%>>>



stgrps = {9:11,2,9,10,11};
stlbls = {'ALL','Rear','Turn','Walk','Pause'};
nIter = 20;
shifts = 0:8:2^8;
tag = DataHash({vlb,stgrps,nbins,nvars,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_2d_ID_',tag,'.mat']);

if exist(filepath,'file'),
    load(filepath);
else

    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2, numel(stgrps),nIter]);
    smJpdf = zeros([3,numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2,numel(stgrps),nIter]);
    sCnt =  zeros([numel(vbc{1}),numel(vbc{2}), 2,numel(stgrps),nIter]);
    for s = 4%:numel(stgrps),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,stgrps{s})),2)            ...
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

    for s = 4%:numel(stgrps)
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




k = 0;
%for s = [1:5];
for s = [4];
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
            title(stlbls{s});
        end
    end
end
end

s = 4;
x = 1;
figure();
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
