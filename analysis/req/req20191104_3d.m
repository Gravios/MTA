
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)            ...
        & duincI                                        ...
        & dpostI                                        ...        
        & dhdist<380;

figure()
hist(dhz(ind&~ismember(dtind,[3,4,5])),100);


cblen = dblen+20*double(ismember(dtind,[3,4,5]));
chroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;
lfErr = ferr([1:2]);
v = cfhrvl;
a = chbang;
chbang  ( a<0 ) = -chbang   ( a<0 );
lfErr{2}( a<0 ) = -lfErr{2} ( a<0 );
cfhrvl  ( a<0 ) = -cfhrvl   ( a<0 );
chroll ( a<0 ) = -chroll  ( a<0 );
chvang  ( a<0 ) = -chvang   ( a<0 );
% $$$ cfhrvl( a<0 & v>0) = -cfhrvl( a<0 & v>0);
% $$$ cfhrvl( a<0 & v<0) = -cfhrvl( a<0 & v<0);
clear('a');


nbins = 3;
vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below

% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
%vbe{end+1} = [0,0.12,0.26,0.44,0.7,1.2];
switch nbins
  case 3, vbe{end+1} = [0,0.4,0.8,1.2];
  case 5,  vbe{end+1} = linspace(0,1.2,6);
  case 7,  vbe{end+1} = linspace(0,1.2,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% 2
vlb{end+1} = 'hva';
vdc{end+1} = 'head anglular velocity';
vnm{end+1} = 'chvang';
switch nbins
  case 3,  vbe{end+1} = [-0.2,-0.012,0.012,0.2];
  case 5,  vbe{end+1} = [-0.2,-0.06,-0.012,0.012,0.06,0.2];    
  case 7,  vbe{end+1} = [-0.2,-0.1,-0.06,-0.012,0.012,0.06,0.1,0.2];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad/sec';

% 3
vlb{end+1} = 'hrl';
vdc{end+1} = 'head roll';
vnm{end+1} = 'chroll';
switch nbins
  case 3,  vbe{end+1} = [-0.5,-0.1,0.1,0.5];
  case 5,  vbe{end+1} = linspace(-0.5,0.5,6);
  case 7,  vbe{end+1} = linspace(-0.5,0.5,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% 4
vlb{end+1} = 'hvl';
vdc{end+1} = 'lateral head speed';
vnm{end+1} = 'cfhrvl';
switch nbins
  case 3,  vbe{end+1} = [-80,-5,5,80];    
  case 5,  vbe{end+1} = [-40,-16,-4,4,16,40];
  case 7,  vbe{end+1} = [-36,-24,-12,-2,2,12,24,36];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';

% 5
vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
switch nbins
  case 3,  vbe{end+1} = [-4,4,12,64];
  case 5,  vbe{end+1} = [-4,4,12,24,36,64];    
  case 7,  vbe{end+1} = [-4,0,4,10,20,30,40,80];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';

% 6
vlb{end+1} = 'hbd';
vdc{end+1} = 'head-body distance';
vnm{end+1} = 'cblen';
switch nbins
  case 3,  vbe{end+1} = [100,120,140,160];
  case 5,  vbe{end+1} = linspace(100,160,6);
  case 7,  vbe{end+1} = linspace(100,170,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'mm';

% 7
vlb{end+1} = 'hbp';
vdc{end+1} = 'head-body pitch';
vnm{end+1} = 'dfet';
switch nbins
  case 3,  vbe{end+1} = [-1.4,-0.8,-0.2,0.2];
  case 5,  vbe{end+1} = linspace(-1.4,0.2,6);
  case 7,  vbe{end+1} = linspace(-1.4,0.2,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';


% $$$ stgrps = {[9:11],4,6,3,5,9,10,11};
% $$$ stlbls = {'all','HP','LP','HW','LW','T','W','P'};
stgrps = {[9:11],9,10,11};
stlbls = {'all','T','W','P'};
nvars = numel(vlb);

tag = DataHash({vlb,stgrps,nbins,nvars,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_3d_',tag,'.mat']);

if exist(filepath,'file'),
    load(filepath);
else
    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),nbins, nbins, nbins, 2, nvars, nvars, nvars, numel(stgrps)]);
    smJpdf = zeros([3,numel(phzBinCenters),nbins, nbins, nbins, 2, nvars, nvars, nvars, numel(stgrps)]);
    sCnt =  zeros([nbins, nbins, nbins, 2,nvars, nvars, nvars, numel(stgrps)]);
    for s = 1:numel(stgrps),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,stgrps{s})),2)            ...
              & duincI                                        ...
              & dpostI                                        ...        
              & dhdist<380;% & chbang <0.25;
        for v = 1:numel(vbe),
            vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
        end
        tdphz = dphz(ind);
        for g = 1:nvars,
            for h = g+1:nvars,
                for d = h+1:nvars,
                    xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
                    yi = vbi{h};        yb = vbc{h};        yc = numel(yb);
                    zi = vbi{d};        zb = vbc{d};        zc = numel(zb);
                    for e = 1:2;
                        tferr = lfErr{e}(ind);
                        % compute JPDF(FERROR, PHZ | HBA, HRVL)
                        for x = 1:xc,
                            for y = 1:yc,
                                for z = 1:zc
disp(['s: ',num2str(s),'   g: ',num2str(g),'   h: ',num2str(h),'   d: ',num2str(d),'   x: ',num2str(x),'  y: ',num2str(y),'  z: ',num2str(z),'  e: ',num2str(e)]);
                                    indb = xi==x & yi==y & zi==z;
                                    sCnt(x,y,z,e,g,h,d,s) = sum(indb);
                                    sJpdf(:,:,x,y,z,e,g,h,d,s) = hist2([tferr(indb),tdphz(indb)],...
                                                                       ferrorBinEdges{e},...
                                                                       phzBins);
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    for s = 1:numel(stgrps),
        for g = 1:nvars,            
            for h = g+1:nvars,                        
                for d = h+1:nvars,
disp(['s: ',num2str(s),'   g: ',num2str(g),'   h: ',num2str(h),'   d: ',num2str(d)]);
                    for x = 1:xc,
                        tic
                        for y = 1:yc
                            for z = 1:zc,
                                for p = 1:numel(phzBinCenters),
                                    for e = 1:2,
                                        try
                                            [c, index] = unique(cumsum(sJpdf(:,p,x,y,z,e,g,h,d,s)./sum(sJpdf(:,p,x,y,z,e,g,h,d,s))),'last');
                                            index(c==0|c==1) = [];
                                            c(c==0|c==1) = [];
                                            smJpdf(:,p,x,y,z,e,g,h,d,s) = interp1(c,ferrorBinCenters{e}(index), ...
                                                                              [0.25,0.5,0.75]);
                                        end
                                    end
                                end
                            end
                        end
                        toc
                    end
                end
            end
        end
    end

    save(fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_3d_',tag,'.mat']),...
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



g = 1;
h = 5;
d = 6;
k = 0;
for s = [2,3];
    %for s = [1,6,7,8];
    k = k+1;        
    hfig = figure();
    hfig.Units = 'centimeters';
    hfig.Position = [(k-1).*18,5,18,24];
    
    for z = 1:3,
        for p = 1:8,
            for e = 1:2
                subplot2(8,nbins*2,p,z+(e-1).*3);
                nmask = sq(double(sCnt(:,:,z,2,g,h,d,s)>10000));
                nmask(~nmask) = nan;
                if e==1,
                    ca = [-60,120];
                else
                    ca = [-60, 60];
                end
% $$$         imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(smJpdf(2,phzOrder(p),:,:,e,g,s,:))'.*nmask'},...
% $$$                    ca,'linear',false,'colorMap',@jet);
% $$$         axis('xy');
                imagesc(vbc{g},vbc{h},sq(smJpdf(2,phzOrder(p),:,:,z,e,g,h,d,s,:))'.*nmask');
                colormap('jet');
                caxis(ca)
                axis('xy');
                if e==1&&p==1&&z==1,
                    title(stlbls{s});
                end
            end
        end
    end
end
%%%>>>
