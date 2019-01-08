% req20180926
%  Tags: 
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Project: MjgER2016:placefields
%  Description: egocentric pfs spk position 
%  Protocol: 
%    1. 
%    2  compute ratemaps


global MTA_PROJECT_PATH;

MjgER2016_load_data();


% SET analysis parameters
sampleRate = 250;   % Hz
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};%,'ripple'};
        
thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;


if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
numComp = size(eigVec{1},2);
pfindex = 1;
MjgER2016_load_bhv_erpPCA_scores();
% output:
%    fsrcz
%    FSrC
%    rmaps
clear('pfd','tags','eigVec','eigVar','eigScore','validDims','zrmMean','zrmStd',...
      'clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','rmapsShuffledMean',...
      'rmapsShuffled','FSrM','FSrS','fsrsMean','fsrsStd','fsrsMean','fsrsStd');


dataFilePath = fullfile(MTA_PROJECT_PATH,'analysis',...
                        ['req20180926-data-',DataHash({[sessionList.sessionName],sampleRate,states}),'.mat']);

if ~exist(dataFilePath,'file')||overwrite,

    % ACCUMULATE spike info
    states = {'theta','rear','walk','turn','pause','hloc','hpause','lloc','lpause'};
    sid = [];
    pfhr = [];
    pfmr = [];
    pfmrb = [];
    phzv = [];
    vxya = [];
    anga = [];
    angab = [];
    mcaa = [];
    mcda = [];
    drza = [];
    draa = [];
    stca = [];
    pfsp = [];
    pfsi = [];
    pfed = [];
    pfmx = [];
    pfui = [];
    unitSessionMap = [];
    ucount = 1;
    for t = 1:23;
        
        Trial = Trials{t}; 
        unitSubset = units{t};        
        subjectId = regexp(Trial.name,'^(\w*)-','tokens');
        subjectId = subjectId{1}{1};

        switch subjectId,
          case 'jg05'        
            pft = MTAApfs(Trial,unitSubset,[],[],'CA1thetaCA3inputPhase');
          otherwise  
            pft = pfs_2d_theta(Trial,unitSubset);
        end

        stc = Trial.load('stc','msnn_ppsvd_raux');
        xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
        fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
        mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
        mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                                  atan2(diff(xyz(:,{'hcom','nose'},2),1,2),...
                                        diff(xyz(:,{'hcom','nose'},1),1,2)));
        lfp = load(Trial,'lfp',sessionList(t).thetaRef);
        phz = lfp.phase([6,12]);
        phz.data = unwrap(phz.data);
        phz.resample(xyz);    
        phz.data = mod(phz.data+pi,2*pi)-pi;
        lfp.resample(xyz);    

        tvec = circshift(fxyz(:,'hcom',[1,2]),-50)-fxyz(:,'hcom',[1,2]);
        tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
        tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);

        tvecb = -circshift(fxyz(:,'hcom',[1,2]),50)-fxyz(:,'hcom',[1,2]);
        tvecb = sq(bsxfun(@rdivide,tvecb,sqrt(sum(tvecb.^2,3))));
        tvecb = cat(3,tvecb,sq(tvecb)*[0,-1;1,0]);
        
        hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
        hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
        hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
        [drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
        
        vxy = vel(filter(copy(xyz),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
        vxy.data(vxy.data<1e-3) = 1e-3;
        vxy.data = log10(vxy.data);

        headAngle = atan2(diff(fxyz(:,{'hcom','nose'},2),1,2),diff(fxyz(:,{'hcom','nose'},1),1,2));
        vang = [circ_dist(circshift(headAngle,-25),headAngle),circ_dist(circshift(headAngle,-75),headAngle),...
                circ_dist(circshift(headAngle,-125),headAngle),circ_dist(circshift(headAngle,-200),headAngle)];
        vangb = [-circ_dist(circshift(headAngle,25),headAngle),-circ_dist(circshift(headAngle,75),headAngle),...
                 -circ_dist(circshift(headAngle,125),headAngle),-circ_dist(circshift(headAngle,200),headAngle)];
        
        stcm = stc2mat(stc,xyz,states);

        nq = get(Trial.load('nq'),'nq');
        edist = nq.eDist(unitSubset)';
        
        spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');    
        
        [~,sind] = sort(pft.data.clu);
        si = subsref(pft.data.si(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));
        spar = subsref(pft.data.spar(ismember(pft.data.clu,unitSubset)),substruct('()',{sind}));    
        
        for unit = unitSubset(edist>20),%&pft.data.spar>0.15&pft.data.spar<0.3),
            [mxr,mxp] = pft.maxRate(unit);
            pfsCenterHR = MTADfet.encapsulate(Trial,                                               ...
                                              multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),...
                                                        hvec,2,[2,3]),                             ...
                                              sampleRate,                                          ...
                                              'placefield_center_referenced_to_head',              ...
                                              'pfsCenterHR',                                       ...
                                              'p'                                                  ...
                                              );
            pfsCenterMR = MTADfet.encapsulate(Trial,                                               ...
                                              multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),...
                                                        tvec,2,[2,3]),                             ...
                                              sampleRate,                                          ...
                                              'placefield_center_referenced_to_head',              ...
                                              'pfsCenterHR',                                       ...
                                              'p'                                                  ...
                                              );
            pfsCenterMRb = MTADfet.encapsulate(Trial,                                               ...
                                               multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),...
                                                         tvecb,2,[2,3]),                             ...
                                               sampleRate,                                          ...
                                               'placefield_center_referenced_to_head',              ...
                                               'pfsCenterHR',                                       ...
                                               'p'                                                  ...
                                               );
            res = spk(unit);
            res(res>size(xyz,1)) = [];
            pfhr = cat(1,pfhr,pfsCenterHR(res,:));
            pfmr = cat(1,pfmr,pfsCenterMR(res,:));        
            pfmrb = cat(1,pfmrb,pfsCenterMRb(res,:));                
            phzv = cat(1,phzv,phz(res,spk.map(spk.map(:,1)==unit,2)));
            vxya = cat(1,vxya,vxy(res,:));
            anga = cat(1,anga,vang(res,:));
            angab = cat(1,angab,vangb(res,:));        
            mcaa = cat(1,mcaa,mazeCenterAng(res));
            mcda = cat(1,mcda,mazeCenterDist(res));        
            drza = cat(1,drza,drz(res,unit==unitSubset));
            draa = cat(1,draa,drang(res,unit==unitSubset));
            stca = cat(1,stca,stcm(res,:));
            pfmx = cat(1,pfmx,mxp);
            pfui = cat(1,pfui,repmat(ucount,[numel(res),1]));
            ucount = ucount+1;
        end 
        unitSessionMap = cat(1,unitSessionMap,[repmat(t,[numel(unitSubset),1]),unitSubset']);
        pfsi = cat(1,pfsi,si');
        pfsp = cat(1,pfsp,spar');
        pfed = cat(1,pfed,edist');
        sid = cat(1,sid,repmat(t,[size(spar',1),1]));
    end
    pfmd = sqrt(sum(pfmx.^2,2));



    save(dataFilePath,...
         'states',          ... 
         'sid'   ,          ...
         'pfhr'  ,          ...      
         'pfmr'  ,          ...
         'pfmrb' ,          ...
         'phzv'  ,          ...
         'vxya'  ,          ...
         'anga'  ,          ...
         'angab' ,          ...
         'mcaa'  ,          ...
         'mcda'  ,          ...
         'drza'  ,          ...
         'draa'  ,          ...
         'stca'  ,          ...
         'pfsp'  ,          ...
         'pfsi'  ,          ...
         'pfed'  ,          ...
         'pfmx'  ,          ...
         'pfui'  ,          ...
         'pfmd'  ,          ...     
         'unitSessionMap',  ...
         'ucount'        ,  ...
         'pfsTag'        ,  ...
         'subjectId');
else,
    load(dataFilePath);
end




% $$$ pfhr = cat(1,pfhr{:});
% $$$ phzv = cat(1,phzv{:});
% $$$ vxya = cat(1,vxya{:});
% $$$ drza = cat(1,drza{:});
% $$$ draa = cat(1,draa{:});
% $$$ stca = cat(1,stca{:});
% $$$ pfsp = cat(1,pfsp{:});
% $$$ pfsi = cat(1,pfsi{:});


bins = linspace(-400,400,15);
binc = (bins(1:end-1)+bins(2:end))./2;

figure();
hind = discretize(pfhr,bins);
ny = size(stca,2);
for s = 1:ny,
ind = stca(:,1)==1                              ...
      &stca(:,s)==s                             ...
      &nniz(phzv)                               ...      
      &nniz(hind)                               ...
      &pfed(pfui)>25                            ...          
      &pfsi(pfui)>1.1                           ...          
      &pfsp(pfui)>0.1                           ...                
      &pfsp(pfui)<0.5                           ...                      
      &pfmd(pfui)<405;
    
    sind = sign(drza)==-1 & ind;

% INBOUND 

    subplot2(ny,7,s,1);
    out = accumarray(hind(sind,:),                                ... subs
                     sind(sind),                                         ... vals
                     repmat(numel(binc),[1,2]),                    ... size
                     @sum);                                    % func
    out(out==0) = nan;                       
    imagescnan({binc,binc,out'},[0,200],'linear',true,'colorMap',@jet);axis('xy');
    ylabel(states{s});
    subplot2(ny,7,s,2);
    scatter(pfhr(sind,1),pfhr(sind,2),10,phzv(sind),'filled');
    colormap(gca,'hsv');
    set(gca,'Color',[0.75,0.75,0.75]);
    xlim(bins([1,end]));ylim(bins([1,end]));
    
    subplot2(ny,7,s,3);    
    outIn = accumarray(hind(sind,:),                                ... subs
                       phzv(sind),                                  ... vals
                       repmat(numel(binc),[1,2]),                    ... size
                       @circ_mean);                                    % func
    outIn(outIn==0) = nan;                   
    imagescnan({binc,binc,outIn'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
    subplot2(ny,7,s,4);
    outIn = accumarray(hind(ind,:),                                ... subs
                       phzv(ind),                                  ... vals
                       repmat(numel(binc),[1,2]),                    ... size
                       @circ_mean);                                    % func
    outIn(outIn==0) = nan;                   
    imagescnan({binc,binc,outIn'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
    
    sind = sign(drza)==1 & ind;    
    subplot2(ny,7,s,5);
    mphz = accumarray(hind(sind,:),                                ... subs
                        phzv(sind),                                  ... vals
                        repmat(numel(binc),[1,2]),                    ... size
                        @circ_mean);                                    % func
    mphz(mphz==0) = nan;                   
    imagescnan({binc,binc,mphz'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
    subplot2(ny,7,s,6);
    scatter(pfhr(sind,1),pfhr(sind,2),10,phzv(sind),'filled');
    colormap(gca,'hsv');
    set(gca,'Color',[0.75,0.75,0.75]);
    xlim(bins([1,end]));ylim(bins([1,end]));
    subplot2(ny,7,s,7);
    mphz = accumarray(hind(sind,:),                                ... subs
                     sind(sind),                                         ... vals
                     repmat(numel(binc),[1,2]),                    ... size
                     @sum);                                    % func
    mphz(mphz==0) = nan;                   
    imagescnan({binc,binc,mphz'}, [0,100],'linear',true,'colorMap',@jet);axis('xy');

end



figure()
bins = linspace(-400,400,20);
binc = (bins(1:end-1)+bins(2:end))./2;

gridBins = cell([1,2]);
[gridBins{:}] = ndgrid(binc,binc);
gridBins = reshape(cat(3,gridBins{:}),[],2);

hind = discretize(pfhr,bins);
binPhz = linspace(-pi,pi,21);
ind = stca(:,1)==1                              ...
      &(stca(:,4)==4|stca(:,4)==4|stca(:,4)==4) ...
      &nniz(phzv)                               ...      
      &nniz(hind)                               ...      
      &anga<-0.1                                ...            
      &pfed(pfui)>25                            ...          
      &pfsi(pfui)>1.1                           ...          
      &pfsp(pfui)>0.1                           ...                
      &pfsp(pfui)<0.5                           ...                      
      &pfmd(pfui)<380;

%      &vxya(:,2)>1.2;
% $$$       &pfsi>1.1                                 ...
% $$$       &pfsp>0.1;

sp = tight_subplot(4,numel(binPhz)-1,0.001,0.1,0.01);
g = fittype( @(A,xa,ya,xya,xo,yo,x,y) A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
             'independent',{'x', 'y'},'dependent', 'z' ); 
fo = {};
for p = 1:numel(binPhz)-1;
    axes(sp(p));
    sind = ind&(binPhz(p)<phzv&phzv<=binPhz(p+1));
    scatter(pfhr(sind,1),pfhr(sind,2),2,phzv(sind),'filled');        
    colormap('hsv');
    caxis([-pi,pi]);
    xlim(bins([1,end]));ylim(bins([1,end]));    
    clear_axes_labels(gca);
    grid('on');    
    axes(sp(p+numel(binPhz)-1));    
    outIn = accumarray(hind(sind,:),                               ... subs
                       sind(sind),                                 ... vals
                       repmat(numel(binc),[1,2]),                  ... size
                       @sum);                                    % func
    outInImage = outIn;
    outInImage(outInImage==0) = nan;                   
    imagescnan({binc,binc,outInImage'./sum(outInImage(:),'omitnan')}, [0,0.03],'colorMap',@parula);axis('xy');
    clear_axes_labels(gca);
    grid('on');
    axes(sp(p+2*(numel(binPhz)-1)));
    fo{p} = fit(gridBins,outIn(:),g,'StartPoint',[100,0.00001,0.0000001,0.00001,0,0]);axis('xy');
    outModel = reshape(fo{p}(gridBins(:,1),gridBins(:,2)),size(outIn));
    imagescnan({binc,binc,outModel'});axis('xy');
    grid('on');    
    axes(sp(p+3*(numel(binPhz)-1)));
    imagescnan({binc,binc,(outIn-outModel)'},[0,20],'colorMap',@jet);axis('xy');
    grid('on');    
end


figure,
hold('on');
scatter(cell2mat(cf(@(f) f.xo, fo)),(binPhz(1:end-1)+binPhz(2:end))./2/pi*180,cell2mat(cf(@(f) f.xa, fo)).*1e6,'g','filled');
scatter(cell2mat(cf(@(f) f.xo, fo)),(binPhz(1:end-1)+binPhz(2:end))./2/pi*180+360,cell2mat(cf(@(f) f.xa, fo)).*1e6,'g','filled');
xlim([-100,100]);
Lines([],360,'k');



sind = ind&abs(pfhr(:,2))<50;
nbins = [20,50];
figure,
imagesc(linspace(-1,1,nbins(1)),...
        linspace(-3,2.5,nbins(2)),... 
        accumarray([discretize(drza(sind),linspace(-1,1,nbins(1))),...
                    discretize(vxya(sind,2),linspace(-3,2.5,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_mean)');
colormap(gca,'hsv');
colorbar();
axis('xy');

%sind = ind;
sind = ind&abs(pfhr(:,2))>150;
nbins = [75,50];
figure,
subplot(121);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-3,2.5,nbins(2)),...
        accumarray([discretize(pfhr(sind,1),linspace(-500,500,nbins(1))),...
                    discretize(vxya(sind,2),linspace(-3,2.5,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_mean)');
colormap(gca,'hsv');
colorbar();
axis('xy');
subplot(122);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-3,2.5,nbins(2)),...
        accumarray([discretize(pfhr(sind,1),linspace(-500,500,nbins(1))),...
                    discretize(vxya(sind,2),linspace(-3,2.5,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_var)');
colormap(gca,'jet');
colorbar();
axis('xy');



%sind = ind;
sind = ind&abs(pfhr(:,1))<50;
nbins = [75,50];
figure,
subplot(121);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-3,2.5,nbins(2)),...
        accumarray([discretize(pfhr(sind,2),linspace(-500,500,nbins(1))),...
                    discretize(vxya(sind,2),linspace(-3,2.5,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_mean)');
colormap(gca,'hsv');
colorbar();
axis('xy');
subplot(122);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-3,2.5,nbins(2)),...
        accumarray([discretize(pfhr(sind,2),linspace(-500,500,nbins(1))),...
                    discretize(vxya(sind,2),linspace(-3,2.5,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_var)');
colormap(gca,'jet');
colorbar();
axis('xy');


%sind = ind;
sind = ind&abs(pfhr(:,1))<150;
nbins = [75,50];
figure,
subplot(121);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-1,1,nbins(2)),...
        accumarray([discretize(pfhr(sind,2),linspace(-500,500,nbins(1))),...
                    discretize(drza(sind),linspace(-1,1,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_mean)');
colormap(gca,'hsv');
colorbar();
axis('xy');
subplot(122);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-1,1,nbins(2)),...
        accumarray([discretize(pfhr(sind,2),linspace(-500,500,nbins(1))),...
                    discretize(drza(sind),linspace(-1,1,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_var)');
colormap(gca,'jet');
colorbar();
axis('xy');


sind = ind;
%sind = ind&abs(pfhr(:,1))<100;
nbins = [75,50];
figure,
subplot(121);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-500,500,nbins(2)),...
        accumarray([discretize(pfhr(sind,1),linspace(-500,500,nbins(1))),...
                    discretize(pfhr(sind,2),linspace(-500,500,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_mean)');
colormap(gca,'hsv');
colorbar();
axis('xy');
subplot(122);
imagesc(linspace(-500,500,nbins(1)),...
        linspace(-500,500,nbins(2)),...
        accumarray([discretize(pfhr(sind,1),linspace(-500,500,nbins(1))),...
                    discretize(pfhr(sind,2),linspace(-500,500,nbins(2)))],...
                    phzv(sind),nbins-1,@circ_var)');
colormap(gca,'jet');
colorbar();
axis('xy');



figure();
ind = stca(:,1)==1                              ...
      &(stca(:,3)==3|stca(:,4)==4|stca(:,5)==5) ...
      &nniz(phzv)                               ...
      &nniz(hind)                               ...
      &pfsi>1.1                                 ...
      &pfsp>0.1;
sind = ind;
%sind = ind&(binPhz(p)<phzv&phzv<binPhz(p+1));
draat = draa;
% $$$ draat(draat<-pi/2) = draat(draat<-pi/2)+pi/2;
% $$$ draat(draat>pi/2) = draat(draat>pi/2)-pi/2;
subplot(331);
hold('on')
hax = scatter(drza(sind,1),phzv(sind),2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(drza(sind,1),phzv(sind)+2*pi,2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-1,1]);ylim([-pi,3*pi]);    
subplot(332);
hold('on')
hax = scatter(pfhr(sind,1),phzv(sind),2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,1),phzv(sind)+2*pi,2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-200,200]);ylim([-pi,3*pi]);    
subplot(333);
hold('on')
hax = scatter(pfhr(sind,2),phzv(sind),2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,2),phzv(sind)+2*pi,2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-200,200]);ylim([-pi,3*pi]);    
subplot(334);
hold('on');
hax = scatter(drza(sind,1),phzv(sind),2,'filled');
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(drza(sind,1),phzv(sind)+2*pi,2,'filled');
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
xlim([-1,1]);ylim([-pi,3*pi]);    
subplot(335);
hold('on');
hax = scatter(pfhr(sind,1),phzv(sind),2,'filled');
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,1),phzv(sind)+2*pi,2,'filled'); 
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
xlim([-200,200]);ylim([-pi,3*pi]);     
subplot(336);
hold('on');
hax = scatter(pfhr(sind,2),phzv(sind),2,'filled');
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,2),phzv(sind)+2*pi,2,'filled');
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
xlim([-200,200]);ylim([-pi,3*pi]);    
subplot(337);
hold('on');
hax = scatter(drza(sind,1),phzv(sind),2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(drza(sind,1),phzv(sind)+2*pi,2,draa(sind,1),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-1,1]);ylim([-pi,3*pi]);    
title('DRZ x Theta phase, c:drzang')
subplot(338);
hold('on');
hax = scatter(pfhr(sind,1),phzv(sind),2,atan2(pfhr(sind,2),pfhr(sind,1)),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,1),phzv(sind)+2*pi,2,atan2(pfhr(sind,2),pfhr(sind,1)),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-200,200]);ylim([-pi,3*pi]);    
title('HRF fwd x Theta phase, c:hrfang')
subplot(339);
hold('on');
hax = scatter(pfhr(sind,2),phzv(sind),2,atan2(pfhr(sind,2),pfhr(sind,1)),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
hax = scatter(pfhr(sind,2),phzv(sind)+2*pi,2,atan2(pfhr(sind,2),pfhr(sind,1)),'filled');        
hax.MarkerFaceAlpha = 0.5;
hax.MarkerEdgeAlpha = 0.5;
colormap('hsv');
caxis([-pi,pi]);
xlim([-200,200]);ylim([-pi,3*pi]);    
title('HRF lat x Theta phase, c:hrfang')



v = 1;
fitRes = 10000;
fitRange = linspace(-2*pi,2*pi,fitRes);

devars = {drza(ind,1),pfhr(ind,1),-abs(drza(ind)).*sign(pfhr(ind,1))};
% COMPUTE projection of the decoded position on longitudinal axis of the head
for v = 1:3
    decError = devars{v};
    decPhase = phzv(ind,1);

    % ASSIGN inputs
    lin = decError(:,1);
    circ = decPhase;
    % TRANSFORM linear and circular components into cos and sin parts
    cosPart = sum(cos(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
    sinPart = sum(sin(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
    R = sqrt((cosPart./length(circ)).^2 + ...
             (sinPart./length(circ)).^2 );
    [lmi,lmv] = LocalMinima(-R',round(fitRes/10),0,1);
    % SAVE model parameters
    modelTraj.R(v) = R(lmi); % R: fit quality 
    modelTraj.parameters(v,:) = [2*pi*fitRange(lmi),... %Regression Parm : Slope
                        atan2(sinPart(lmi),cosPart(lmi))];%: Offset
                                                          % COMPUTE circular-linear correlation coefficient
    linC = mod(abs(modelTraj.parameters(v,1))*lin,2*pi);
    circMean = atan2(sum(sin(circ)),sum(cos(circ)));
    linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
    modelTraj.rho(v) = sum(sin(circ-circMean).*sin(linC-linCMean))...
        ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
end




pfchrI = cat(1,pfchrI{:});
pfchrO = cat(1,pfchrO{:});
vxyO =  cat(1,vxyO{:});
vxyI =  cat(1,vxyI{:});
hindO = cat(1,hindO{:});
hindI = cat(1,hindI{:});
phzvO = cat(1,phzvO{:});
phzvI = cat(1,phzvI{:});

    

nind = nniz(hindI)&nniz(phzvI);
outIn = accumarray(hindI(nind,:),                                ... subs
                   phzvI(nind),                                  ... vals
                   repmat(numel(binc),[1,2]),                    ... size
                   @circ_mean);                                    % func
outIn(outIn==0) = nan;                   
nind = nniz(hindO)&nniz(phzvO);
outOut = accumarray(hindO(nind,:),                               ... subs
                    phzvO(nind),                                 ... vals
                    repmat(numel(binc),[1,2]),                   ... size
                    @circ_mean);                                   % func
outOut(outOut==0) = nan;

subplot(223);
imagescnan({binc,binc,outIn'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
subplot(224);
imagescnan({binc,binc,outOut'},[-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');

nindO = nniz(hindO)&nniz(phzvO);
nindI = nniz(hindI)&nniz(phzvI);
outAll = accumarray([hindO(nindO,:);hindI(nindI,:)],               ... subs
                    [phzvO(nindO);phzvI(nindI)],                   ... vals
                    repmat(numel(binc),[1,2]),                   ... size
                    @circ_mean);                                   % func
outAllS = accumarray([hindO(nindO,:);hindI(nindI,:)],               ... subs
                    [phzvO(nindO);phzvI(nindI)],                   ... vals
                    repmat(numel(binc),[1,2]),                   ... size
                    @circ_var);                                   % func
figure();
subplot(121);
imagescnan({binc,binc,outAll'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
subplot(122);
imagescnan({binc,binc,outAllS'}, [0,pi/2],'circular',true,'colorMap',@jet);axis('xy');


hindO = discretize(pfchrO,bins);
hindI = discretize(pfchrI,bins);

figure();
v = 1;
vbins = prctile([vxyO(vxyO(:,v)>-3,v);vxyI(vxyI(:,v)>-3,v)],linspace(0.001,99.999,7))';
vbinc = (vbins(1:end-1)+vbins(2:end))./2;
vindO = discretize(vxyO(:,v),vbins);
vindI = discretize(vxyI(:,v),vbins);
ny = numel(vbinc);
for vv = 1:ny,
    nindI = nniz(vindI)&nniz(pfchrI)&nniz(phzvI)&nniz(hindI);
    nindO = nniz(vindO)&nniz(pfchrO)&nniz(phzvO)&nniz(hindO);
% INBOUND 
    sind = nindI&(vindI==vv);
    subplot2(ny,6,ny+1-vv,1);
    out = accumarray(hindI(sind,:),                                ... subs
                     sind(sind),                                         ... vals
                     repmat(numel(binc),[1,2]),                    ... size
                     @sum);                                    % func
    out(out==0) = nan;                   
    imagescnan({binc,binc,out'},[0,100],'linear',true,'colorMap',@jet);axis('xy');

    subplot2(ny,6,ny+1-vv,2);
    scatter(pfchrI(sind,1),pfchrI(sind,2),10,phzvI(sind),'filled');
    colormap(gca,'hsv');
    set(gca,'Color',[0.75,0.75,0.75]);
    xlim([-300,300]);ylim([-300,300]);
    
    subplot2(ny,6,ny+1-vv,3);
    outIn = accumarray(hindI(sind,:),                                ... subs
                       phzvI(sind),                                  ... vals
                       repmat(numel(binc),[1,2]),                    ... size
                       @circ_mean);                                    % func
    outIn(outIn==0) = nan;                   
    imagescnan({binc,binc,outIn'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
    sind = nindO&(vindO==vv);    
    subplot2(ny,6,ny+1-vv,4);
    outOut = accumarray(hindO(sind,:),                                ... subs
                        phzvO(sind),                                  ... vals
                        repmat(numel(binc),[1,2]),                    ... size
                        @circ_mean);                                    % func
    outOut(outOut==0) = nan;                   
    imagescnan({binc,binc,outOut'}, [-pi,pi],'circular',true,'colorMap',@hsv);axis('xy');
    subplot2(ny,6,ny+1-vv,5);
    scatter(pfchrO(sind,1),pfchrO(sind,2),10,phzvO(sind),'filled');
    colormap(gca,'hsv');
    set(gca,'Color',[0.75,0.75,0.75]);
    xlim([-300,300]);ylim([-300,300]);
    subplot2(ny,6,ny+1-vv,6);
    out = accumarray(hindO(sind,:),                                ... subs
                     sind(sind),                                         ... vals
                     repmat(numel(binc),[1,2]),                    ... size
                     @sum);                                    % func
    out(out==0) = nan;                   
    imagescnan({binc,binc,out'}, [0,100],'linear',true,'colorMap',@jet);axis('xy');

end


% COMPUTE circular-linear correlation
% $$$ clear('modelTraj'); 
% $$$ modelTraj.partitions = prctile(vxy(sind,speedInd),linspace(0,100,7));
% $$$ % VISUALDOC partitioning of head speed space
% $$$ % $$$     figure,
% $$$ % $$$     hist(vxy(sind,speedInd),1000)
% $$$ % $$$     Lines(modelTraj.partitions,[],'r');
% $$$ % SET fiting resolution and range
% $$$ fitRes = 10000;
% $$$ fitRange = linspace(-2*pi,2*pi,fitRes);
% $$$ % FIT circular linear model for each speed partition 
% $$$ for v = 1:numel(modelTraj.partitions)-1,
% $$$     ind = cind  &  sind      ...
% $$$           &  vxy(:,speedInd)>modelTraj.partitions(v)       ...
% $$$           &  vxy(:,speedInd)<modelTraj.partitions(v+1);
% $$$     % COMPUTE projection of the decoded position on longitudinal axis of the head
% $$$     decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
% $$$     decPhase = phz(ind,1);
% $$$     % ASSIGN inputs
% $$$     lin = decError(:,1);
% $$$     circ = decPhase;
% $$$     % TRANSFORM linear and circular components into cos and sin parts
% $$$     cosPart = sum(cos(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
% $$$     sinPart = sum(sin(bsxfun(@minus,circ,2*pi*lin*fitRange)),1);
% $$$     R = sqrt((cosPart./length(circ)).^2 + ...
% $$$              (sinPart./length(circ)).^2 );
% $$$     [lmi,lmv] = LocalMinima(-R',round(fitRes/10),0,1);
% $$$     % SAVE model parameters
% $$$     modelTraj.R(v) = R(lmi); % R: fit quality 
% $$$     modelTraj.parameters(v,:) = [2*pi*fitRange(lmi),... %Regression Parm : Slope
% $$$                         atan2(sinPart(lmi),cosPart(lmi))];%: Offset
% $$$                                                           % COMPUTE circular-linear correlation coefficient
% $$$     linC = mod(abs(modelTraj.parameters(v,1))*lin,2*pi);
% $$$     circMean = atan2(sum(sin(circ)),sum(cos(circ)));
% $$$     linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
% $$$     modelTraj.rho(v) = sum(sin(circ-circMean).*sin(linC-linCMean))...
% $$$         ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
% $$$ end
% $$$ 



ufr = load(Trial,'ufr',xyz,[],unitSubset,0.240,true,'gauss');    

[mins,mvals] = LocalMinima(-sum(ufr.data(:,5),2),10,-0.5);
mvals = mvals(stcm(mins,1)==1);
mins = mins(stcm(mins,1)==1);


figure,
plot(vxy(mins,2),-mvals,'.')




                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
% Subject Population egocentric phase precession
hfigPop = figure();
hfig = figure();
hfig.Units = 'Centimeters';
hfig.Position = [0,0,40,10];


bins = linspace(-400,400,20);
binc = (bins(1:end-1)+bins(2:end))./2;

gridBins = cell([1,2]);
[gridBins{:}] = ndgrid(binc,binc);
gridBins = reshape(cat(3,gridBins{:}),[],2);

hind = discretize(pfhr,bins);
binPhz = linspace(-pi,pi,21);
binPhzc = (binPhz(1:end-1)+binPhz(2:end))./2;

g = fittype( @(A,xa,ya,xya,xo,yo,x,y) A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
             'independent',{'x', 'y'},'dependent', 'z' ); 

hfigPopTitle = 'Egocentric Phase Precession: General Behaviors';
stateColors = 'krbgm';
hfigPopTag = 'general'; 
shift = 0;
% $$$ 
% $$$ stateColors = 'rmcb';
% $$$ hfigPopTitle = 'Egocentric Phase Precession: Behavior Subtypes';
% $$$ hfigPopTag = 'subtypes'; 
% $$$ shift = 5;

clf(hfigPop);
for s = 1:numel(stateColors);
    datFileName = ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
                   ['egocentric_pp_fit_',pfsTag,subjectId,'_',states{s+shift},'.mat']];

    if 1,%~exist(datFileName,'file'),
        ind = stca(:,1)==1                              ...
              &(stca(:,s+shift)==s+shift)                           ...
              &nniz(phzv)                               ...      
              &nniz(hind)                               ...      
              &pfed(pfui)>20                            ...          
              &pfsi(pfui)>1.1                           ...          
              &pfsp(pfui)>0.1                           ...                
              &pfsp(pfui)<0.5                           ...                      
              &pfmd(pfui)<405;

        figure(hfig);
        clf();
        sp = tight_subplot(4,numel(binPhz)-1,0.001,0.1,0.05);
        fo = {};
        for p = 1:numel(binPhz)-1;
            axes(sp(p));
            sind = ind&(binPhz(p)<phzv&phzv<=binPhz(p+1));
            scatter(pfhr(sind,1),pfhr(sind,2),2,phzv(sind),'filled');        
            colormap('hsv');
            caxis([-pi,pi]);
            xlim(bins([1,end]));ylim(bins([1,end]));    
            clear_axes_labels(gca);
            grid('on');    
            Lines(0,[],'r');      
            if p==1,
                ylabel({'pfs center','given spike'}); 
            end
            % PLOT JPDF of egocentric head position
            axes(sp(p+numel(binPhz)-1));    
            outIn = accumarray(hind(sind,:),                               ... subs
            sind(sind),                                 ... vals
            repmat(numel(binc),[1,2]),                  ... size
            @sum);                                    % func
            outInImage = outIn;
            outInImage(outInImage==0) = nan;                   
            imagescnan({binc,binc,outInImage'./sum(outInImage(:),'omitnan')}, [0,0.03],'colorMap',@parula);axis('xy');
            clear_axes_labels(gca);
            grid('on');
            Lines(0,[],'r');            
            if p==1,ylabel({'JPDF','ca:[0,0.03]'});end
            % FIT gaussian to phase bin
            axes(sp(p+2*(numel(binPhz)-1)));
            fo{p} = fit(gridBins,outIn(:),g,'StartPoint',[100,0.00001,0.0000001,0.00001,0,0]);axis('xy');
            outModel = reshape(fo{p}(gridBins(:,1),gridBins(:,2)),size(outIn));
            imagescnan({binc,binc,(outModel./sum(outModel(:),'omitnan'))'},[0,0.03]);axis('xy');
            clear_axes_labels(gca);
            grid('on');
            Lines(0,[],'r');            
            if p==1,ylabel({'Gaussian','Fit'});end
            % PLOT residuals
            axes(sp(p+3*(numel(binPhz)-1)));
            imagescnan({binc,binc,(outIn-outModel)'},[0,20],'colorMap',@jet);axis('xy');
            clear_axes_labels(gca);
            grid('on');    
            Lines(0,[],'r');
            if p==1,ylabel({'Residuals','ca:[0,20]'});end
        end
        
        fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
        xlim([0,1]);
        ylim([0,1]);
        line([0.05,0.05+0.04405], [0.08,0.08]);


        save(datFileName,...
             'fo','binPhz','bins','binPhzc','binc');
        print(hfig,'-depsc2',...
              ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
               ['egocentric_pp_fit_',pfsTag,subjectId,'_',states{s+shift},'.eps']]);
        print(hfig,'-dpng',...
              ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
               ['egocentric_pp_fit_',pfsTag,subjectId,'_',states{s+shift},'.png']]);       
        
    else
        load(datFileName);
    end
    
        

    figure(hfigPop);
    hold('on');
    scatter(repmat(cell2mat(cf(@(f) f.xo, fo)),[1,2]),...
            [binPhzc/pi*180,binPhzc/pi*180+360],...
            repmat(cell2mat(cf(@(f) f.xa, fo)).*1e6,[1,2]),...
                   stateColors(s),'filled');

end
legend(states(shift+[1:numel(stateColors)]));
xlim([-100,100]);
Lines([],360,'k');
Lines([],0,'k');    
Lines([],180,'k','--');        


title({hfigPopTitle,'mean placefield position projected onto the rostro-caudal axis of the head'})

ylabel('CA1pyr Theta Phase (degrees)');
xlabel('Position (mm)');


print(hfigPop,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['egocentric_pp_fit_',pfsTag,subjectId,'_',hfigPopTag,'.eps']]);
print(hfigPop,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['egocentric_pp_fit_',pfsTag,subjectId,'_',hfigPopTag,'.png']]);       


                     
                     
                      
                     
                     
                     

                     
                     
% Subject Population egocentric phase precession binned by speed
hfigPop = figure();
hfig = figure();
hfig.Units = 'Centimeters';
hfig.Position = [0,0,40,10];


bins = linspace(-400,400,20);
binc = (bins(1:end-1)+bins(2:end))./2;

gridBins = cell([1,2]);
[gridBins{:}] = ndgrid(binc,binc);
gridBins = reshape(cat(3,gridBins{:}),[],2);

hind = discretize(pfhr,bins);
binPhz = linspace(-pi,pi,31);
binPhzc = (binPhz(1:end-1)+binPhz(2:end))./2;

g = fittype( @(A,xa,ya,xya,xo,yo,x,y) A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
             'independent',{'x', 'y'},'dependent', 'z' ); 

hfigPopTitle = 'Egocentric Phase Precession: speed';



clf(hfigPop);
ind = stca(:,1)==1                              ...
      &(stca(:,3)==3|stca(:,4)==4|stca(:,5)==5) ...
      &nniz(phzv)                               ...      
      &nniz(hind)                               ...      
      &pfed(pfui)>20                            ...          
      &pfsi(pfui)>1.1                           ...          
      &pfsp(pfui)>0.1                           ...                
      &pfsp(pfui)<0.5                           ...                      
      &pfmd(pfui)<405;

speedBins = [prctile(vxya(ind&vxya(:,2)>-2&vxya(:,2)<1.9,2),linspace(0,100,9))];
speedBinc = (speedBins(1:end-1)+speedBins(2:end))./2;
speedColors = jet(numel(speedBinc));
vind = discretize(vxya(:,2),speedBins); 


for s = 1:size(speedColors,1);
    datFileName = ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
                   ['egocentric_pp_fit_',pfsTag,subjectId,'_speed-',num2str(round(speedBinc(s),3)),'.mat']];
    if 0%~exist(datFileName,'file'),
        figure(hfig);
        clf();
        sp = tight_subplot(4,numel(binPhz)-1,0.001,0.1,0.01);
        fo = {};
        for p = 1:numel(binPhz)-1;
            axes(sp(p));
            sind = ind&(binPhz(p)<phzv&phzv<=binPhz(p+1))&vind==s;
            scatter(pfhr(sind,1),pfhr(sind,2),2,phzv(sind),'filled');        
            colormap('hsv');
            caxis([-pi,pi]);
            xlim(bins([1,end]));ylim(bins([1,end]));    
            clear_axes_labels(gca);
            grid('on');    
            % PLOT JPDF of egocentric head position
            axes(sp(p+numel(binPhz)-1));    
            outIn = accumarray(hind(sind,:),                               ... subs
            sind(sind),                                 ... vals
            repmat(numel(binc),[1,2]),                  ... size
            @sum);                                    % func
            outInImage = outIn;
            outInImage(outInImage==0) = nan;                   
            imagescnan({binc,binc,outInImage'./sum(outInImage(:),'omitnan')}, [0,0.03],'colorMap',@parula);axis('xy');
            clear_axes_labels(gca);
            grid('on');
            % FIT gaussian to phase bin
            axes(sp(p+2*(numel(binPhz)-1)));
            fo{p} = fit(gridBins,outIn(:),g,'StartPoint',[100,0.00001,0.0000001,0.00001,0,0]);axis('xy');
            outModel = reshape(fo{p}(gridBins(:,1),gridBins(:,2)),size(outIn));
            imagescnan({binc,binc,outModel'});axis('xy');
            clear_axes_labels(gca);
            grid('on');
            % PLOT residuals
            axes(sp(p+3*(numel(binPhz)-1)));
            imagescnan({binc,binc,(outIn-outModel)'},[0,20],'colorMap',@jet);axis('xy');
            clear_axes_labels(gca);
            grid('on');    
        end
% $$$         fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
% $$$         xlim([0,1]);
% $$$         ylim([0,1]);
% $$$         Lines(hfig.Position(3)/2+0.2,[],'k');
        
        save(datFileName,...
             'fo','binPhz','bins','binPhzc','binc');
        print(hfig,'-depsc2',...
              ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
               ['egocentric_pp_fit_',pfsTag,subjectId,'_speed-',num2str(round(speedBinc(s),3)),'.eps']]);
        print(hfig,'-dpng',...
              ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
               ['egocentric_pp_fit_',pfsTag,subjectId,'_speed-',num2str(round(speedBinc(s),3)),'.png']]);       
        
    else
        load(datFileName);
    end
    
        
    figure(hfigPop);
    hold('on');
    scatter(repmat(cell2mat(cf(@(f) f.xo, fo)),[1,2]),...
            [binPhzc/pi*180,binPhzc/pi*180+360],...
            repmat(cell2mat(cf(@(f) f.xa, fo)).*1e6/2,[1,2]),...
                   repmat(speedColors(s,:),[numel(binPhzc).*2,1]),'filled');
end

vlabel = {};
for s = 1:numel(speedBinc),
    vlabel{s} = num2str(median(10.^nonzeros(vxya(ind&vind==s,2))));
end

set(gca,'XTick',-50:10:100)
grid('on');

%legend(cf(@num2str,num2cell(speedBinc)));
%legend(cf(@num2str,num2cell(10.^speedBinc)));
legend(vlabel);
xlim([-50,100]);
Lines([],360,'k');
Lines([],0,'k');    
Lines([],180,'k','--');        

title({hfigPopTitle,'mean placefield position projected onto the rostro-caudal axis of the head'})

ylabel('CA1pyr Theta Phase (degrees)');
xlabel('Position (mm)');

print(hfigPop,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['egocentric_pp_fit_',pfsTag,subjectId,'_speed.eps']]);
print(hfigPop,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['egocentric_pp_fit_',pfsTag,subjectId,'_speed.png']]);       







                     
                     
% Subject Population egocentric phase precession binned by angular speed and speed
plotFitFlag = false;

subjectId = 'er01';
subjectInds = [1:2];
phzBinCnt = 8;
spcBinCnt = 20;
angInd = 4;
angBinCnt = 0;
angBins = [-3,-0.7,-0.25,0,0.25,0.7,3];
% $$$ angBins = [-3,-0.75,-0.20,-0.05,0,0.05,0.20,0.75,3];
% $$$ angColors = [0,0,1;0,0,0.85;0,0,0.75;0,0,0.5;0.5,0,0;0.75,0,0;0.85,0,0;1,0,0];


subjectId = 'ER06';
subjectInds = [3:5];
phzBinCnt = 10;
spcBinCnt = 20;
angInd = 4;
angBinCnt = 0;
angBins = [-2.99,-0.57,-0.14,0,0.14,0.57,2.99];


 
subjectId = 'Ed10';
subjectInds = [6:7];
phzBinCnt = 8;
spcBinCnt = 20;
angInd = 4;
angBinCnt = 0;
angBins = [-2.99,-0.57,-0.14,0,0.14,0.57,2.99];


subjectId = 'jg04';
subjectInds = [8:12];
subjectId = 'jg04';
subjectInds = [13:16];
phzBinCnt = 8;
spcBinCnt = 20;
angInd = 4;
angBinCnt = 0;
angBins = [-2.99,-0.57,-0.14,0,0.14,0.57,2.99];

subjectId = 'jg05';
subjectInds = [17:23];
angBins = [-2.99,-0.57,-0.14,0,0.14,0.57,2.99];
%angBins = [-2.99,-1.13,-0.74,-0.51,-0.34,-0.21,-0.11,-0.03,0,0.03,0.11,0.21,0.34,0.51,0.74,1.13,2.99];
angBinCnt = 2;
spcBinCnt = 20;
phzBinCnt = 20;
spdBinCnt = 12;
spdInd = 2;
colorMode = 'speed';
colorMap = @jet;
speedTags = {'body speed','head speed'};
%colorMode = 'angvel';
%colorMap = @cool;
partitionTag = 'Partitioned by:';


% $$$ stateInd = [4];
stateInd = [5];
stateInd = [3,4,5];

thresholds.mazeCenterDist = 400;
thresholds.mazeCenterAng = pi/2;


hfigPop = figure();
hfig = figure();
hfig.Units = 'Centimeters';
hfig.Position = [0,0,40,10];


% SET bins edges of 2d space
% DISCRETIZE egocentric space
bins = linspace(-400,400,spcBinCnt);
binc = (bins(1:end-1)+bins(2:end))./2;
gridBins = cell([1,2]);
[gridBins{:}] = ndgrid(binc,binc);
gridBins = reshape(cat(3,gridBins{:}),[],2);
hind = discretize(pfhr,bins);

% SET bins edges of theta phase
phzBins = linspace(-pi,pi,phzBinCnt+1);
phzBinc = (phzBins(1:end-1)+phzBins(2:end))./2;
if phzBinCnt>1,  partitionTag = [partitionTag,' theta phase;'];  end

hfigPopTitle = 'Egocentric Phase Precession: speed';

clf(hfigPop);


% SELECT spikes
cind =    mcda < thresholds.mazeCenterDist  ...
        | abs(mcaa) < thresholds.mazeCenterAng;
sind = any(logical(stca(:,stateInd)),2) & logical(stca(:,1)); % select spikes, theta and bhv state union
tind = ismember(sid(pfui),subjectInds);
pind =  pfed(pfui)>20                            ... select units by Edist, unit qualitiy
      & pfsp(pfui)>0.05                          ... select units by sparsity, placefileds size
      & pfsp(pfui)<0.4                           ... select units by sparsity, placefileds size
      & pfmd(pfui)<390;                            % select units by spatial position 


ind =   cind & tind & sind & pind & nniz(phzv) & nniz(hind);
        

% SET speed bins
% DISCRETIZE speed
spdBins = [prctile(vxya(ind&cind&vxya(:,spdInd)>-2&vxya(:,spdInd)<1.9,2),linspace(0,100,spdBinCnt+1))];
spdBinc = (spdBins(1:end-1)+spdBins(2:end))./2;
vind = discretize(vxya(:,spdInd),spdBins); 
if spdBinCnt>1,  partitionTag = [partitionTag,' ',speedTags{speedInd},';']; end

% SET angular speed bins
% DISCRETIZE Angular speed
if angBinCnt,
    angBins = reshape([prctile(anga(ind&anga(:,angInd)>-3&anga(:,angInd)<3,angInd),...
                               linspace(0,100,angBinCnt))],[],2);
    angBins = mean([abs(angBins(1:end,1)),flipud(angBins(1:end,2))]');
    angBins = [-angBins,0,fliplr(angBins)];
end
angBinc = (angBins(1:end-1)+angBins(2:end))./2;
angBinCnt = numel(angBinc);
aind = discretize(anga(:,angInd),angBins); 
if angBinCnt>1,  partitionTag = [partitionTag,' head angvel;'];  end

switch colorMode,
    case 'speed'
      markerColorMap = colorMap(numel(spdBinc));
      colorMapIteratorName = 's';
      colorMapTag = ['Color map: speed [',num2str(spdBins(1)),',',num2str(spdBins(end)),']'];
    case 'angvel'
      markerColorMap = colorMap(numel(angBinc));
      colorMapIteratorName = 'a';
      colorMapTag = ['Color map: angvel [',num2str(angBins(1)),',',num2str(angBins(end)),']'];
end

tag = DataHash({subjectId,stateInd,angInd,angBins,spdInd,spdBins,phzBins});
figPath = '/storage/share/Projects/BehaviorPlaceCode/phase_precession/';
fitModel = fittype( @(A,xa,ya,xya,xo,yo,x,y) ...
                    A.*exp(-(xa.*(x-xo).^2+xya.*(x-xo).*(y-yo)+ya.*(y-yo).^2)),...
                   'independent',{'x', 'y'},'dependent', 'z' ); 
fitStartPoint = [100,0.00001,0.0000001,0.00001,0,0];
for s = 1:numel(spdBinc),
    for a = 1:numel(angBinc),
        figName = ['egocentric_pp_fit_',tag];
        fitDataPath = [dataFilePath(1:end-4),'-',tag,'-',num2str(s),'-',num2str(a),'.mat'];
        if ~exist(fitDataPath,'file')|overwrite,
            figure(hfig);
            clf();
            sp = tight_subplot(4,numel(phzBins)-1,0.001,0.1,0.01);
            fo = {};
            for p = 1:numel(phzBins)-1;                
                % SELECT for phase
                sind = ind&(phzBins(p)<phzv&phzv<=phzBins(p+1))&vind==s&aind==a;
                % BIN spike egocentric placefield position
                outIn = accumarray(hind(sind,:),                               ... subs
                                   sind(sind),                                 ... vals
                                   repmat(numel(binc),[1,2]),                  ... size
                                   @sum);                                        % func
                % FIT gaussian 
                fo{p} = fit(gridBins,outIn(:),fitModel,'StartPoint',fitStartPoint);axis('xy');
                outModel = reshape(fo{p}(gridBins(:,1),gridBins(:,2)),size(outIn));

                % DISPLAY 
                if plotFitFlag
                    % PLOT spike egocentric placefield position
                    axes(sp(p));
                    scatter(pfhr(sind,1),pfhr(sind,2),2,phzv(sind),'filled');        
                    %scatter(drza(sind,1),abs(drza(sind,1)).*sin(draa(sind,1)),2,phzv(sind),'filled');
                    colormap('hsv');
                    caxis([-pi,pi]);
                    xlim(bins([1,end]));ylim(bins([1,end]));    
                    clear_axes_labels(gca);
                    grid('on');    
                    % PLOT JPDF of spike egocentric placefield position
                    axes(sp(p+numel(phzBins)-1));
                    outInImage = outIn;
                    outInImage(outInImage==0) = nan;                   
                    imagescnan({binc,binc,outInImage'./sum(outInImage(:),'omitnan')}, ...
                               [0,0.03],'colorMap',@parula);axis('xy');
                    clear_axes_labels(gca);
                    grid('on');
                    % PLOT Gaussian fit
                    axes(sp(p+2*(numel(phzBins)-1)));                
                    imagescnan({binc,binc,outModel'});axis('xy');
                    clear_axes_labels(gca);
                    grid('on');
                    % PLOT residuals
                    axes(sp(p+3*(numel(phzBins)-1)));
                    imagescnan({binc,binc,(outIn-outModel)'},[0,20],'colorMap',@jet);axis('xy');
                    clear_axes_labels(gca);
                    grid('on');    
                end
            end

            
            save(fitDataPath,'fo','phzBins','bins','phzBinc','binc');
            if plotFitFlag            
                fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
                xlim([0,1]);ylim([0,1]);
                
                %print(hfig,'-depsc2',fullfile(figPath,[figName,'.eps']));
                print(hfig,'-dpng',  fullfile(figPath,[figName,'.png']));
            end
            clear('fax');
        else
            load(fitDataPath);
        end
        
        figure(hfigPop);
        hold('on');
        sax = scatter3(repmat(cell2mat(cf(@(f) f.xo, fo)),[1,2]),...
                       repmat(cell2mat(cf(@(f) f.yo, fo)),[1,2]),...
                       [phzBinc/pi*180,phzBinc/pi*180+360],...
                       repmat(cell2mat(cf(@(f) f.xa, fo)).*1e6/2,[1,2]),...
                       repmat(markerColorMap(eval(colorMapIteratorName),:),[numel(phzBinc).*2,1]),'filled');
        % Shhhhh,... I don't EVER want to hear a word about the previous line.
    end
end
title({'Head Centered - Head Coordinates',partitionTag,colorMapTag});
zlabel('theta phase (deg)');
xlabel('logitudinal (mm)');
ylabel('lateral (mm)');
% $$$         Lines(hfig.Position(3)/2+0.2,[],'k');
% $$$ vlabel = {};
% $$$ for s = 1:numel(speedBinc),
% $$$     vlabel{s} = num2str(median(10.^nonzeros(vxya(ind&vind==s,2))));
% $$$ end
% $$$ 
% $$$ set(gca,'XTick',-50:10:100)
% $$$ grid('on');
% $$$ 
% $$$ %legend(cf(@num2str,num2cell(speedBinc)));
% $$$ %legend(cf(@num2str,num2cell(10.^speedBinc)));
% $$$ legend(vlabel);
% $$$ xlim([-50,100]);
% $$$ % $$$ Lines([],360,'k');
% $$$ % $$$ Lines([],0,'k');    
% $$$ % $$$ Lines([],180,'k','--');        
% $$$ 
% $$$ title({hfigPopTitle,'mean placefield position projected onto the rostro-caudal axis of the head'})
% $$$ 
% $$$ ylabel('CA1pyr Theta Phase (degrees)');
% $$$ xlabel('Position (mm)');
% $$$ 
% $$$ print(hfigPop,'-depsc2',...
% $$$       ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
% $$$        ['egocentric_pp_fit_',subjectId,'_speed.eps']]);
% $$$ print(hfigPop,'-dpng',...
% $$$       ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
% $$$        ['egocentric_pp_fit_',subjectId,'_speed.png']]);       


                     
                     
                     
                     
                     
                     

                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     
                     

                        
