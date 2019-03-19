% MjgER2016_figure5
% Behavior Theta Decomposition

% DECOMPOSITION of behavior space across theta phases


% FIGURE : Behavioral theta phase decomposition Examples -------------------------------------------%

% bhvMulti  mod 52, 60, 61, 80,
% bhvSingle mod 81, 86, 103, 104 
% bhvProg   mod 63, 34, 83, 72, 145

MjgER2016_load_data();
%  Variables:
%      Trials
%      units
%      cluSessionMap
%      pitchReferenceTrial
%      FigDir
%      sessionListName
%      sessionList
%      states
%      numStates
%      interpParPfsp
%      interpParDfs
%
%  Functions:
%      reshape_eigen_vector

t = 20;
Trial = Trials{t};
unitSubset = units{t};

% LOAD placefields and theta restricted behavior fields
[pft]          = pfs_2d_theta(Trial,unitSubset);
[pfd]          = req20180123_ver5(Trial);
[pfs,metaData] = req20181106(Trial,[],[],unitSubset);

% LOAD group factor analysis of theta-space restricted behavior fields
[~,tags,eigVec,eigVar,eigScore,validDims,...
 unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials,'loadPFDFlag',false);

% LOAD group factor analysis of theta-space theta decomposition of each unit
[eigVecs, eigScrs, eigVars, eSpi, FSrBhvThp, metaData] = req20181119();


MjgER2016_load_bhv_erpPCA_scores();
%     rmaps    - matrix[D x S](numeric); rate maps corresponding to the valid eigenvector dims
%     FSrC     - matrix[U x V](Numeric); fscores
%     fsrcz    - matrix[U x V](Numeric); normalized fscores
clear(varNonEss{:},varNonAux{:});


 %req20190122();% -> bhvInfoTS

unitsExampleSubset = [81,86,104];
unitsExampleSubset = [52,61,72];
unitsExampleSubset = [63, 34, 83];
unitsExampleSubset = [34, 83, 63, ... % rearing transitions
                      52, 61, 72, ... % nonselective transitions
                      81, 86,104];    % selective stationary

unitsExampleSubset = [34, 83, 63, ... % rearing transitions
                      110, 79, 111, ... % nonselective transitions
                      35,139,103];    % selective stationary


nx = pfs.adata.binSizes(2)/2+4;
ny = numel(unitsExampleSubset);


cond_round = @(rate) max([round(rate,0),round(rate,1)].*[rate>=10,rate<10]);
nanColor = [0.25,0.25,0.25];

% SET figure options
fig.page.units    = 'centimeters';
fig.page.width    = 21.0;
fig.page.height   = 29.7;
fig.page.PaperPositionMode = 'auto';
fig.page.marginLeft  = 2.5;
fig.page.marginTop   = 3.7;
fig.subplot.width    = 1.15;
fig.subplot.height   = 1.15;
fig.subplot.horizontalPadding = 0.0;
fig.subplot.verticalPadding   = 0.15;

% SETUP subplot grid positions
xpos = fig.page.marginLeft ...
       : (fig.subplot.width  + fig.subplot.horizontalPadding) ...
       : fig.page.width;
ypos = fliplr(0.5 ...
       : (fig.subplot.height + fig.subplot.verticalPadding) ...
       : fig.page.height - fig.page.marginTop);

% SETUP Main figure
hfig = figure(666005);
clf(hfig);
hfig.Units = fig.page.units;
hfig.Position = [1,                               ...
                 1,                               ...
                 fig.page.width,                  ...
                 fig.page.height];
hfig.PaperPositionMode = fig.page.PaperPositionMode;
sp = gobjects([0,1]);

% SETUP Background axes
sp = gobjects([1,0]);
fax = axes('Position',[0,0,1,1],'Visible','off','Units','centimeters');
xlim([0,hfig.Position(3)]);
ylim([0,hfig.Position(4)]);


% BEGIN Main figure
phzOrder = [9:16,1:8];
%phzOrder = [10,12,14,16,2,4,6,8];
phzLables = {'','90','','180','','270','',''};
clear('i');
for y = 1:ny,
    yind = y;

% GET space and theta phase jointly restricted bhv rate maps
    rmap = plot(pfs,unitsExampleSubset(y),1,'colorbar',[],false);
    rmapNind = sq(sum(reshape(~isnan(rmap(:)),size(rmap)),2))>=pfs.adata.binSizes(2);
    zmap = reshape(permute(rmap,[1,3,2]),[],pfs.adata.binSizes(2));
    zmap(repmat(~rmapNind(:),[1,pfs.adata.binSizes(2)])) = 0;            
    zmap(zmap(:)<0) = 0;            
    rmapNind = rmapNind&reshape(sum(reshape(zmap,[],pfs.adata.binSizes(2))','omitnan')~=0,...
                                size(rmap,1),size(rmap,3));      
    mrate= prctile(rmap(:),99.9);    
    
% PLOT theta restricted spatial rate maps
    xind = 1;
    sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),fig.subplot.width,fig.subplot.height],...
                 'FontSize', 8);
    plot(pft,unitsExampleSubset(y),'mean','text',mrate,'nanColor',nanColor,'colorMap',@jet);
    %ylabel(['Unit: ',num2str(unitsExampleSubset(y))]);
    ylabel([num2str(unitsExampleSubset(y))]);
    set(sp(end),'YTickLabels',{});
    set(sp(end),'XTickLabels',{});    
    if y == 1,
        title({'PFS'});
    end
    
% PLOT spatially restricted bhv rate maps
    xind = 2;
    sp(end+1) = axes('Units','centimeters',...
                 'Position',[xpos(xind),ypos(yind),fig.subplot.width,fig.subplot.height],...
                 'FontSize', 8);
    plot(pfd{1},unitsExampleSubset(y),'mean','text',mrate,false,'nanColor',nanColor,'colorMap',@jet);
    set(sp(end),'YTickLabels',{});
    set(sp(end),'XTickLabels',{});    
    sp(end).Position(1) = sp(end).Position(1)+0.1;    
    if y == 1,
% $$$         ylabel('Body Pitch');
% $$$         xlabel('Head Pitch');
        title({'BHV'});
    end
    
% PLOT space and theta phase jointly restricted bhv rate maps
    for p = 1:(pfs.adata.binSizes(2)/2),
        xind = 2+p;        
        sp(end+1) = axes('Units','centimeters',...
                         'Position',[xpos(xind),ypos(yind),fig.subplot.width,fig.subplot.height],...
                         'FontSize', 8);
        hold(sp(end),'on');
        imagescnan({pfd{1}.adata.bins{:},sq(sum(rmap(:,phzOrder([p*2-1,p*2]),:),2)./2)'},...
                   [0,mrate],[],false,nanColor,'colorMap',@jet);
        sp(end).Position(1) = sp(end).Position(1)+0.2;
        axis('xy');
        axis('tight');
        set(sp(end),'YTickLabels',{});
        set(sp(end),'XTickLabels',{});    
    end 
    text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
         pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
         sprintf('%2.0f',mrate),...
        'Color','w','FontWeight','bold','FontSize',8)
   
   
% COMPUTE circular statistics of rate changes within the theta cycle
    mask = double(sq(max(rmap,[],2))>4)';
    mask(mask==0) = nan;
    rmapComplex = sum(reshape(bsxfun(@times,                                                     ...
                                     reshape(permute(rmap,[1,3,2]),[],size(rmap,2)),             ...
                                     exp(-i.*pfs.adata.bins{2}(phzOrder))'),                     ...
                              [size(rmap,1),size(rmap,3),size(rmap,2)]),                         ...
                      3)                                                                         ...
        ./sq(sum(rmap,2));

    xind = 3+p;        
    sp(end+1) = axes('Units','centimeters',                                                      ...
                     'Position',[xpos(xind),                                                     ...
                        ypos(yind),                                                     ...
                        fig.subplot.width,                                              ...
                        fig.subplot.height],                                            ...
                     'FontSize', 8);
    hold(sp(end),'on');
    imagescnan({pfd{1}.adata.bins{:},angle(rmapComplex)'.*mask},...
               [-pi,pi],[],false,nanColor,'colorMap',@hsv);
    sp(end).Position(1) = sp(end).Position(1)+0.3;
    axis('xy');
    axis('tight');
    set(sp(end),'YTickLabels',{});
    set(sp(end),'XTickLabels',{});    
    

% RESULTANT length of theta phase dependent rate changes
   xind = 4+p;        
   sp(end+1) = axes('Units','centimeters',                                                      ...
                    'Position',[xpos(xind),                                                     ...
                                ypos(yind),                                                     ...
                                fig.subplot.width,                                              ...
                                fig.subplot.height],                                            ...
                    'FontSize', 8);
   hold(sp(end),'on');
   imagescnan({pfd{1}.adata.bins{:},abs(rmapComplex)'.*mask},...
              [0,0.5],[],false,nanColor,'colorMap',@copper);
   sp(end).Position(1) = sp(end).Position(1)+0.3;
   axis('xy');
   axis('tight');
   set(sp(end),'YTickLabels',{});
   set(sp(end),'XTickLabels',{});    
   

   

% GET unit row within population data
   uind = find(ismember(cluSessionMap,[t,unitsExampleSubset(y)],'rows'));   
   for e = 1:2
       if eigVars(uind,e)>15,
           xind = 4+p+e;        
% SETUP subplot
           sp(end+1) = axes('Units','centimeters',                                             ...
                            'Position',[xpos(xind),                                            ...
                               ypos(yind),                                                     ...
                               fig.subplot.width,                                              ...
                               fig.subplot.height],                                            ...
                            'FontSize', 8);
           hold(sp(end),'on');
% PLOT eigenvector
           imagescnan({pfd{1}.adata.bins{:},sq(eigVecs(uind,e,:,:))'},...
                      [],[],false,nanColor,'colorMap',@parula);
           sp(end).Position(1) = sp(end).Position(1)+0.4;
           axis('xy');
           axis('tight');
           set(sp(end),'YTickLabels',{});
           set(sp(end),'XTickLabels',{});    
           
           text(pfs.adata.bins{1}(end)-0.45*diff(pfs.adata.bins{1}([1,end])),...
                pfs.adata.bins{3}(end)-0.10*diff(pfs.adata.bins{3}([1,end])),...
                sprintf('%2.0f',eigVars(uind,e)),...
                'Color','w','FontWeight','bold','FontSize',8)
       else
           sp(end+1) = gobjects([1,1]);
       end
   end%for e
end%for y


% CREATE reference bars
axes(fax)
xlim([0,fax.Position(3)]);
ylim([0,fax.Position(4)]);

% SUPERAXIS Theta phase 
supAxisThpOffset = -(ny+1)*fig.subplot.height-(ny).*fig.subplot.verticalPadding+0.7;
% $$$ line([sp(3).Position(1),sum(sp(nx).Position([1,3]))],...
% $$$      sum(sp(3).Position([2,4])).*[1,1]+0.5,'LineWidth',1)
text(mean([sp(3).Position(1),sum(sp(nx-2).Position([1,3]))]),...
     sum(sp(3).Position([2,4]))-0.1+supAxisThpOffset,'Theta Phase','HorizontalAlignment','center','FontSize',8)
% TICKS theta phase
text(sp(3).Position(1),sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'0^{o}','HorizontalAlignment','center','FontSize',8)
text(mean([sp(3).Position(1),sum(sp(nx-2).Position([1,3]))]),...
     sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'180^{o}','HorizontalAlignment','center','FontSize',8)
text(sum(sp(nx-2).Position([1,3])),...
     sum(sp(3).Position([2,4]))+0.25+supAxisThpOffset,'360^{o}','HorizontalAlignment','center','FontSize',8)

% ADD phase colorbar to the top of the examples
hcb = colorbar(sp(3),'NorthOutside');
hcb.Units = 'centimeters';
hcb.Position(2) = sp(1).Position(2) + fig.subplot.height + 0.5 + supAxisThpOffset;
hcb.Position(3) = fig.subplot.width*round(numel(phzOrder)/2);
colormap(hcb,'hsv');
caxis([0,360]);
hcb.YTick = [];

% REFERENCE bar behavior space
line([sp(end-nx).Position(1),sum(sp(end-nx).Position([1,3]).*[1,11/28])],...
     sp(end-nx).Position(2).*[1,1]-0.1,'LineWidth',1);
text(sp(end-nx).Position(1),sp(end-nx).Position(2)-0.3,'1 rad','FontSize',8);

% REFERENCE bar physical space
line([sp(end-nx-1).Position(1),sum(sp(end-nx-1).Position([1,3]).*[1,0.5])],...
     sp(end-nx-1).Position(2).*[1,1]-0.1,'LineWidth',1);
text(sp(end-nx-1).Position(1),sp(end-nx-1).Position(2)-0.3,'50 cm','FontSize',8);


%fin








csthresh = 2;
csoffset = 1;
cluSessionMapSubset = cluSessionMap(unitSubsets{1},:);
cluSessionMapSubset_L = cluSessionMapSubset(fsrcz(:,1)>2.1&fsrcz(:,2)<(csthresh-csoffset)&fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_H = cluSessionMapSubset(fsrcz(:,3)>2.1&fsrcz(:,1)<(csthresh-csoffset)&fsrcz(:,2)<(csthresh-csoffset),:);
cluSessionMapSubset_R = cluSessionMapSubset(fsrcz(:,2)>2.1&fsrcz(:,1)<(csthresh-csoffset)&fsrcz(:,3)<(csthresh-csoffset),:);
cluSessionMapSubset_N = cluSessionMapSubset(fsrcz(:,2)<csthresh&fsrcz(:,1)<csthresh&fsrcz(:,3)<csthresh,:);
cluSessionMapSubset_C = cluSessionMapSubset(~ismember(cluSessionMapSubset,...
                                                  [cluSessionMapSubset_L;...
                                                   cluSessionMapSubset_H;...
                                                   cluSessionMapSubset_R;...
                                                   cluSessionMapSubset_N],'rows'),:);






%csmUnitSubset = cluSessionMap(unitSubsets{1},:);
csmEigScrs = eigScrs;
csmEigVecs = eigVecs;
csmEigVars = eigVars;

clear('i');
[~,csmPhzPrefMax] = max(sq(csmEigScrs(:,:,:)),[],3);


%csmPhzPrefMean = angle(sum(reshape(bsxfun(@times,reshape(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),[],size(csmEigScrs,3)),exp(-i.*pfs.adata.bins{2})'),size(csmEigScrs)),3)./sum(bsxfun(@plus,csmEigScrs,abs(min(csmEigScrs,[],3))),3));





figure,plot(pfs.adata.bins{2}(abs(csmPhzPrefMax(:,1)))+randn(size(csmPhzPrefMax(:,1)))/10,...
            csmPhzPrefMean(:,1),'.');
% $$$ figure,plot(sq(csmEigScrs(:,1,:))')

figure,plot(log10(csmEigVars(:,1)),log10(csmEigVars(:,2)),'.')

figure();
hold('on');
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,1),'.b')
plot(circ_dist(csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)),csmEigVars(:,2),'.r')


figure,
plot(circ_mean([csmPhzPrefMean(:,1),csmPhzPrefMean(:,2)]')+pi,csmEigVars(:,1),'.b')



p = pfs.adata.bins{2};
for h = 1:size(rmap,1),
    for b = 1:size(rmap,3),
        rc = rmap(h,:,b)'.*exp(i*p);
        Rmat = rc*conj(rc)';
        R(h,b) = sum(sum(triu(Rmat,1)))./sum(sum(triu(rmap(h,:,b)'*rmap(h,:,b),1)));
    end
end

figure();
subplot(121);
imagesc(abs(R)'),axis('xy');
subplot(122);
imagescnan({pfd{1}.adata.bins{:},abs(rmapComplex)'},...
           [],[],false,nanColor,'colorMap',@parula);
axis('xy')




figure(), 
pBins = linspace(0,2*pi,17);
ks = [0.1:0.2:2];
sp = tight_subplot(2,5,0.05,0.05);
for k = 1:numel(ks),
    axes(sp(k));
    hold('on');
    for b = [10,100,1000],
        for c = 1:1000,
            th = VonMisesRnd(0,ks(k),b,1);
            pBinsInds = discretize(th,pBins);    
            thCount = accumarray(pBinsInds,ones(size(th)),[numel(pBins)-1,1],@sum);
            rc = thCount.*exp(i*p);
            Rmat = rc*conj(rc)';
            R = sum(sum(triu(Rmat,1)))./sum(sum(triu(thCount*thCount',1)));
            Rasm(c) = abs(R);
            Rppc(c) = PPC(th);
        end
        plot(Rppc,Rasm,'.');
    end
end

af(@(s) xlim(s,[-0.2,1]),sp);
af(@(s) ylim(s,[-0.2,1]),sp);
af(@(s) grid(s,'on'),sp);