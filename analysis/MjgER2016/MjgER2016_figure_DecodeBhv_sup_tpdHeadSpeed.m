
% COMPUTE ------------------------------------------------------------------------------

%%%<<< COMPUTE error for TDP decoding per head speed bin

% COMPUTE head speed
hvxyBinEdges = [-2,-1,0.2,0.75,1.2,1.6,2];
hvxy = cf(@(x)  vel(filter(copy(x),'ButFilter',3,2.5,'low'),'hcom',[1,2]),  dc4.xyz);
vtid = {};
for t = 1:numTrials,
    hvxy{t}.data(hvxy{t}.data<1e-3) = 1e-3;
    hvxy{t}.data(stcm{t}(:,1)~=1) = 1e-3;    
    hvxy{t} = log10(hvxy{t}.data);
    vtid{t} = t.*ones([size(hvxy{t},1),1]);
end
vtid =  cat(1,vtid{:});
chvxy = cat(1,hvxy{:});
cstcm = cat(1,stcm{:});
uchvxy = chvxy;
% Convert to uniform distribution
uind = -2 < chvxy ...
       & chvxy < 2 ...
       & logical(cstcm(:,1)) ...
       & ~any(logical(cstcm(:,[2,7,8])),2);
[uchvxy(uind),xx] = MakeUniformDistr(chvxy(uind));
uhvxy = {};
for t = 1:numTrials,
    uhvxy{t} = uchvxy(vtid==t);
end

figure();
hold('on');
for sts = 3:6,
    subplot(4,1,sts-2);
    uind = -2 < chvxy ...
           & chvxy < 2 ...
           & logical(cstcm(:,1)) ...
           & logical(cstcm(:,sts)) ...       
           & ~any(logical(cstcm(:,[2,7,8])),2);
    histogram(chvxy(uind),vxyBinEdges);
end



% GET the velocity bin edges of the uniform distribution
vxyBinEdges = linspace(-2,2,13);
vxyBinOri = [];
for vi = 1:numel(vxyBinCenters),
   ss = chvxy(cstcm(:,1)==1 & WithinRanges(uchvxy,vxyBinEdges([vi,vi+1])));
   vxyBinOri(vi,:) = [min(ss),mean(ss),max(ss)];
end

vxyBinCenters = mean(GetSegs(vxyBinEdges,1:numel(vxyBinEdges)-1,2));
dcTDPfrontalV = repmat({cell([1,numTrials])},[1,numel(vxyBinCenters)]);
dcTDPlateralV = repmat({cell([1,numTrials])},[1,numel(vxyBinCenters)]);
indTDPV = cell([1,numel(stid)]);
errorBinEdges = linspace(-500,500,100);      
errorBinCenters = mean(GetSegs(errorBinEdges,1:numel(errorBinEdges)-1,2));
for vi = 1:numel(vxyBinCenters),
    indTDPV{vi}  = cf(@(p,u,s,v)                                            ...
                      all(p>0.0005,2)                                       ...
                      & sum(double(u>=2),2)>6                               ...
                      & s(:,1)==1                                           ...
                      & any(logical(s(:,[5,6])),2)                         ...                      
                      & ~any(logical(s(:,[2,7,8])),2)                         ...
                      & WithinRanges(v,vxyBinEdges([vi,vi+1])),             ...
                      posteriorMaxTPD,unitInclusionTPD,stcm,uhvxy);
    for t = 1:numel(Trials),
        for j = 1:numel(phzBins)-1;
            dcTDPfrontalV{vi}{t} =                                             ...
                cat(2,                                                          ...
                    dcTDPfrontalV{vi}{t},                                      ...
                    histcounts(sq(dErr{t}(indTDPV{vi}{t},1,j)),errorBinEdges)');
            dcTDPlateralV{vi}{t} =                                             ...
                cat(2,                                                          ...
                    dcTDPlateralV{vi}{t},                                      ...
                    histcounts(sq(dErr{t}(indTDPV{vi}{t},2,j)),errorBinEdges)');
        end
    end
end
%%%>>>

%%%<<< COMPUTATION : theta phase ridge values for each state's dcTDPfrontal error
dcTpdFrntVSum = {};
dcTpdVInd = [];
dcTpdVMax = [];
dcTpdVWmd = [];
dcTpdVWmn = [];
dcTpdVVar = [];
dcTpdVEnt = [];
ebcSub = errorBinCenters;
for sts = 1:numel(vxyBinCenters),
    dcTpdFrntVSum{sts} = RectFilter(sum(cat(3,dcTDPfrontalV{sts}{:}),3),5,5);
    dcTpdFrntVSum{sts} = bsxfun(@rdivide,dcTpdFrntVSum{sts},sum(dcTpdFrntVSum{sts}));
    [~, dcTpdVInd(end+1,:)] = max(dcTpdFrntVSum{sts});
    dcTpdVMax(end+1,:) = ebcSub(dcTpdVInd(end,:)');
    dcTpdVWmd(end+1,:) = errorBinCenters(sum(double(cumsum(dcTpdFrntVSum{sts}) < 0.5)));        
    dcTpdVWmn(end+1,:) = sum(bsxfun(@times,dcTpdFrntVSum{sts},ebcSub'))';
    dcTpdVVar(end+1,:) = 1./mean(sqrt(reshape(bsxfun(@minus,...
                                                     permute(dcTpdFrntVSum{sts},[1,3,2]),...
                                                     permute(dcTpdFrntVSum{sts},[3,1,2])).^2,...
                                              [],size(dcTpdFrntVSum{sts},2))));
    dcTpdVEnt(end+1,:) = -sum(dcTpdFrntVSum{sts}.*log2(dcTpdFrntVSum{sts}));
end
%%%>>>




% PLOT ---------------------------------------------------------------------------------

%%%<<< PLOT frontal error for TDP decoding
% ADJUST subplot coordinates

figure
xcoords = linspace(-500,500,250);
ycoords = circ_rad2ang([phzBinCenters;phzBinCenters+2*pi]);
for sts = 1:numel(vxyBinCenters),

    dcTDPfrontalVSum{sts} = RectFilter(sum(cat(3,dcTDPfrontalV{sts}{:}),3),5,1);
    subplot(numel(vxyBinCenters),1,sts);
    hold('on');
    imagesc(xcoords,                               ...
            ycoords,                               ...
            repmat(dcTDPfrontalVSum{sts},[1,2])');
    axis('xy');
    xlim([-200,200]);
    ylim([-60,420]);            
    % MASK areas outside single cycle
    patch([-300,-300,300,300],...
          [ 360, 420,420,360],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    patch([-300,-300,300,300],...
          [ -60,   0,  0,-60],...
          [0.2,0.2,0.2],...
          'EdgeAlpha',0,...
          'FaceAlpha',0.4);
    sax(end).XTick = [-150,0,150];
    sax(end).XTickLabels = [-15,0,15];
    Lines(0,[],'k');

    %plot(repmat(dcTpdVWmd(sts,:),[1,2])',ycoords,'g','LineWidth',1);    
    plot(repmat(dcTpdVWmn(sts,:),[1,2])',ycoords,'k','LineWidth',1);    
    %plot(repmat(dcTpdVMax(sts,:),[1,2])',ycoords,'k','LineWidth',1);
end
colormap('jet');


vcmap = jet(numel(vxyBinCenters));
figure();
hold('on');    
for vi = 1:numel(vxyBinCenters),
    plot(repmat(dcTpdVWmn(vi,:),[1,2])',...    
         [phzBinCenters;phzBinCenters+2.*pi],'Color',vcmap(vi,:));
end



vcmap = jet(numel(vxyBinCenters));
figure();
hold('on');    
for vi = 1:numel(vxyBinCenters),
    plot(repmat(dcTpdVVar(vi,:),[1,2])',...    
         [phzBinCenters;phzBinCenters+2.*pi],'Color',vcmap(vi,:));
end


vcmap = jet(numel(vxyBinCenters));
figure();
hold('on');    
for vi = 1:numel(vxyBinCenters),
    plot(repmat(dcTpdVEnt(vi,:),[1,2])',...    
         ycoords,'Color',vcmap(vi,:));
end
xl = [5.25,6.5];
xlim(xl);
ylim([-60,420]);
% MASK areas outside single cycle
patch([xl(1),xl(1),xl(2),xl(2)],...
      [ 360, 420,420,360],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
patch([xl(1),xl(1),xl(2),xl(2)],...
      [ -60,   0,  0,-60],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
Lines([],180,'k');

[vb,pb] = ndgrid(10.^[vxyBinOri(:,2)],circ_rad2ang([phzBinCenters;phzBinCenters+2.*pi]));
figure
pcolor(vb, pb, repmat(dcTpdVWmn,[1,2]))
colormap(gca(),'jet');    
xl = vxyBinOri([1,end],2);
xlim(xl);
ylim([-60,420]);            
% MASK areas outside single cycle
patch([xl(1),xl(1),xl(2),xl(2)],...
      [ 360, 420,420,360],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);
patch([xl(1),xl(1),xl(2),xl(2)],...
      [ -60,   0,  0,-60],...
      [0.2,0.2,0.2],...
      'EdgeAlpha',0,...
      'FaceAlpha',0.4);


figure,
    imagesc(xcoords,                               ...
            ycoords,                               ...
            repmat(dcTDPfrontalVSum{sts},[1,2])');
Lines(errorBinCenters(find(cumsum(dcTpdFrntVSum{sts}(:,2))>0.5,1,'first')),[],'k')
dcTpdVWmd =[];
axis('xy')
    dcTpdFrntVSum{sts} = RectFilter(sum(cat(3,dcTDPfrontalV{sts}{:}),3),5,1);
    dcTpdFrntVSum{sts} = bsxfun(@rdivide,dcTpdFrntVSum{sts},sum(dcTpdFrntVSum{sts}));
    dcTpdVWmd(end+1,:) = errorBinCenters(sum(double(cumsum(dcTpdFrntVSum{sts}) < 0.5)));    
    dcTpdVWmn(end+1,:) = sum(bsxfun(@times,dcTpdFrntVSum{sts},errorBinCenters'))';
    
    
hold('on');
plot(repmat(sum(bsxfun(@times,dcTpdFrntVSum{sts},errorBinCenters'))',[2,1]),...
     ycoords,'k')
