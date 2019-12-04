;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.


pxz = [];
shifts = 0:8:2^8;
for sts = 1:6,
    ind =   logical(dstcm(:,1))                             ...
            & logical(dstcm(:,sts))                         ...
            & ~any(logical(dstcm(:,[7,8])),2)               ...
            & dpostI                                        ...
            & duincI;        
    for t = [1:numel(ttid),numel(ttid)+1],
        if t <= numel(ttid)
            tid = tind(ttid(t));
        end
        for i = 1:2
            if i == 1,
                if t == numel(ttid)+1,                
                    indb = ind & smMask & duinc < 7;
                else
                    indb = ind & smMask & duinc < 7 & ( dtind == tid );
                end
            else
                if t == numel(ttid)+1,                                
                    indb = ind & smMask & duinc >= 7;
                else
                    indb = ind & smMask & duinc >= 7 & ( dtind == tid );
                end
            end
                
            for f = 1:4
                pxz(:,:,i,sts,t,f) = ...
                    histcounts2(ferr{f}(indb),                  ...
                                dphz(indb),                     ...
                                ferrorBinEdges{f},              ...
                                phzBins,                        ...
                                'Normalization','probability');
            end
        end
    end
end


figure();
for sts = 1:6,
    subplot2(6,2,sts,1);
        imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,1),[1,2])');
        axis('xy')
    subplot2(6,2,sts,2);
        imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,1),[1,2])');
        axis('xy')
end
colormap('jet');
ForAllSubplots('Lines(0,[],''k'');');
ForAllSubplots('Lines([],2*pi,''k'');');
ForAllSubplots('Lines([],pi,''k'');');
ForAllSubplots('Lines([],0,''k'');');
ForAllSubplots('caxis([0,0.01]);');
ForAllSubplots('caxis([0,0.01]);');
ForAllSubplots('xlim([-300,300]);');

figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,2),[1,2])');
axis('xy')
Lines(0,[],'k');
subplot2(6,2,sts,2);
imagesc(errorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,2),[1,2])');
axis('xy')
Lines(0,[],'k');
end
colormap('jet');


figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,4),[1,2])');
axis('xy')
subplot2(6,2,sts,2);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,4),[1,2])');
axis('xy')
end

figure,
for sts = 1:6,
subplot2(6,2,sts,1);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,1,sts,8,3),[1,2])');
axis('xy')
subplot2(6,2,sts,2);
imagesc(perrorBinCenters,[phzBinCenters;phzBinCenters+2.*pi],repmat(pxz(:,:,2,sts,8,3),[1,2])');
axis('xy')
end


figure,
for sts = 1:6,
    subplot(6,1,sts);
    ind =   logical(dstcm(:,1))                             ...
            & logical(dstcm(:,sts))                         ...
            & ~any(logical(dstcm(:,[7,8])),2)               ...
            & dpostI                                        ...
            & duincI;        
    out = hist2([duinc(ind),dphz(ind)],1:16,phzBins);
    imagesc(1:16,[phzBinCenters;phzBinCenters+2*pi],repmat(bsxfun(@rdivide,out,sum(out)),[1,2])');
    axis('xy');
end
