% new place field calculations

%Trial = MTATrial.validate('jg05-20120312.cof.all');
%Trial = Trials{20};

sampleRate = 32;
state = 'theta-groom-sit';
stcMode = 'msnn_ppsvd_raux';

latticeInterval = 20;
xpos = -500:latticeInterval:500;
ypos = -500:latticeInterval:500;
latticeSize = [numel(xpos),numel(ypos)];
xbins = xpos;
ybins = ypos;

sigma = 40;
sigmaD = 5*sigma;
sigmaDS = 2*sigma^2;

rmap = {};
rmapShuff = {};


for tind = 1:numel(Trials)

    Trial = Trials{tind};
    Trial.load('stc',stcMode);

    tper = [Trial.stc{[state]}];
    tper.resample(sampleRate);

    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    hba = fet_hba(Trial,'xyz',xyz);
    thba = sq(hba(tper));

    rot = transform_vector_to_rotation_matrix( ...
            xyz,{'hcom','nose'}, Trial.meta.correction.headYaw);

    units = unitsEgo{tind};

    spk = Trial.load('spk',sampleRate,state,units);

    rmap{tind} = nan([numel(xpos),numel(ypos),numel(units),numel(hbaBinEdg)-1]);
    %% cell array version
    for u = 1:numel(units)
        tic
        [mxr,mxp] = pft{tind}.maxRate(units(u));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                      bsxfun(                                         ...
                                          @plus,                                       ...
                                          multiprod(bsxfun(@minus,                     ...
                                                      mxp,                             ...
                                                      sq(xyz(:,'hcom',[1,2]))),        ...
                                                    rot,2,[2,3]),                        ...
                                          [0,0]),                       ...   
                                      sampleRate,                                      ...
                                      'egocentric_placefield',                         ...
                                      'egopfs',                                        ...
                                      'p'                                              ...
                                      );
        tpos = pfsCenterHR(tper,:);    
        mapSpk = {};
        mapOcc = {};
        mapInd = {};
        mapW   = {};
        mapHba = {};
        
        % Collect the occupancy 
        for xind = 1:latticeSize(1)
            for yind = 1:latticeSize(2)
                tempPos = bsxfun(@minus,tpos,[xpos(xind),ypos(yind)]);
                tempDst = sqrt(sum(tempPos.^2,2));
                pind = tempDst < sigmaD;
                mapOcc{xind,yind} = tempPos(pind,:);
                mapInd{xind,yind} = find(pind);
                mapHba{xind,yind} = discretize(thba(pind),hbaBinEdg);
            end
        end
        
        occ = zeros([numel(mapOcc),3]);
        for lind = 1:numel(mapOcc)
            for hbaInd = 1:numel(hbaBinEdg)-1
                occ(lind,hbaInd) = sum(exp(-sum(mapOcc{lind}(mapHba{lind}==hbaInd,:).^2,2)./(sigmaDS)));
            end
        end
        
        mapSpk = {};
        mapW = {};
        tres = SelectPeriods(spk(units(u)),tper.data,'d',1,1);
        [ures,IA,IC] = unique(tres);
        % GET spike degeneracy weights
        dic = [diff([IC(1)-1;IC])];    
        
        % map degenerate spikes
        wpos = zeros([size(tpos,1),1]);
        k = 0;
        for di = dic'
            if di==0
                wpos(ures(k)) = wpos(ures(k))+1;
            else
                k = k+1;                        
                if k == numel(ures)+1
                    break;
                end
                wpos(ures(k)) = 1;
            end
        end
        %cwpos = wpos(wpos~=0);
        spos = nan([size(tpos,1),2]);
        spos(tres,:) = tpos(tres,:);
        % assume tres is monotonically increasing
        for hbaInd = 1:numel(hbaBinEdg)-1
            for xind = 1:latticeSize(1)
                for yind = 1:latticeSize(2)
                    tempPos = bsxfun(@minus,spos,[xpos(xind),ypos(yind)]);
                    mapSpk{xind,yind} = tempPos(mapInd{xind,yind}(mapHba{xind,yind}==hbaInd),:);
                    mapW{xind,yind} = wpos(mapInd{xind,yind}(mapHba{xind,yind}==hbaInd));
                end
            end
            scc = zeros([numel(mapOcc),1]);
            for lind = 1:numel(mapOcc)
                scc(lind,hbaInd) = sum(mapW{lind}.*exp(-sum(mapSpk{lind}.^2,2)./(sigmaDS)),'omitnan');
            end
            rmap{tind}(:,:,u,hbaInd) = reshape(scc(:,hbaInd)./(occ(:,hbaInd)./sampleRate),latticeSize);
        end
        toc
    end



    %%%<<< permuted hba

    rmapShuff{tind} = nan([numel(xpos),numel(ypos),numel(units),numel(hbaBinEdg)-1,100]);
    %% cell array version
    for u = 1:numel(units)
        tic

        [mxr,mxp] = pft{tind}.maxRate(units(u));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                          bsxfun(                                         ...
                                              @plus,                                       ...
                                              multiprod(bsxfun(@minus,                     ...
                                                          mxp,                             ...
                                                          sq(xyz(:,'hcom',[1,2]))),        ...
                                                        rot,2,[2,3]),                        ...
                                              [0,0]),                       ...   
                                          sampleRate,                                      ...
                                          'egocentric_placefield',                         ...
                                          'egopfs',                                        ...
                                          'p'                                              ...
                                          );
        tpos = pfsCenterHR(tper,:);    
        mapSpk = {};
        mapOcc = {};
        mapDst = {};
        mapInd = {};
        mapW   = {};

        mapHba = {};

        % Collect the occupancy 
        for xind = 1:latticeSize(1)
            for yind = 1:latticeSize(2)
                %lpdist = sqrt(sum(bsxfun(@minus,tpos,[xpos(xind),ypos(yind)]).^2,2)) < 2*sigma;
                tempPos = bsxfun(@minus,tpos,[xpos(xind),ypos(yind)]);
                tempDst = sqrt(sum(tempPos.^2,2));
% $$$         mapDst(sub2ind(latticeSize,xind,yind)) = tempDst(tempDst < sigmaD);
                pind = tempDst < sigmaD;
                mapOcc{xind,yind} = tempPos(pind,:);
                mapInd{xind,yind} = find(pind);
                mapHba{xind,yind} = discretize(thba(pind),hbaBinEdg);
            end
        end

        for iter = 1:100
            tMHba = mapHba;
            occ = zeros([numel(mapOcc),3]);
            for lind = 1:numel(mapOcc)
                tMHba{lind} = mapHba{lind}(randperm(numel(mapHba{lind})));
                for hbaInd = 1:numel(hbaBinEdg)-1
                    occ(lind,hbaInd) = sum(exp(-sum(mapOcc{lind}(tMHba{lind}==hbaInd,:).^2,2)./(sigmaDS)));
                end
            end



            mapSpk = {};
            mapW = {};
            tres = SelectPeriods(spk(units(u)),tper.data,'d',1,1);
            [ures,IA,IC] = unique(tres);
            % GET spike degeneracy weights
            dic = [diff([IC(1)-1;IC])];    
            
            % map degenerate spikes
            wpos = zeros([size(tpos,1),1]);
            k = 0;
            for di = dic'
                if di==0
                    wpos(ures(k)) = wpos(ures(k))+1;
                else
                    k = k+1;                        
                    if k == numel(ures)+1
                        break;
                    end
                    wpos(ures(k)) = 1;
                end
            end
            %cwpos = wpos(wpos~=0);
            spos = nan([size(tpos,1),2]);
            spos(tres,:) = tpos(tres,:);
            % assume tres is monotonically increasing
            for hbaInd = 1:numel(hbaBinEdg)-1
                for xind = 1:latticeSize(1)
                    for yind = 1:latticeSize(2)
                        tempPos = bsxfun(@minus,spos,[xpos(xind),ypos(yind)]);
                        mapSpk{xind,yind} = tempPos(mapInd{xind,yind}(tMHba{xind,yind}==hbaInd),:);
                        mapW{xind,yind} = wpos(mapInd{xind,yind}(tMHba{xind,yind}==hbaInd));
                    end
                end
                scc = zeros([numel(mapOcc),1]);
                for lind = 1:numel(mapOcc)
                    scc(lind,hbaInd) = sum(mapW{lind}.*exp(-sum(mapSpk{lind}.^2,2)./(sigmaDS)),'omitnan');
                end
                rmapShuff{tind}(:,:,u,hbaInd,iter) = reshape(scc(:,hbaInd)./(occ(:,hbaInd)./sampleRate),latticeSize);
            end
        end
    end
end

%%%>>>
% END permuted hba

% $$$ %mask = create_tensor_mask({xpos,ypos})
% $$$ mask = double(sqrt(bsxfun(@plus,xbins.^2,ybins'.^2)') < 445);
% $$$ mask(~mask) = nan;
% $$$ 
% $$$ u = find(units==21);
% $$$ iter = 15;
% $$$ %u = find(units==17);
% $$$ figure,
% $$$ for hbaInd = 1:numel(hbaBinEdg)-1
% $$$ subplot2(2,numel(hbaBinEdg)-1,1,hbaInd);
% $$$ shading(gca(),'flat');
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,fliplr(rot90(rmap(:,:,u,hbaInd)',-1)).*mask),'EdgeColor','none');
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ Lines([],0,'k');
% $$$ Lines(0,[],'k');
% $$$ subplot2(2,numel(hbaBinEdg)-1,2,hbaInd);
% $$$ shading(gca(),'flat');
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2, ...
% $$$            fliplr(rot90(rmapShuff(:,:,u,hbaInd,iter)',-1)).*mask),'EdgeColor','none');
% $$$ caxis([2,12])
% $$$ axis('xy');
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ Lines([],0,'k');
% $$$ Lines(0,[],'k');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ nIter= 1000;
% $$$ occBS = nan([numel(mapOcc),nIter]);
% $$$ sccBS = nan([numel(mapOcc),nIter]);
% $$$ tic
% $$$ for lind = 1:numel(mapOcc)
% $$$     if numel(mapSpk{lind})>256
% $$$         for iter = 1:nIter
% $$$             [spkSubsample,rid] = datasample(mapSpk{lind},128,1);
% $$$             sccBS(lind,iter)  = sum(mapW{lind}(rid).*exp(-sum(spkSubsample.^2,2)./(sigmaDS)),'omitnan');
% $$$             occBS(lind,iter) = sum(exp(-sum(mapOcc{lind}(rid,:).^2,2)./(sigmaDS)));
% $$$         end
% $$$     end
% $$$ end
% $$$ toc
% $$$ 
% $$$ 
% $$$ rmapBS = reshape(mean(sccBS./(occBS./sampleRate),2,'omitnan'),latticeSize);
% $$$ smapBS =  reshape(std(sccBS./(occBS./sampleRate),[],2,'omitnan'),latticeSize);
% $$$ 
% $$$ figure();
% $$$ subplot(131);
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,rmapBS'.*mask),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ subplot(132);
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,smapBS'.*mask),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ subplot(133);
% $$$ %set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,rmapBS'.^2./smapBS'.^2.*mask),'EdgeColor','none');
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,smapBS'./rmapBS'.*mask),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ 
% $$$ %% Block shuffling
% $$$ 
% $$$ nIter= 1000;
% $$$ occb = nan([numel(mapOcc),nIter]);
% $$$ sccb = nan([numel(mapOcc),nIter]);
% $$$ blockSize = 8;
% $$$ for iter = 1:nIter
% $$$     tic
% $$$ mapSpkBlock = cell([numel(mapSpk),1]);
% $$$ mapOccBlock = cell([numel(mapOcc),1]);
% $$$ mapWBlock = cell([numel(mapW),1]);
% $$$ 
% $$$ for lind = 1:numel(mapOcc)
% $$$     if numel(mapSpk{lind})>256
% $$$         randshift = randi([1,round(size(mapSpk{lind},1)/2)],1)-round(size(mapSpk{lind},1)/4);
% $$$         tms = reshape(circshift(mapSpk{lind}(1:end-mod(size(mapSpk{lind},1),blockSize),:),randshift),[],8,2);
% $$$         tmo = reshape(circshift(mapOcc{lind}(1:end-mod(size(mapOcc{lind},1),blockSize),:),randshift),[],8,2);
% $$$         tmw = reshape(circshift(  mapW{lind}(1:end-mod(size(mapW{lind}  ,1),blockSize),:),randshift),[],8,1);
% $$$         [mapSpkBlock{lind},rid] = datasample(tms,16,1,'Replace',true);
% $$$         mapOccBlock{lind} = tmo(rid,:,:);
% $$$         mapWBlock{lind} = tmw(rid,:,:);
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ mask = double(sqrt(bsxfun(@plus,xbins.^2,ybins'.^2)') < 445);
% $$$ blockIndex = ~cellfun(@isempty,mapSpkBlock)&mask(:);
% $$$ 
% $$$ mapSpkBlockMat = reshape(permute(cat(4,mapSpkBlock{blockIndex}),[4,1,2,3]),[],8,2);
% $$$ rperm = randperm(size(mapSpkBlockMat,1));
% $$$ mapSpkBlockMat = reshape(mapSpkBlockMat(rperm,:,:),[],16*8,2);
% $$$ mapWBlockMat = reshape(permute(cat(4,mapWBlock{blockIndex}),[4,1,2,3]),[],8,1);
% $$$ mapWBlockMat = reshape(mapWBlockMat(rperm,:,:),[],16*8,1);
% $$$ 
% $$$ mapOccBlockMat = reshape(permute(cat(4,mapOccBlock{blockIndex}),[4,1,2,3]),[],128,2);
% $$$ 
% $$$ 
% $$$ 
% $$$ %sccb = nan(latticeSize);
% $$$ %sccb = nan([prod(latticeSize),1]);
% $$$ %occb = nan(latticeSize);
% $$$ %occb = nan([prod(latticeSize),1]);
% $$$ sccb(blockIndex,iter) = sum(mapWBlockMat.*exp(-sum(mapSpkBlockMat.^2,3)./(sigmaDS)),2,'omitnan');
% $$$ occb(blockIndex,iter) = sum(exp(-sum(mapOccBlockMat.^2,3)./(sigmaDS)),2);
% $$$ 
% $$$ 
% $$$ toc
% $$$ end
% $$$ 
% $$$ rmapb = sccb./(occb./16);
% $$$ 
% $$$ 
% $$$ mask = double(sqrt(bsxfun(@plus,xbins.^2,ybins'.^2)') < 445);
% $$$ mask(~mask) = nan;
% $$$ rmapb = sccb./(occb./16);
% $$$ sax = gobjects([1,0]);
% $$$ figure();
% $$$ sax(end+1) = subplot(151);
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,rmap'.*mask),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ sax(end+1) = subplot(152);
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,reshape(mean(rmapb,2),latticeSize)'.*mask),'EdgeColor','none');
% $$$ title('mean')
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ sax(end+1) = subplot(153);
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,reshape(std(rmapb,[],2),latticeSize)'.*mask),'EdgeColor','none');
% $$$ title('std')
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ sax(end+1) = subplot(154);
% $$$ rstd = reshape(std(rmapb,[],2),latticeSize)';
% $$$ rmean = reshape(mean(rmapb,2),latticeSize)';
% $$$ zmap = (rmap'-rmean)./rstd.*mask;
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,zmap),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ hold(gca(),'on');
% $$$ contour(xpos,ypos,(rmap'-rmean)./rstd.*mask,[3,3],'-m','LineWidth',2);
% $$$ copyobj(sax(end).Children(1),sax(1));
% $$$ sax(end+1) = subplot(155);
% $$$ rzmap = zmap;
% $$$ rzmap(rzmap<1) = 1;
% $$$ rzmap(nniz(rzmap)) = 1;
% $$$ set(pcolor(xpos-diff(xpos(1:2))/2,ypos-diff(ypos(1:2))/2,(rmap'.*mask).*(log10(rzmap)./max(log10(rzmap(:))))),'EdgeColor','none');
% $$$ ylim([ypos([1,end])+[-1,1].*diff(ypos(1:2))/2])
% $$$ xlim([xpos([1,end])+[-1,1].*diff(xpos(1:2))/2])
% $$$ colormap('jet');
% $$$ colorbar();
% $$$ axis('xy');
% $$$ 
% $$$ 
% $$$ rcenter = [160,-180];
% $$$ zcenter = [250,-210];
% $$$ 
% $$$ sampleRate = 250;
% $$$ 
% $$$ spk = Trial.load('spk',sampleRate,state,units);
% $$$ 
% $$$ xyz = preproc_xyz(Trial,'trb');
% $$$ xyz.resample(sampleRate);
% $$$ 
% $$$ %ghz = compute_ghz(Trial,units);
% $$$ % TODO Compute head referenced field center coordinates for 2d field and plot spikes
% $$$ 
% $$$ 
% $$$ 
% $$$ 
