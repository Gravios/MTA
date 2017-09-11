function xyz = repair_rigidbody_partial_reconstruction(Session,varargin)
%function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model,gind)
% 
% Tries to correct marker swaps and markers which wander from
% a ridgid body constrains
%
% Caution: This function is quick and dirty ...
%
% Note: This function is interactive and depends on the function ClusterPP.
% 
% variables:
%
%    xyz: MTADxyz, motion capture "mocap" data 
%
%  varargin:
%
%    goodIndex: double, the index which servers as a "good" template.
%                        if goodIndex is empty a gui will help you out :)
%
%    rb_model: MTAModel, model which contains the markers with a
%                        ridgid body constraint.
%
%    MarkerSwapErrorThreshold: double, Threshold for switching to auxilary
%                                      Error correction methods.
%
%    bufferSize: double, For machines with low memory (e.g. <4gb)
%
%  Note: the number of markers in rb_model must be <= 6.
%  
%  TODO: Make the assignment of MarkerSwapErrorThreshold automatic
%        Make memory use automatic (will probably occur upon the rapture :() 
%
%      case 'BEST_SWAP_PERMUTATION'
%      case 'RIGIDBODY_PARTIAL_RECONSTRUCTION'        


% TESTING vars
% Session = MTASession('jg04-20120213');
% varargin = [];


% DEFARGS ------------------------------------------------------------------------------------------
defArgs = struct('goodIndicies', [],                                                             ...
                 'rigidBodyMarkers',  {{'head_back','head_left','head_front','head_right'}},     ...
                 'errorThreshold',0.01,                                                          ...
                 'bufferSize',              2^16                                                 ...
);

[goodIndicies,rigidBodyMarkers,MarkerSwapErrorThreshold,bufferSize] = DefaultArgs(varargin,defArgs,'--struct');

% END DEFARGS ---------------------------------------------------------------------------------------



% DEFVARS ----------------------------------------------------------------------
if ischar(Session)||isa(Session,'MTATrial'),
    Session = MTASession(Session);
end
xyz = Session.load('xyz');

numChunks = floor(xyz.size(1)./bufferSize);
lastPiece = (numChunks*bufferSize+1):xyz.size(1);

% END DEFVARS ----------------------------------------------------------------------



% MAIN -------------------------------------------------------------------------    


% CHECK and SANITIZE rigidBodyMarkers input type
if iscell(rigidBodyMarkers)&&~isa(rigidBodyMarkers,'MTAModel'),
    rigidBody = xyz.model.rb(rigidBodyMarkers);
end
assert(rigidBody.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);
assert(rigidBody.N>3,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be gt 3']);


% CREATE a rigidbody model and xyz object
rxyz = xyz.copy;
trxyz = rxyz.copy;
rxyz.data = xyz(:,rigidBody.ml,:);
rxyz.model = rigidBody;


%
rxyz = xyz.copy;
rxyz.model = rigidBody;
rxyz.data = xyz(:,rigidBodyMarkers,:);

% CREATE an orthogonal coordinate system base on each Triad
markerTriad = compute_marker_triads(rxyz);


% SELECT good indicies if none are specified: select with the aid of a gui
if isempty(goodIndicies),
    %pfig = PlotSessionErrors(Session);
    pfig = figure(gen_figure_id);
    plot([reshape(markerTriad.imd,size(markerTriad.imd,1),[]),markerTriad.imo.*10]);
    goodIndicies = select_ranges(pfig);    
    fprintf(['\nNumber of indicies selected: = ',num2str(numel(goodIndicies)),'\n\n'])
    delete(pfig);
end




% COMPUTE the orthonormal triad basis for each frame
[~,~,mtBasis] = cellfun(@(x) bsxfun(@rdivide,x,sqrt(sum(x.^2,2))),...
                        cellfun(@squeeze,mat2cell(markerTriad.coor,...
                                                  ones([size(markerTriad.coor,1),1]),...
                                                  ones([size(markerTriad.coor,2),1]),3,3),...
                                'UniformOutput',false),...
                        'UniformOutput',false);

reconstructedSolutions  = nan([size(markerTriad.nck,1),rigidBody.N-3,xyz.size(3),numel(goodIndicies)]);
reconstructedSolutionsMarkerInds  = nan([size(markerTriad.nck,1),rigidBody.N-3]);
gcount = 1;
for goodIndex = goodIndicies(:)'
    for nck = 1:size(markerTriad.nck,1),
        markerToReconstruct = find(~ismember(1:numel(rigidBodyMarkers),markerTriad.nck(nck,:)));
        reconstructedSolutionsMarkerInds(nck,:) = markerToReconstruct;
        for marker = 1:numel(markerToReconstruct)
            goodTargetVector = sq(rxyz(goodIndex,markerToReconstruct(marker),:)...
                                  -rxyz(goodIndex,markerTriad.nck(nck,2),:));
            solutionBasisCoordinates = rref(cat(2,mtBasis{goodIndex,nck}',goodTargetVector));
            reconstructedSolutions(nck,marker,:,gcount) = solutionBasisCoordinates(:,4);
        end
    end
    gcount = gcount+1;
end


% TRAIN a gaussian mixture model to sort out outliers in the goodIdicies
[reid,remgmModel,remgmLlh] = mixGaussEm([reshape(markerTriad.imd(goodIndicies,:,:),numel(goodIndicies),[]),markerTriad.imo(goodIndicies,:)]',12);

%[reid,remgmModel,remgmLlh] = mixGaussEm(reshape(permute(reconstructedSolutions,[4,1,2,3]),numel(goodIndicies),[])',6);

meanGoodErrorDist = [];
varGoodErrorDist  = [];
goodErrorGroupCount = [];
for i = unique(reid),
    goodErrorGroupCount(i) = sum(reid==i);    
    for nck = 1:size(markerTriad.nck,1),            

    meanGoodErrorDist(nck,i) = mean(...
        nonzeros(...
            triu(...
                sq(sqrt(sum(bsxfun(@minus,...
                                   sq(reconstructedSolutions(nck,1,:,reid==i)),...
                                   permute(sq(reconstructedSolutions(nck,1,:,reid==i)),[1,3,2])).^2)...
                        )...
                   ),...
                1)...
            )...
        );
    varGoodErrorDist(nck,i) = var(...
        nonzeros(...
            triu(...
                sq(sqrt(sum(bsxfun(@minus,...
                                   sq(reconstructedSolutions(nck,1,:,reid==i)),...
                                   permute(sq(reconstructedSolutions(nck,1,:,reid==i)),[1,3,2])).^2)...
                        )...
                   ),...
                1)...
            )...
        );    
    end
end

meanGoodErrorDist = mean(meanGoodErrorDist);
varGoodErrorDist = mean(varGoodErrorDist);
gegid = unique(reid)';
validGegIds = gegid(goodErrorGroupCount>1000);
[meanGoodError ,gid] = max(meanGoodErrorDist(validGegIds)./varGoodErrorDist(validGegIds));
goodErrorEmGmId = validGegIds(gid);
goodErrorDistErrorMean = meanGoodErrorDist(goodErrorEmGmId);
goodErrorDistErrorVar  = varGoodErrorDist(goodErrorEmGmId);
%goodErrorGroupCount(goodErrorEmGmId)




% COMPUTE xyz mean matrix of each nck subset and solution target
reconstructedSolution.mean = cellfun(@(x) mean(sq(x)'),... 
                                     mat2cell(reconstructedSolutions(:,:,:,reid==goodErrorEmGmId),...
                                              ones([size(markerTriad.nck,1),1]),...
                                              ones([numel(markerToReconstruct),1]),...
                                              size(xyz,3),goodErrorGroupCount(goodErrorEmGmId)),...
                                     'UniformOutput',false...
                                     );

% COMPUTE xyz covariance matrix of each nchoosek subset and solution target
reconstructedSolution.var = cellfun(@(x) cov(sq(x)'),... 
                                    mat2cell(reconstructedSolutions(:,:,:,reid==goodErrorEmGmId),...
                                             ones([size(markerTriad.nck,1),1]),...
                                             ones([numel(markerToReconstruct),1]),...
                                             size(xyz,3),goodErrorGroupCount(goodErrorEmGmId)),...
                                    'UniformOutput',false...
                                    );





% $$$ nck = 3; 
% $$$ testInd = 166500; 
% $$$ testInd = 128164; 
% $$$ testInd = 123; 
% $$$ [mtBasis{testInd,nck}*reconstructedSolution.mean{nck}']'
% $$$ [sq(rxyz(testInd,reconstructedSolutionsMarkerInds(nck,1),:)-rxyz(testInd,markerTriad.nck(nck,2),:))]'
% $$$ 
% $$$ rcms =cellfun(@mtimes,mtBasis,repmat(cellfun(@transpose,reconstructedSolution.mean,'UniformOutput',false)',size(mtBasis,1),1),'UniformOutput',false);
% $$$ rcms = cellfun(@transpose,rcms,'UniformOutput',false);
% $$$ 
% $$$ nck = 2;
% $$$ figure,plot( sqrt(sum([cell2mat(rcms(:,nck))...
% $$$                        -sq(rxyz(:,reconstructedSolutionsMarkerInds(nck,1),:)...
% $$$                            -rxyz(:,markerTriad.nck(nck,2),:))].^2,2)))
% $$$ 



% COMPUTE smoothed basis for reconstruction
segSize = 11; % must be odd
embMtBasis = nan([segSize,size(mtBasis)]);
embMtBasis = reshape(mtBasis(bsxfun(@plus,[1:size(mtBasis,1)-segSize]',0:segSize-1),:),...
                     [segSize,size(mtBasis,1)-segSize,size(mtBasis,2)]);
embMtBasis = cellfun(@shiftdim,embMtBasis,repmat({-1},size(embMtBasis)),'UniformOutput',false);

smtBasis = repmat({nan([size(xyz,3),size(xyz,3)])},size(mtBasis));
halfSeg = round(segSize/2);
for i = halfSeg:size(mtBasis)-halfSeg,
    for nck = 1:size(markerTriad.nck,1),
        smtBasis{i,nck} = sq(mean(cell2mat(embMtBasis(:,i-halfSeg+1,nck))));
    end
end


% COMUPUTE the distance b
mrdist = nan([size(mtBasis,1),size(markerTriad.nck,1),1,2]);
sign = [1,-1];
markerTriad.solutions = nan([size(mtBasis,1),size(markerTriad.nck,1),1,size(rxyz,3)]);
markerTriad.oriMarkerNckCoor =nan([size(mtBasis,1),size(markerTriad.nck,1),1,size(rxyz,3)]);
m = 1;
% COMPUTE distance between reconstrution and original marker
for nck = 1:size(markerTriad.nck,1)
    markerTriad.solutions(:,nck,m,:) = ...
        cell2mat(cellfun(@transpose,...
                         cellfun(@(x,y) x*y', smtBasis(:,nck),... %mtBasis(:,nck),...
                                 repmat(reconstructedSolution.mean(nck,m),...
                                        size(mtBasis,1),1),...
                                 'UniformOutput',false),...
                         'UniformOutput',false)...
                 );
    markerTriad.oriMarkerNckCoor(:,nck,m,:) = sq(rxyz(:,reconstructedSolutionsMarkerInds(nck,1),:)...
                                     -rxyz(:,markerTriad.nck(nck,2),:));    
    for s = 1:2        
        mrdist(:,nck,m,s) = sqrt(sum((markerTriad.solutions(:,nck,1,:)...
                                     -sign(s).*markerTriad.oriMarkerNckCoor(:,nck,1,:)).^2,4));
    end
end

figure,
sp=[];
for nck = 1:size(markerTriad.nck,1)
    sp(nck) = subplot(6,1,nck);
    hold on;    
    plot(min(mrdist(:,nck,:),[],3))
end
sp(5) = subplot(6,1,5); plot(reshape(zimdo,size(zimdo,1),[]));
sp(5) = subplot(6,1,5); plot(mean(abs(zimdo),3))
sp(6) = subplot(6,1,6); plot(ang(:,5,7,2)*20)
linkaxes(sp,'xy')


imdo = cat(3,markerTriad.imd,markerTriad.imo);
zimdo = bsxfun(@rdivide,bsxfun(@minus,imdo,mean(imdo(goodIndicies,:,:))),std(imdo(goodIndicies,:,:)));


[~,solSign] = min(mrdist(:,:,1,:),[],4);
[bestTriadMeanZscore,bestTriadInd] = min(max(abs(zimdo),[],3),[],2);
[errorZscore] = max(max(abs(zimdo),[],3),[],2);


% Need head com
rcom = mean(rxyz.data,2);
rcs = sqrt(sum(diff(rcom).^2,3));
% replace marker with best Triad solution
% compair derivatives of oriHcom with solHcom
nxyz = rxyz.copy;
m = 1;
for i = 2:size(rxyz,1),
    bti = bestTriadInd(i);
    %for m = 1,
    if errorZscore(i)>5,
        txyz = sign(solSign(i,bti)).*shiftdim(markerTriad.solutions(i,bti,m,:),1)...
               +rxyz(i,markerTriad.nck(bti,2),:);

        %if sqrt(sum(diff(mean(txyz-nxyz(i-1,reconstructedSolutionsMarkerInds(bti,m),:),2)).^2,3))<rcs(i),
            nxyz.data(i,reconstructedSolutionsMarkerInds(bti,m),:) = txyz;            
            %end


    end
%   end    
end

ncom = mean(nxyz.data,2);
ncs = sqrt(sum(diff(ncom).^2,3));

figure,plot([rcs,ncs])



rcom = mean(rxyz.data,2);
rcs = sqrt(sum(diff(rcom).^2,3));
nxyz = rxyz.copy;
m = 1;
for i = 2:size(rxyz,1),
    for bti = 1:size(markerTriad.nck,1),
       nxyz.data(i,reconstructedSolutionsMarkerInds(bti,m),:) = sign(solSign(i,bti)).*shiftdim(markerTriad.solutions(i,bti,m,:),1)+rxyz(i,markerTriad.nck(bti,2),:);
    end
end
ncom = mean(nxyz.data,2);
ncs = sqrt(sum(diff(ncom).^2,3));

figure,plot([rcs,ncs])




% CREATE an orthogonal coordinate system base on each Triad
newMarkerTriad = compute_marker_triads(nxyz);
% COMPUTE new basis
[~,~,newmtBasis] = cellfun(@svd,cellfun(@squeeze,mat2cell(newMarkerTriad.coor,ones([size(newMarkerTriad.coor,1),1]),ones([size(newMarkerTriad.coor,2),1]),3,3),'UniformOutput',false),'UniformOutput',false);




%% For each error group find best marker swap solution
% test on subset


errorIds = unique(eid)';
errorIds(errorIds==0)=[];

figure,hold on
c = jet(numel(errorIds));
for u = errorIds
    scatter(efet(eid==u,1),efet(eid==u,5),4,c(u,:));
    eidCount(u) = sum(eid==u);
end

disp('')
disp(['Number of error groups: ' num2str(numel(errorIds))]);
disp('')

for i = errorIds,
    switch method
      case 'BEST_SWAP_PERMUTATION'
        % marker swap - all permutations 
        ectry = zeros([eidCount(i),nperm]);
        tectry = zeros([eidCount(i),nperm]);
        bids = eid==i;    
        % Compute residual for each marker swap permutation
        k = 1;
        for c = bperm',
            trxyz.data = rxyz(bids,c,:);
            imori = imo(trxyz);
            
            j = 1;
            for p = fperms,
                tectry(bids,j) = imori(:,p(1),p(2),p(3),p(4),1);
                j = j+1;
            end           
            ectry(bids,k) = var(bsxfun(@minus,tectry(bids,:),mefet),[],2);
            k = k+1;
        end

        
        % Best permutation
        [bpVal,bpInd] =min(mean(ectry));
        
        %%check this bpVal since it changes with each model.

        if bpVal < MarkerSwapErrorThreshold,
            disp('')
            disp(['Error Group: ' num2str(i) ' , bpVal = ' num2str(bpVal)])
            smar = subsref(rigidBody.ml,substruct('()',{bperm(bpInd,:)}));
            [~,rind] = sort(xyz.model.gmi(rigidBody.ml));
            marInd = xyz.model.gmi(smar(rind));
            corInd = sort(marInd);
            xyz.data(eid==i,corInd,:) = xyz(eid==i,marInd,:);
            xyz.save;
        else
            disp('')
            disp(['Error Group: ' num2str(i) ' , bpVal = ' num2str(bpVal)])
            %disp('No options, ignoring error')
            disp(['Ridgid body intermarker distance optimization'])
        end
        
      case 'RIGIDBODY_PARTIAL_RECONSTRUCTION'
    

        
    end

    
end    

disp('Press <any key> to initiate another round of corrections')
disp('Exiting...')
return
%disp('Press <q> to quit')
%waitfor(hfig,'CurrentCharacter')
    





% $$$     lrxyz = rxyz.copy;
% $$$     lrxyz.data = bsxfun(@minus,lrxyz.data,lrxyz(:,2,:));    
% $$$     figure,hold on
% $$$     plotSkeleton(Session,lrxyz,goodIndex);    
% $$$     plot3([0,goodTargetVector(1)],...
% $$$           [0,goodTargetVector(2)],...
% $$$           [0,goodTargetVector(3)])
% $$$     
% $$$     for c = 1:3,
% $$$         ln = plot3(mtBasis{goodIndex,1}(1,c).*5,...
% $$$               mtBasis{goodIndex,1}(2,c).*5,...
% $$$               mtBasis{goodIndex,1}(3,c).*5,...
% $$$               ['.',clrs(c)]);
% $$$         ln.MarkerSize=20;
% $$$     end    
    
