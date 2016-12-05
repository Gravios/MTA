function xyz = ERCOR_fillgaps_RidgidBody(Session,varargin)
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



% DEFARGS ------------------------------------------------------------------------------------------
defArgs = struct('mode',      'EMGM',                                                            ...
                 'method',    'BEST_SWAP_PERMUTATION',                                           ...
                 'goodIndicies', [],                                                             ...
                 'rb_model',  {{'head_back','head_left','head_right','head_front','head_top'}},  ...
                 'MarkerSwapErrorThreshold',0.01,                                                ...
                 'bufferSize',              2^16                                                 ...
);

[mode,method,goodIndex,rb_model,MarkerSwapErrorThreshold,bufferSize] = DefaultArgs(varargin,defArgs,'--struct');

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
% if no good index is specified select one with the aid of a gui
if isempty(goodIndicies),
    pfig = PlotSessionErrors(Session);
    disp(['Please select a point on the x-axis with the data cursor ' ...
          'and then press enter']);
    ylm = ylim;
    while ~strcmp(get(pfig,'CurrentCharacter'),char(27));
        set(pfig,'CurrentCharacter',char(32));        
        while ~strcmp(get(pfig,'CurrentCharacter'),char(13))
            waitforbuttonpress
            if strcmp(get(pfig,'CurrentCharacter'),char(27))
                bflag = 1;
                break
            end
        end
        if bflag, break,end
        rect = getrect(pfig);
        patch([rect(1),rect(1),rect(1)+rect(3),rect(1)+rect(3)],...
              [ylm(1),ylm(2),ylm(2),ylm(1)],...
              [-1,-1,-1,-1],[.9,.9,.9]);
        goodIndicies = cat(1,goodIndicies,[round(rect(1)):round(rect(1)+rect(3))]');
    end
    goodIndicies = unique(goodIndicies);
    
    disp(['Number of indicies selected: = ' num2str(numel(goodIndicies))])
    delete(pfig);
end

% Check and sanitize rb_model input type
if iscell(rb_model)&&~isa(rb_model,'MTAModel'),
    rb_model = xyz.model.rb(rb_model);
end
assert(rb_model.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);
assert(rb_model.N>3,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be gt 3']);



rb_xyz = xyz.copy;
trb_xyz = rb_xyz.copy;

rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;

nperm = factorial(rb_xyz.model.N);

% Get all permutations of 4 markers
bperm = 1:rb_xyz.model.N;
bperm = perms(bperm);
fperms = bperm(:,1:4)';



switch method
  case 'RIGIDBODY_PARTIAL_RECONSTRUCTION'
    headRigidBodyMarkers={'head_back','head_left','head_front','head_right'};
    headRigidBody = xyz.model.rb(headRigidBodyMarkers);
    hcom = xyz.com(headRigidBody);
    hxyz = xyz.copy;
    hxyz.model = headRigidBody;
    hxyz.data = xyz(:,headRigidBodyMarkers,:);
    markerIndNCK = nchoosek(1:size(hxyz,2),3);
    markerTrioCOM = nan([size(hxyz,1),1,size(hxyz,3)]);
    markerTrioCOOR = nan([size(hxyz,1),size(markerIndNCK,1),size(hxyz,3),size(hxyz,3)]);
    for nck = 1:size(markerIndNCK,1),        
        markerTrioCOOR(:,nck,[1,2],:) = permute(bsxfun(@minus,hxyz(:,markerIndNCK(nck,[1,3]),:),hxyz(:,markerIndNCK(nck,2),:)),[1,4,2,3]);
        markerTrioCOOR(:,nck,3,:) = cross(markerTrioCOOR(:,nck,1,:),markerTrioCOOR(:,nck,2,:));
        markerTrioCOOR(:,nck,2,:) = cross(markerTrioCOOR(:,nck,1,:),markerTrioCOOR(:,nck,3,:));
        markerTrioCOM(:,nck,:) = mean(hxyz(:,markerIndNCK(nck,:),:),2);
    end
    
    [~,~,mtEigenVector] = cellfun(@svd,cellfun(@squeeze,mat2cell(markerTrioCOOR,ones([size(markerTrioCOOR,1),1]),ones([size(markerTrioCOOR,2),1]),3,3),'UniformOutput',false),'UniformOutput',false);

    reconstructedSolutions  = nan([size(markerIndNCK,1),3,headRigidBody.N-3,numel(goodIndicies)]);
    reconstructedSolutionsMarkerInds  = nan([size(markerIndNCK,1),headRigidBody.N-3]);
    gcount = 1;
    for goodIndex = goodIndicies(:)'
        for nck = 1:size(markerIndNCK,1),            
            markerToReconstruct = find(~ismember(1:numel(headRigidBodyMarkers),markerIndNCK(nck,:)));
            reconstructedSolutionsMarkerInds(nck,:) = markerToReconstruct;
            for marker = 1:numel(markerToReconstruct)
                goodTargetVector = sq(hxyz(goodIndex,markerToReconstruct(marker),:)...
                                      -hxyz(goodIndex,markerIndNCK(nck,2),:));
                solutionBasisCoordinates = rref(cat(2,mtEigenVector{goodIndex,1},goodTargetVector));
                reconstructedSolutions(nck,:,marker,gcount) = solutionBasisCoordinates(:,4);
            end
        end
        gcount = gcount+1;
    end
    
    
    
    
    
end


% Calculate intermarker orientaions
efet = zeros([rb_xyz.size(1),nperm]);
for c = 0:numChunks,
    if c~=numChunks,
        ind = (c*bufferSize+1):(c+1)*bufferSize;
    else
        ind = lastPiece;
    end
    trb_xyz.data = rb_xyz(ind,:,:);
    imori = imo(trb_xyz);

    i = 1;
    for p = fperms,
        efet(ind,i) = imori(:,p(1),p(2),p(3),p(4),1);
        i = i+1;
    end
end

dtgmori = var(bsxfun(@minus,efet,efet(goodIndex,:)),[],2);

switch mode
  case 'EMGM'
    % em-gaussian mixture model
    nind = nniz(efet)&log10(dtgmori)>-0.5;
    [teid,emgmModel,emgmLlh] = mixGaussEm(efet(nind,:)',30);
    eid = zeros(size(nind));
    eid(nind)=teid;
  case 'MANUAL'
    % Use gui to select error periods and attempt to correct
    % Create 1-d distance metric for detecting errors
    % Manually select error groups
    hfig = figure(3848283);
    hfig.CurrentCharacter = ' ';
    
    clf(hfig)
    plot(dtgmori)
    disp('******************************');
    disp('* Select each error group:   *')
    disp('*   add point:   left click  *')
    disp('*   close group: right click *')
    disp('*   quit:        escape key  *')
    disp('******************************'); 
    
    eid = ClusterPP(hfig);
        
  otherwise
    error('MTA:utilities:ERCOR_fillgaps_RidgidBody:ModeNotFound');
end






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
            trb_xyz.data = rb_xyz(bids,c,:);
            imori = imo(trb_xyz);
            
            j = 1;
            for p = fperms,
                tectry(bids,j) = imori(:,p(1),p(2),p(3),p(4),1);
                j = j+1;
            end           
            ectry(bids,k) = var(bsxfun(@minus,tectry(bids,:),efet(goodIndex,:)),[],2);
            k = k+1;
        end

        
        % Best permutation
        [bpVal,bpInd] =min(mean(ectry));
        
        %%check this bpVal since it changes with each model.

        if bpVal < MarkerSwapErrorThreshold,
            disp('')
            disp(['Error Group: ' num2str(i) ' , bpVal = ' num2str(bpVal)])
            smar = subsref(rb_model.ml,substruct('()',{bperm(bpInd,:)}));
            [~,rind] = sort(xyz.model.gmi(rb_model.ml));
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
    





% $$$     lhxyz = hxyz.copy;
% $$$     lhxyz.data = bsxfun(@minus,lhxyz.data,lhxyz(:,2,:));    
% $$$     figure,hold on
% $$$     plotSkeleton(Session,lhxyz,goodIndex);    
% $$$     plot3([0,goodTargetVector(1)],...
% $$$           [0,goodTargetVector(2)],...
% $$$           [0,goodTargetVector(3)])
% $$$     
% $$$     for c = 1:3,
% $$$         ln = plot3(mtEigenVector{goodIndex,1}(1,c).*5,...
% $$$               mtEigenVector{goodIndex,1}(2,c).*5,...
% $$$               mtEigenVector{goodIndex,1}(3,c).*5,...
% $$$               ['.',clrs(c)]);
% $$$         ln.MarkerSize=20;
% $$$     end    
    
