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
                 'rigidBodyMarkers',  {{'head_back','head_left','head_right','head_front','head_top'}},  ...
                 'MarkerSwapErrorThreshold',0.01,                                                ...
                 'bufferSize',              2^16                                                 ...
);

[mode,method,goodIndicies,rigidBodyMarkers,MarkerSwapErrorThreshold,bufferSize] = DefaultArgs(varargin,defArgs,'--struct');

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
    % Manually select periods where rigid body is correct
    pfig = PlotSessionErrors(Session);
    disp(['Please select a point on the x-axis with the data cursor ' ...
          'and then press enter']);
    ylm = ylim;
    bflag = false;
    while ~strcmp(get(pfig,'CurrentCharacter'),char(27));
        set(pfig,'CurrentCharacter',char(32));        
        while ~strcmp(get(pfig,'CurrentCharacter'),char(13))
            waitforbuttonpress
            if strcmp(get(pfig,'CurrentCharacter'),char(27))
                bflag = true;
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



% Check and sanitize rigidBodyMarkers input type
if iscell(rigidBodyMarkers)&&~isa(rigidBodyMarkers,'MTAModel'),
    rigidBody = xyz.model.rb(rigidBodyMarkers);
end
assert(rigidBody.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);
assert(rigidBody.N>3,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be gt 3']);


rxyz = xyz.copy;
trxyz = rxyz.copy;

rxyz.data = xyz(:,rigidBody.ml,:);
rxyz.model = rigidBody;

nperm = factorial(rxyz.model.N);

% Get all permutations of 4 markers
bperm = 1:rxyz.model.N;
bperm = perms(bperm);
fperms = bperm(:,1:4)';



% CREATE 1d error feature based on intermarker orientaions
efet = zeros([rxyz.size(1),nperm]);
for c = 0:numChunks,
    if c~=numChunks,
        ind = (c*bufferSize+1):(c+1)*bufferSize;
    else
        ind = lastPiece;
    end
    trxyz.data = rxyz(ind,:,:);
    imori = imo(trxyz);

    i = 1;
    for p = fperms,
        efet(ind,i) = imori(:,p(1),p(2),p(3),p(4),1);
        i = i+1;
    end
end
mefet = mean(efet(goodIndicies,:));
vefet = var(efet(goodIndicies,:));
dtgmori = var(bsxfun(@minus,efet,mefet),[],2);
%dtgmori = log10(sum(abs(bsxfun(@ldivide,bsxfun(@minus,efet,mefet),vefet)),2));

% SELECT Good periods if none are provided
switch mode
  case 'EMGM'
    % em-gaussian mixture model
    nind = nniz(efet)&log10(dtgmori)>-0.5;
    [teid,emgmModel,emgmLlh] = mixGaussEm(efet(nind,:)',40);
    eid = zeros(size(nind))';
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
for u = errorIds'
    scatter(efet(eid==u,1),efet(eid==u,5),4,c(u,:));
    eidCount(u) = sum(eid==u);
end

disp('')
disp(['Number of error groups: ' num2str(numel(errorIds))]);
disp('')

for i = errorIds',
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
    
