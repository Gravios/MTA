function [xyz,eid,bpIndC] = correct_rigidbody_errors(Session,varargin)
%function xyz = correct_rigidbody_errors(Session, varargin)
% 
% Tries to correct marker swaps and markers which wander from
% a ridgid body constrains
%
% Caution: This function is quick and dirty ...
%
% Note: This function is interactive and depends on the function ClusterPP.
% 
%  VARARGIN 
%
%    Session: MTASession, Session to be processed
%
%    xyz: MTADxyz, xyz object to be processed
%
%    errorSelectionMethod: string, method of selecting period requiring correction
%                                   VALS : 'EMGM'        : gaussian mixture model
%                                          'MANUAL'      : manual error annotation
%                                          'AUTO_THRESH' : automatic detection of error periods
%    
%    correctionMethod: string, method of error correction
%                              VALS - 'BEST_SWAP_PERMUTATION'    : find best fit permutation
%                                     'RIGIDBODY_RECONSTRUCTION' : reconstruct 'noisy' marker 
%                                                        
%
%    goodIndices: double, the indicies used to create a good rigidbody template
%                 empty,  Index selection GUI will appear
%
%    rigidBodyMarkers: cellstr, Markers of rigidbody
%
%    errorThreshold: double, Threshold for switching to auxilary
%                                      Error correction methods.
%
%    bufferSize: double, For machines with low memory (e.g. <4gb)
%
%  Note: the number of markers in rb_model must be <= 6.
%  
%  TODO: Make the assignment of errorThreshold automatic
%        Make memory use automatic (will probably occur upon the rapture :() 
%
%      case 'BEST_SWAP_PERMUTATION'
%      case 'RIGIDBODY_PARTIAL_RECONSTRUCTION'        
%
%  ADD: Input for error periods
%       diagnostic output variables such as error periods and ids
% 
%
% SECTIONS 
%
%   DEFARGS 
%
%   MAIN 
%     S1 : SELECT representative periods
%     S2 : VALIDATE rigidBodyMarkers input type and count
%     S3 : SELECT error periods
%     S4 : CORRECT rigid body errors
%



Session = MTASession.validate(Session);



% DEFARGS ------------------------------------------------------------------------------------------
defArgs = struct('xyz',                   [],                                                    ...
                 'errorSelectionMode',    'MANUAL',                                              ...
                 'correctionMethod',      'BEST_SWAP_PERMUTATION',                               ...
                 'goodIndicies',          [],                                                    ...
                 'rigidBodyMarkers',      {{'head_back','head_left','head_right',                ...
                                            'head_front','head_top'}},                           ...
                 'errorThreshold',        0.01,                                                  ...
                 'bufferSize',            2^16                                                   ...
);
[xyz,errorSelectionMode,correctionMethod,goodIndicies,rigidBodyMarkers,errorThreshold,           ...
 bufferSize] = DefaultArgs(varargin,defArgs,'--struct');

if isempty(xyz),  xyz = Session.load('xyz');  end
numChunks = floor(size(xyz,1)./bufferSize);
lastPiece = (numChunks*bufferSize+1):size(xyz,1);
% END DEFARGS --------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------
% S1 : SELECT representative periods ---------------------------------------------------------------
if isempty(goodIndicies),
% SELECT periods without errors if no good index is specified
    pfig = PlotSessionErrors(Session);
    disp(['Please select a point on the x-axis with the data cursor ' ...
          'and then press enter']);
    ylm = ylim;
    bflag = false;
% MANUALLY select periods where rigid body is correct    
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
% END S1 -------------------------------------------------------------------------------------------




% S2 : VALIDATE rigidBodyMarkers input type and count ----------------------------------------------
if iscell(rigidBodyMarkers)&&~isa(rigidBodyMarkers,'MTAModel'),
    rigidBody = xyz.model.rb(rigidBodyMarkers);
end
assert(rigidBody.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);
assert(rigidBody.N>3,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be gt 3']);
% END S2 -------------------------------------------------------------------------------------------




% S3 : SELECT error periods ------------------------------------------------------------------------

% COPY xyz object
rxyz = xyz.copy;

% SELECT rigidbody marker subset
rxyz.data = xyz(:,rigidBody.ml,:);
rxyz.model = rigidBody;

% GET all permutations of rigid body markers
nperms = factorial(rxyz.model.N);
bperms  = mat2cell(perms(1:rxyz.model.N),ones([nperms,1]),ones([rxyz.model.N,1]));
fperms = bperms(:,1:rxyz.model.N-1);

% CREATE 1d error feature based on intermarker orientaions
efet = zeros([size(rxyz,1),nperms]);
for c = 0:numChunks,
    if c~=numChunks,
        ind = (c*bufferSize+1):(c+1)*bufferSize;
    else
        ind = lastPiece;
    end
    imori = imo(rxyz(ind,:,:));
    for p = 1:size(fperms,1),
        efet(ind,p) = imori(:,fperms{p,:},1);
    end
end
referenceErrorFeatureMean = mean(efet(goodIndicies,:));
referenceErrorFeatureVar = var(efet(goodIndicies,:));
dtgmori = var(bsxfun(@minus,efet,referenceErrorFeatureMean),[],2);
%dtgmori = log10(sum(abs(bsxfun(@ldivide,bsxfun(@minus,efet,referenceErrorFeatureMean),referenceErrorFeatureVar)),2));

% SELECT Good periods if none are provided
switch errorSelectionMode
  case 'AUTO_THRESH'
    nind = nniz(efet);
    eper = ThreshCross(dtgmori,0.2,round(xyz.sampleRate/8));    
    eid = zeros(size(nind))';    
    for e = 1:size(eper,1),
        errorInd = zeros(size(nind));
        errorInd(eper(e,1):eper(e,2)) = true;
        eid(nind&errorInd)=e;        
    end
    
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
% END S3 -------------------------------------------------------------------------------------------



% S4 : CORRECT rigid body errors -------------------------------------------------------------------
% For each error group find best marker swap solution

errorIds = unique(eid)';
errorIds(errorIds==0)=[];

% CREATE figure to visualize error groups
figure();
hold('on');
c = jet(numel(errorIds));
emptyErrIds = [];
for u = errorIds'
    if sum(efet(eid==u,1))>1,
        scatter(efet(eid==u,1),efet(eid==u,5),4,c(u,:));
        eidCount(u) = sum(eid==u);
    else
        emptyErrIds(end+1) = u;
    end
end
errorIds(ismember(errorIds,emptyErrIds)) = [];

disp('')
disp(['Number of error groups: ' num2str(numel(errorIds))]);
disp('')
bpIndC = [];
for i = errorIds',
    tic        
% $$$     switch correctionMethod
% $$$       case 'BEST_SWAP_PERMUTATION', % marker swap - all permutations 
        if eidCount(i)==1, continue, end
        ectry  = zeros([eidCount(i),nperms]);
        tectry = zeros([eidCount(i),nperms]);
        bids = eid==i;    
% COMPUTE residual for each marker swap permutation

        for b = 1:size(bperms,1),
            imori = imo(rxyz(bids,[bperms{b,:}],:));
            for p = 1:size(fperms,1),
                tectry(bids,p) = imori(:,fperms{p,:},1);
            end           
            ectry(bids,b) = var(bsxfun(@minus,...
                                       tectry(bids,:),...
                                       referenceErrorFeatureMean),[],2);
        end

        % Best permutation
        [bpVal,bpInd] =min(mean(ectry));
        
        bpIndC(i) = bpInd;
        if bpVal < errorThreshold,
        %check this bpVal since it changes with each model.
            disp('')
            disp(['Error Group: ' num2str(i) ' , bpVal = ' num2str(bpVal)])
            smar = subsref(rigidBody.ml,substruct('()',{[bperms{bpInd,:}]}));
            [~,rind] = sort(xyz.model.gmi(rigidBody.ml));
            marInd = xyz.model.gmi(smar(rind));
            corInd = sort(marInd);
            xyz.data(eid==i,corInd,:) = xyz(eid==i,marInd,:);
        else
            disp('')
            disp(['Error Group: ' num2str(i) ' , bpVal = ' num2str(bpVal)])
            %disp('No options, ignoring error')
            disp(['Ridgid body intermarker distance optimization'])
        end
        
% $$$       case 'RIGIDBODY_RECONSTRUCTION'
% $$$           disp('RIGIDBODY_RECONSTRUCTION under contruction')
        
% $$$     end% switch correctionMethod
    toc        
end% for i

xyz.save;


% END S4 -------------------------------------------------------------------------------------------



disp('Press <any key> to initiate another round of corrections')
disp('Exiting...')
return
% END MAIN -----------------------------------------------------------------------------------------




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
    
