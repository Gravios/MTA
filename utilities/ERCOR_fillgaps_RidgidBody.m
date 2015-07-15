function xyz = ERCOR_fillgaps_RidgidBody(Session,varargin)
%function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model,gind)
% 
% Tries to correct marker swaps and markers which wander from
% a ridgid body constrains
%
% Note: This function is interactive and depends on the function ClusterPP.
%
% variables:
%
%    xyz: MTADxyz, motion capture "mocap" data 
%
%  varargin:
%
%    good_index: double, the index which servers as a "good" template.
%                        if good_index is empty a gui will help you out :)
%
%    rb_model: MTAModel, model which contains the markers with a
%                        ridgid body constraint.
%
%    MarkerSwapErrorThreshold: double, Threshold for switching to auxilary
%                                      Error correction methods.
%
%    BufferSize: double, For machines with low memory (e.g. <4gb)
%
%  Note: the number of markers in rb_model must be <= 6.
%  
%  TODO: Make the assignment of MarkerSwapErrorThreshold automatic
%        Make memory use automatic (will probably occur upon the rapture :() 
%

xyz = Session.load('xyz');

defArgs = ...
{... good_index
     [],...
 ... rb_model
     xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'}),...
 ... MarkerSwapErrorThreshold
     0.01,...
 ... BufferSize
     2^16 ...
};
 
% Load Default Arguments
[good_index,rb_model,MarkerSwapErrorThreshold,BufferSize] = ...
    DefaultArgs(varargin,defArgs);

NumChunks = floor(xyz.size(1)./BufferSize);
LastPiece = (NumChunks*BufferSize+1):xyz.size(1);

% Manual specify a good index with the aid of a gui
if isempty(good_index),
    pfig = PlotSessionErrors(Session);
    dcm_obj = datacursormode(gcf);
    waitfor(pfig,'CurrentCharacter',char(13));
    good_index = dcm_obj.getCursorInfo.Position(1);
    delete(pfig);
end

% Check and sanitize rb_model input type
if iscell(rb_model)&&~isa(rb_model,'MTAModel'),
    rb_model = xyz.model.rb(rb_model);
end

assert(rb_model.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);

hfig = figure(3848283);

while ~strcmp(get(hfig,'CurrentCharacter'),'q'),
    
    rb_xyz = xyz.copy;
    trb_xyz = rb_xyz.copy;
    
    rb_xyz.data = xyz(:,rb_model.ml,:);
    rb_xyz.model = rb_model;
        
    nperm = factorial(rb_xyz.model.N);
    bperm = 1:rb_xyz.model.N;
    bperm = perms(bperm);
    fperms = bperm(:,1:rb_xyz.model.N-1)';
    
    
    efet = zeros([rb_xyz.size(1),nperm]);
    for c = 0:NumChunks,
        if c~=NumChunks,
            ind = (c*BufferSize+1):(c+1)*BufferSize;
        else
            ind = LastPiece;
        end
        trb_xyz.data = rb_xyz(ind,:,:);
        imori = imo(trb_xyz);

        i = 1;
        for p = fperms,
            efet(ind,i) = imori(:,p(1),p(2),p(3),p(4),1);
            i = i+1;
        end
    end
    dtgmori = var(bsxfun(@minus,efet,efet(good_index,:)),[],2);
    
    %% Manually select error groups
    plot(dtgmori)
    eid = ClusterPP(hfig);
    
    
    %% given errors find best marker swap solution
    % test on subset
    nSamp = 100;
    
    errors = unique(eid);
    ectry = zeros([nSamp,nperm]);
    tectry = zeros([nSamp,nperm]);
    
    for i = errors(2:end),
        bids = find(eid==i,nSamp,'first');
        
        k = 1;
        for c = bperm',
            trb_xyz.data = rb_xyz(bids,c,:);
            imori = imo(trb_xyz);
            
            j = 1;
            for p = fperms,
                tectry(:,j) = imori(:,p(1),p(2),p(3),p(4),1);
                j = j+1;
            end
            ectry(:,k) = var(bsxfun(@minus,tectry,efet(good_index,:)),[],2);
            k = k+1;
        end
        
        % Best permutation
        [bpVal,bpInd] =min(mean(ectry));
        
        %%check this bpVal since it changes with each model.
        if bpVal < MarkerSwapErrorThreshold,
            smar = subsref(rb_model.ml,substruct('()',{bperm(bpInd,:)}));
            [~,rind] = sort(xyz.model.gmi(rb_model.ml));
            marInd = xyz.model.gmi(smar(rind));
            corInd = sort(marInd);
            xyz.data(eid==i,corInd,:) = xyz(eid==i,marInd,:);
            xyz.save;
        else
            warning('ERCOR_fillgaps_RidgidBody:MarkerSwapFailed, procceding to attemtep marker reconstruction from rigid body');
            
            % $$$         drb_xyz = imd(rb_xyz);
            % $$$         ndrb_xyz = bsxfun(@minus,drb_xyz,drb_xyz(good_index,:,:));
            % $$$
            % $$$         log10(abs(mean(drb_xyz(eid==2,:,:))))
            % $$$         gperms = find(abs(mean(efet(bids,:))-gpat)<.2);
            
            % Try some other correction method
            % probably by reconstructing the bad
            % markers position using the other markers
        end
        
    end
    
end

