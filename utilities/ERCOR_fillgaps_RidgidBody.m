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

if ischar(Session),
    Session = MTASession(Session);
end


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
    disp(['Please select a point on the x-axis with the data cursor ' ...
          'and then press enter']);
    waitfor(pfig,'CurrentCharacter',char(13));
    good_index = dcm_obj.getCursorInfo.Position(1);
    disp(['Index = ' num2str(good_index)])
    delete(pfig);
end

% Check and sanitize rb_model input type
if iscell(rb_model)&&~isa(rb_model,'MTAModel'),
    rb_model = xyz.model.rb(rb_model);
end

assert(rb_model.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);

hfig = figure(3848283);
hfig.CurrentCharacter = ' ';

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
    clf(hfig)
    plot(dtgmori)
    disp('******************************');
    disp('* Select each error group:   *')
    disp('*   add point:   left click  *')
    disp('*   close group: right click *')
    disp('*   quit:        escape key  *')
    disp('******************************'); 
    
    eid = ClusterPP(hfig);
    
    
    %% given errors find best marker swap solution
    % test on subset
    nSamp = 100;
    
    errors = unique(eid);
    ectry = zeros([nSamp,nperm]);
    tectry = zeros([nSamp,nperm]);
   
    disp('')
    disp(['Number of error groups: ' num2str(numel(errors)-1)]);
    disp('')
    
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
            
            %warning(['ERCOR_fillgaps_RidgidBody:MarkerSwapFailed, ' ...
            %         'procceding to attemtep marker reconstruction from rigid body']);
            
            
            % $$$         drb_xyz = imd(rb_xyz);
            % $$$         ndrb_xyz = bsxfun(@minus,drb_xyz,drb_xyz(good_index,:,:));
            % $$$
            % $$$         log10(abs(mean(drb_xyz(eid==2,:,:))))
            % $$$         gperms = find(abs(mean(efet(bids,:))-gpat)<.2);
            
            % Try some other correction method
            % probably by reconstructing the bad
            % markers position using the other markers

% $$$             dectry = sum(abs(bsxfun(@minus,reshape(im,size(im,1),[]),reshape(im(68000,:,:),1,[]))),2);
% $$$             
% $$$             im =imd(rb_xyz);
% $$$             imm = bsxfun(@minus,im,im(good_index,:,:));
% $$$             gm = {};
% $$$             erind = find(eid==i);
% $$$             min_dist_thresh = 1;
% $$$             for c = erind,
% $$$             for m= 1:5,
% $$$                 gm{m} = find(abs(imm(erind,m,[1:m-1,m+1:5]))<min_dist_thresh);
% $$$             end
% $$$             gm = unique(cell2mat(gm(~cellfun(@isempty,gm))'))
% $$$             % gm 
% $$$             end
            
        end
        
    end
    disp('Press <any key> to initiate another round of corrections')
    disp('Press <q> to quit')
    waitfor(hfig,'CurrentCharacter')
    
end

