function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model,gind)
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
%    rb_model: MTAModel, model which contains the markers with a
%                        ridgid body constraint.
%
%    gind: double, the index which servers as a "good" template.
%
% Note: the number of markers in rb_model must be <= 6.
%
%% Testing Vars
Session = MTASession('jg05-20120310');
gind = 8000;
% $$$ Session = MTASession('jg05-20120317');
% $$$ gind = 314500;
xyz = Session.load('xyz');
rb_model = xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'});
%%

if iscell(rb_model)&&~isa(rb_model,'MTAModel'),
    rb_model = xyz.model.rb(rb_model);
end


assert(rb_model.N<=6,['MTA:utilities:ERCOR_fillgaps_RidgidBody, ' ...
                      'number of markers in rb correction must be eq or lt 6']);

rb_xyz = xyz.copy;
trb_xyz = rb_xyz.copy;

rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;

BufferSize = 2^16;
NumChunks = floor(rb_xyz.size(1)./BufferSize);
LastPiece = (NumChunks*BufferSize+1):rb_xyz.size(1);

nperm = factorial(rb_xyz.model.N);
bperm = [1:rb_xyz.model.N];
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

dfet = bsxfun(@minus,efet,efet(gind,:));
dtgmori = var(bsxfun(@minus,efet,efet(gind,:)),[],2);



figure,imagesc(dfet(1:100000,:)');
figure,  plot(dfet( 7000:10000,68),dfet( 7000:10000,69),'b.');
hold on, plot(dfet(13000:14000,68),dfet(13000:14000,69),'r.');

%% Manually select error groups
hfig = figure(39293);
plot(dtgmori)
eid = ClusterPP(hfig);


%% given errors find best marker swap solution
% test on subset
nSamp = 100;

nerrors = unique(eid)-1;
ectry = zeros([nSamp,nperm]);
tectry = zeros([nSamp,nperm]);

for i = 1:nerrors,
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
        ectry(:,k) = var(bsxfun(@minus,tectry,efet(gind,:)),[],2);
        k = k+1;
    end

    % Best permutation
    [bpVal,bpInd] =min(mean(ectry));
    
    %%check this bpVal since it changes with each model.
    if bpVal<0.01,
        smar = subsref(rb_model.ml,substruct('()',{bperm(bpInd,:)}));
        [~,rind] = sort(xyz.model.gmi(rb_model.ml));
        marInd = xyz.model.gmi(smar(rind));
        corInd = sort(marInd);
        xyz.data(eid==i,corInd,:) = xyz(eid==i,marInd,:);
    else
        drb_xyz = imd(rb_xyz);
        ndrb_xyz = bsxfun(@minus,drb_xyz,drb_xyz(gind,:,:));

        log10(abs(mean(drb_xyz(eid==2,:,:))))
        gperms = find(abs(mean(efet(bids,:))-gpat)<.2);

        % Try some other correction method
        % probably by reconstructing the bad 
        % markers position using the other markers
    end
    
end

