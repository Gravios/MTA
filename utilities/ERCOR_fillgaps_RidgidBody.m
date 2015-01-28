function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model,gind)

%% Testing Vars
Trial = MTATrial('jg05-20120310');
gind = 8000;
%%


xyz = Trial.load('xyz');
rb_model = xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'});
assert(rb_model.N<=6,'MTA:utilities:ERCOR_fillgaps_RidgidBody, number of markers in rb correction must be eq or lt 6');
rb_xyz = xyz.copy;
rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;


BufferSize = 2^16;
NumChunks = mod(rb_xyz.size(1),BufferSize);
LastPiece = (NumChunks*BufferSize+1):rb_xyz.size(1);


imo(rb_xyz(gind,:,:));
bperm = [1:rb_xyz.model.N];
fperms = perms(bperm);
fperms = fperms(:,1:rb_xyz.model.N-1)';


dtgmori = zeros([rb_xyz.size(1),1]);

for c = 0:NumChunks,
    if c~=NumChunks,
        ind = (c*BufferSize+1):(c+1)*BufferSize;
    else
        ind = LastPiece;
    end
    imori = imo(rb_xyz(ind,:,:));

    smori = zeros([BufferSize,size(fperms,2)]);
    i = 1;
    for p = fperms,
        smori(:,i) = imori(:,p(1),p(2),p(3),p(4),1);
        i = i+1;
    end
    
    dtgmori(ind) = var(bsxfun(@minus,smori,gpat),[],2);
end

figure,plot(dtgmori)

badper = ThreshCross(dtgmori,2,5);



%for i = 1:size(badper,1),
%badseg = 
%chkmori = sqrt(sum(bsxfun(@minus,smori,gpat).^2,2));


%interMarDist = imd(xyz);
%figure,imagesc(log10(reshape(interMarDist,[],81)'))

