function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model)


Trial = MTATrial('jg05-20120310');
xyz = Trial.load('xyz');
rb_model = xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'});
assert(rb_xyz.model.N<=6,'MTA:utilities:ERCOR_fillgaps_RidgidBody, number of markers in rb correction must be eq or lt 6');
rb_xyz = xyz.copy;
rb_xyz.data = xyz(1:50000,rb_model.ml,:);
rb_xyz.model = rb_model;


imori = imo(rb_xyz);



bperm = [1:rb_xyz.model.N];
fperms = perms(bperm);
fperms = fperms(:,1:4)';

smori = zeros([size(imori,1),size(fperms,2)]);
for p = fperms,
smori(:,i) = imori(:,p(1),p(2),p(3),p(4),1);
i = i+1;
end
clear('imori');

gind = 8000;
gpat = smori(gind,:);
epat = smori(8001,:);

dtgmori = sqrt(sum(bsxfun(@minus,smori,gpat).^2,2));

figure,plot((dtgmori))

badper = ThreshCross(dtgmori,2,5);



for i = 1:size(badper,1),
badseg = 
chkmori = sqrt(sum(bsxfun(@minus,smori,gpat).^2,2));


interMarDist = imd(xyz);
figure,imagesc(log10(reshape(interMarDist,[],81)'))
