function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model)


Trial = MTATrial('Ed10-20140816');
xyz = Trial.xyz.copy;xyz.load(Trial);
rb_model = xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'});

rb_xyz = xyz.copy;
rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;

imori = imo(rb_xyz);

i = 1;
fperms = perms(1:4)';
smori = zeros([size(imori,1),size(fperms,2)]);
for p = fperms,
smori(:,i) = imori(:,p(1),p(2),p(3),p(4),1);
i = i+1;
end
clear('imori');

%figure,imagesc(unity(smori)');
%caxis([-.6,.6])

smori = unity(smori);
gind = 8000;
gpat = smori(gind,:);
epat = smori(8001,:);
dtgmori = sqrt(sum(bsxfun(@minus,smori,gpat).^2,2));

interMarDist = imd(xyz);
figure,imagesc(log10(reshape(interMarDist,[],81)'))
