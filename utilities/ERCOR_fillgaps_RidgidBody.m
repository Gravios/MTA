function xyz = ERCOR_fillgaps_RidgidBody(xyz,rb_model)


Trial = MTATrial('Ed10-20140816');
xyz = Trial.xyz.copy;xyz.load(Trial);
rb_model = xyz.model.rb({'head_back','head_left','head_right','head_front','head_top'});

rb_xyz = xyz.copy;
rb_xyz.data = xyz(:,rb_model.ml,:);
rb_xyz.model = rb_model;

imori = imo(rb_xyz);
figure,plot(imori(:,1,2,2,4,2))