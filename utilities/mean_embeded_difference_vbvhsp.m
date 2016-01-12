function [mz,vz] = mean_embeded_difference_vbvhsp(Data,Trial)

NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
ANG_HISTOGRAM_BOUNDARIES = linspace(-pi/2,pi/2,NBINS);

xyz = Trial.load('xyz');
xyz.addMarker('bcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

ang = create(MTADang,Trial,xyz);
vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
ind = Trial.stc('a-m-r').cast('TimeSeries');
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_p_s] = histc(ang(:,3,4,2),ANG_HISTOGRAM_BOUNDARIES);
mind = nniz([ind_v_b,ind_v_h,ind_p_s])&ind.data;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind),ind_p_s(mind)];

acvar = Data(mind,:,:,:);
mz = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanmean);
vz = accumarray(manifoldIndex,acvar,[NBINS,NBINS,NBINS],@nanstd);

