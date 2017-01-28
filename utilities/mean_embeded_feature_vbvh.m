function [mz,vz,cz,czDom] = mean_embeded_feature_vbvh(Data,Trial,varargin)

[dims,featureDomainBoundaries] = DefaultArgs(varargin,{[],[]});

NBINS = 25;
VEL_HISTOGRAM_BOUNDARIES = linspace(-2,2,NBINS);

xyz = Trial.load('xyz');
xyz.resample(Data);
xyz.addMarker('bcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

% Get the bin index for each marginal distribution
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);

% Select only non-zero points whithout rearing and grooming
if ~isempty(regexpi(Trial.stc.mode,'^hand_labeled.*')),
    ind = Trial.stc('a-m-r').cast('TimeSeries');
    ind.resample(Data);
else
    ang = create(MTADang,Trial,xyz);    
    ind = Trial.stc('a').cast('TimeSeries');
    ind.resample(Data);
    ind.data(isnan(ind.data))=0;
    % Need some hard thresholds for subspace selection 
    bang = ang(:,'spine_lower','spine_upper',3).*cos(ang(:,'spine_lower','spine_upper',2));
    bang = MTADxyz('data',bang/prctile(bang(nniz(bang)),95),'sampleRate',xyz.sampleRate);
    bang.data(~nniz(bang.data))=0;
    ind.data = ind.data&bang>.8;
end

mind = nniz([ind_v_b,ind_v_h])&ind.data==1;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind)];

mz = cell([1,Data.size(2)]);
vz = cell([1,Data.size(2)]);
cz = cell([1,Data.size(2)]);
czDom = cell([1,Data.size(2)]);
for dind = dims(:)',
% $$$     if isempty(featureDomainBoundaries),
% $$$         dataHistogramEdges = linspace([prctile(Data(logical(ind.data),dind),[5,95]),50]);
% $$$     else        
% $$$         dataHistogramEdges = linspace([featureDomainBoundaries(dind==dims,:),50]);
% $$$     end
% $$$     czDom{dind} = dataHistogramEdges;
% $$$     [~,ind_f_d] = histc(Data(:,dind),dataHistogramEdges);    
% $$$     mind = nniz([ind_v_b,ind_v_h,ind_f_d])&ind.data==1;
% $$$     manifoldIndex = [ind_v_b(mind),ind_v_h(mind)];

    acvar = Data(mind,dind);    
% $$$     cz{dind} = accumarray([manifoldIndex,ind_f_d(mind)],ones([sum(mind),1]),[repmat(NBINS,[1,size(manifoldIndex,2)]),50],@sum);
    mz{dind} = accumarray(manifoldIndex,acvar,repmat(NBINS,[1,size(manifoldIndex,2)]),@nanmean);
    vz{dind} = accumarray(manifoldIndex,acvar,repmat(NBINS,[1,size(manifoldIndex,2)]),@nanstd);
end
