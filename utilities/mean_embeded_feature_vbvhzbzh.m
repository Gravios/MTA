function [mz,vz,cz,czDom] = mean_embeded_feature_vbvhzbzh(Data,Trial,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
[dims,featureDomainBoundaries,verbose] = DefaultArgs(varargin,{[],[],false},true);

NBINS = 20;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
HEIGHT_HISTOGRAM_BOUNDARIES = linspace(0,200,NBINS);

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');
xyz.resample(Data);

vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
% Get the bin index for each marginal distribution
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_z_b] = histc(xyz(:,'bcom',3),HEIGHT_HISTOGRAM_BOUNDARIES);
[~,ind_z_h] = histc(xyz(:,'hcom',3),HEIGHT_HISTOGRAM_BOUNDARIES);

% SELECT only non-zero points whithout rearing and grooming
if ~isempty(regexpi(Trial.stc.mode,'^hand_labeled.*')),
    ind = Trial.stc('a-m-r').cast('TimeSeries');
    ind.resample(Data);
else
    ang = create(MTADang,Trial,xyz);    
    ind = Trial.stc('a').cast('TimeSeries');
    ind.resample(Data);
    ind.data(isnan(ind.data))=0;
% REQUIRE a hard thresholds to remove groom and rear like periods
    bang = ang(:,'spine_lower','spine_upper',3).*cos(ang(:,'spine_lower','spine_upper',2));
    bang = MTADxyz('data',bang/prctile(bang(nniz(bang)),95),'sampleRate',xyz.sampleRate);
    bang.data(~nniz(bang.data))=0;
    ind.data = ind.data&bang>.8;
end

mind = nniz([ind_v_b,ind_v_h,ind_z_b,ind_z_h])&ind.data==1;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind),ind_z_b(mind),ind_z_h(mind)];

mz = cell([1,Data.size(2)]);
vz = cell([1,Data.size(2)]);
cz = cell([1,Data.size(2)]);
czDom = cell([1,Data.size(2)]);
if verbose, hfig = figure(gen_figure_id); end
for dind = dims(:)',
    acvar = Data(mind,dind);    
    mz{dind} = accumarray(manifoldIndex,acvar,repmat(NBINS,[1,size(manifoldIndex,2)]),@nanmean);
    vz{dind} = accumarray(manifoldIndex,acvar,repmat(NBINS,[1,size(manifoldIndex,2)]),@nanstd);
    cz{dind} = accumarray(manifoldIndex,acvar,repmat(NBINS,[1,size(manifoldIndex,2)]),@kurtosis);
    %if verbose,
    %    imagesc(VEL_HISTOGRAM_BOUNDARIES,VEL_HISTOGRAM_BOUNDARIES,mz{dind}')
    %end
end

% END MAIN -----------------------------------------------------------------------------------------