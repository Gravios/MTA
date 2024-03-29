function [mz,vz,cz,czDom] = mean_embeded_feature_vbvhza(Data,Trial,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
[dims,featureDomainBoundaries,verbose] = DefaultArgs(varargin,{[],[],false},true);

NBINS = 20;
VEL_HISTOGRAM_BOUNDARIES = linspace(-3,2,NBINS);
HEIGHT_HISTOGRAM_BOUNDARIES = linspace(0,200,NBINS);


% INITIALIZE output
mz = {};
vz = {};
cz = {};
czDom = {};

%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

try,
    xyz = preproc_xyz(Trial,'trb');
catch err
    disp(err)
    xyz = preproc_xyz(Trial);
end

xyz.resample(Data);

vxy = xyz.copy;
vxy.filter('ButFilter',3,2.4,'low');
vxy = vxy.vel({'bcom','hcom'},[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

% FIND the bin indicies for each marginal distribution
[~,ind_v_b] = histc(vxy(:,'bcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_v_h] = histc(vxy(:,'hcom'),VEL_HISTOGRAM_BOUNDARIES);
[~,ind_z_a] = histc(xyz(:,'acom',3),HEIGHT_HISTOGRAM_BOUNDARIES);

% SELECT only non-zero points whithout rearing and grooming
if ~isempty(regexpi(Trial.stc.mode,'^hand_labeled.*')),
    ind = Trial.stc('a-m-r').cast('TimeSeries');
    ind.resample(Data);
elseif ~isempty(regexpi(Trial.stc.mode,'^mswnn_ppsvd$')),    
    ind = Trial.stc('a-m-r').cast('TimeSeries');
    ind.resample(Data);
else,
    ang = create(MTADang,Trial,xyz);    
    ind = Trial.stc('a').cast('TimeSeries');
    ind.resample(Data);
    ind.data(isnan(ind.data))=0;
% REQUIRE a hard thresholds to remove groom and rear like periods
    bang = ang(:,'spine_lower','hcom',3).*cos(ang(:,'spine_lower','hcom',2));
    %bang = ang(:,'spine_lower','spine_upper',3).*cos(ang(:,'spine_lower','spine_upper',2));    
    bang = MTADxyz('data',bang/prctile(bang(nniz(bang)),95),'sampleRate',xyz.sampleRate);
    bang.data(~nniz(bang.data))=0;
    ind.data = ind.data&bang>.8&bang<1.1;
end

mind = nniz([ind_v_b,ind_v_h,ind_z_a])&ind.data==1;
manifoldIndex = [ind_v_b(mind),ind_v_h(mind),ind_z_a(mind)];

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
    czDom{dind} = accumarray(manifoldIndex,ones(size(acvar)),repmat(NBINS,[1,size(manifoldIndex,2)]),@sum);
    %if verbose,
    %    imagesc(VEL_HISTOGRAM_BOUNDARIES,VEL_HISTOGRAM_BOUNDARIES,mz{dind}')
    %end
end

% END MAIN -----------------------------------------------------------------------------------------