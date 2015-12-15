function feature = normalize_to_reference(feature,referenceFeature)
% function feature = normalize_features_to_reference_trial(Trial,feature,referenceTrial)
% This Function is only meant for use with fet_tsne at the moment.

NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-1,2,NBINS);

% Match current feature sampleRate
referenceFeature.resample(feature.sampleRate);


for f = [1:4,11:14];
    % oh god change this... change it now!!
    % Future MTAData objects should have model which has a isCirc
    % property to determin how to handle the vars
    if f<=5,
        minus_fun = @minus;
        median_fun = @median;
        std_fun = @std;
        stdThreshold = 10;
    else
        minus_fun = @circ_dist;
        median_fun = @circ_median;
        std_fun = @circ_std; 
        stdThreshold = .2;
    end

    
    % Once for the target Trial
    lind = feature.data(:,12)<0;    
    [~,ind_v_b] = histc(log10(feature.data(:,6)),VEL_HISTOGRAM_BOUNDARIES);
    [~,ind_v_h] = histc(log10(feature.data(:,8)),VEL_HISTOGRAM_BOUNDARIES);
    mind = nniz([ind_v_b,ind_v_h]);
    % Mean value       _____________________________________
    %                  |body speed index   |head speed index
    mz{1} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       feature.data(mind&lind,f),...
                       [NBINS,NBINS],...
                       median_fun);
    % Standard dev     _____________________________________
    %                  |body speed index   |head speed index
    sz{1} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],... 
                       feature.data(mind&lind,f),...
                       [NBINS,NBINS],...                           
                       std_fun);

    % And once for the reference Trial
    lind = referenceFeature.data(:,12)<0;    
    [~,ind_v_b] = histc(log10(referenceFeature.data(:,6)),VEL_HISTOGRAM_BOUNDARIES);
    [~,ind_v_h] = histc(log10(referenceFeature.data(:,8)),VEL_HISTOGRAM_BOUNDARIES);
    mind = nniz([ind_v_b,ind_v_h]);
    % Mean value       _____________________________________
    %                  |body speed index   |head speed index
    mz{2} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],...
                       referenceFeature.data(mind&lind,f),...
                       [NBINS,NBINS],...
                       median_fun);
    % Standard dev     _____________________________________
    %                  |body speed index   |head speed index
    sz{2} = accumarray([ind_v_b(mind&lind),ind_v_h(mind&lind)],... 
                       referenceFeature.data(mind&lind,f),...
                       [NBINS,NBINS],...                           
                       std_fun);

    
    nnz = mz{2}~=0&mz{1}~=0;
    mzd = minus_fun(mz{1},mz{2});
    szd = sz{2}(:)<10&sz{1}(:)<stdThreshold;
    feature.data(:,f) = minus_fun(feature.data(:,f),median(mzd(nnz(:)&szd)));

end


