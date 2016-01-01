function feature = normalize_to_reference(feature,referenceFeature,varargin)
% function feature = normalize_features_to_reference_trial(Trial,feature,referenceTrial)
% This Function is only meant for use with fet_tsne at the moment.

[normEachSyncEpoch] = DefaultArgs(varargin,{false},1);

NBINS = 100;
VEL_HISTOGRAM_BOUNDARIES = linspace(-1,2,NBINS);

% Match current feature sampleRate
referenceFeature.resample(feature.sampleRate);

if normEachSyncEpoch,
    inSync = feature.sync&feature.sync.sync;
    inSync.resample(feature);
    inSync = inSync.data'-inSync.data(1)+1;
    inSync(end) = feature.size(1);
else
    inSync = [1;feature.size(1)];
end

switch feature.label
  case 'fet_tsne',
    fets = [1:4];
    fref1 = 6;
    fref2 = 8;
  case 'fet_tsne_rev1',
    fets = [1:4];
    fref1 = 6;
    fref2 = 8;
  case 'fet_tsne_rev2',
    nf = 18;
    fref1 = 6+nf;
    fref2 = 8+nf;
    fets = [1:4]+nf;
  case 'fet_tsne_rev3',
    fets = [1:4];
    fref1 = 6;
    fref2 = 8;
  otherwise ,
    fets = [];
end  

for ind = inSync,
    for f = fets;%,11:14
                 % oh god change this... change it now!!
                 % Future MTAData objects should have model which has a isCirc
                 % property to determin how to handle the vars
% $$$     if f<=5,
        minus_fun = @minus;
        median_fun = @nanmedian;
        std_fun = @nanstd;
        stdThreshold = 10;
% $$$     else
% $$$         minus_fun = @circ_dist;
% $$$         median_fun = @circ_median;
% $$$         std_fun = @circ_std; 
% $$$         stdThreshold = .2;
% $$$     end

        
        % Once for the target Trial
        lind = feature.data(ind(1):ind(2),12)<0;    
        [~,ind_v_b] = histc(feature.data(ind(1):ind(2),fref1),VEL_HISTOGRAM_BOUNDARIES);
        [~,ind_v_h] = histc(feature.data(ind(1):ind(2),fref2),VEL_HISTOGRAM_BOUNDARIES);
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
        lind = referenceFeature.data(ind(1):ind(2),12)<0;    
        [~,ind_v_b] = histc(referenceFeature.data(ind(1):ind(2),6),VEL_HISTOGRAM_BOUNDARIES);
        [~,ind_v_h] = histc(referenceFeature.data(ind(1):ind(2),8),VEL_HISTOGRAM_BOUNDARIES);
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
        switch feature.label
          case 'fet_tsne_rev2',            
            f = [f-nf,f,f+nf];
        end
        feature.data(ind(1):ind(2),f) = bsxfun(minus_fun,...
                            feature.data(ind(1):ind(2),f),median(mzd(nnz(:)&szd)));

    end


end