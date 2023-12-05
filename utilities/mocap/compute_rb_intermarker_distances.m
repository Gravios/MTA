function [distances] = compute_rb_intermarker_distances(rb)

% MAIN ---------------------------------------------------------------------------------------------

pairCount = 1;
distances = zeros( [ size(rb,1), (rb.model.N^2-rb.model.N)/2] );

for markerA = 1 : rb.model.N
    for markerB = markerA+1 : rb.model.N
        %distances(:,pairCount) = rb.dist(markerA,markerB);
        distances(:,pairCount) = sqrt(sum(diff(rb(:,[markerA,markerB],:),1,2).^2,3));
        pairCount = pairCount+1;
    end
end

% END MAIN -----------------------------------------------------------------------------------------


