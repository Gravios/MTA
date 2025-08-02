function [distance] = dist(Data,markerA,markerB)


% MAIN ---------------------------------------------------------------------------------------------
if ~isnumber(markerA)
    markerA = Data.model.gmi(markerA);
end
if ~isnumber(markerB)
    markerB = Data.model.gmi(markerB);
end

distance = sqrt(sum(diff(Data.data(:,[markerA,markerB],:),1,2).^2,3));

% END MAIN -----------------------------------------------------------------------------------------



