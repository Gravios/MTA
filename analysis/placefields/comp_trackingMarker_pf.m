function comp_trackingMarker_pf(Session)
% function comp_trackingMarker_pf(Session)
% Comparison of MTAApfs.m placefield while iterating through
% trackingMarkers

markerList = Session.Model.ml;

Pfs = {};
for marker = 1:length(markerList),
   Session.trackingMarker = markerList{marker};
   Pfs{marker} = MTAPlaceField(Session,[],{'head','theta'},1);
end

pfview(Pfs)

mrp = zeros(9,124);
for marker = 1:length(markerList),
    for unit = 1:124,
        tmr = Pfs{marker}.maxRate{unit}(Pfs{marker}.maxRateMax{unit});
        if isempty(tmr),continue,end
        mpr(marker,unit) = tmr;
    end
end

[~,mmr] = max(mpr);

mmr(find(mmr<5))=0;
mmr(find(mmr>4))=1;

