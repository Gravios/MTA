function comp_trackingMarker_pf(Session)

ml = Session.Model.ml;

Pfs = {};
for marker = 1:length(ml),
   Session.trackingMarker = ml{marker};
   Pfs{marker} = MTAPlaceField(Session,[],{'head','theta'},1);
end

pfview(Pfs)

mrp = zeros(9,124);
for marker = 1:length(ml),
for unit = 1:124,
tmr = Pfs{marker}.maxRate{unit}(Pfs{marker}.maxRateMax{unit});
if isempty(tmr),continue,end
   mpr(marker,unit) = tmr;
end
end

[~,mmr] = max(mpr);

mmr(find(mmr<5))=0;
mmr(find(mmr>4))=1;

