function comp_smoothness_pf(Session)

Session = MTASession('jg05-20120310');
Session = MTATrial(Session,'all');

Session.trackingMarker = 'head_front';

smooth = 0.01:0.005:.05;

Pfsmooth = {};
for sm = 1:length(smooth),
   Pfsmooth{sm} = MTAPlaceField(Session,[],{'head','theta'},1,[],[],[],smooth(sm));
end

mrs = zeros(9,124);
for sm = 1:length(smooth),
for unit = 1:124,
tmr = Pfsmooth{sm}.maxRate{unit}(Pfsmooth{sm}.maxRateMax{unit});
if isempty(tmr),continue,end
   mrs(sm,unit) = tmr;
end
end


p = [];
for unit = 1:124,
p(unit,:) = polyfit(smooth,unity(mrs(:,unit))',2);
end

pfview(Pfsmooth)





pf_search = MTAPlaceField([]);

pf_search.stateLabel = 'head.theta';
pf_search.trackingMarker = 'head_front';
Session = Session.load_Pfs(pf_search);
pfview(Session.Pfs)
