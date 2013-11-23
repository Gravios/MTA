function compare_bhv_pfs(Session)

Session = MTASession('jg05-20120315');

trial_set =    {'crt1',     'alt1',     'crt2',     'alt2',     'crt3'};
ts_bhv_names = {'m20130304','m20130305','m20130307','m20130307','m20130307'};
state_set = {{'walk'},{'walk','theta'}};


pfs = {};
for t = 1:length(trial_set),
Trial = MTATrial(Session,[],trial_set{t});
for s = 1:length(state_set),
pfs{1,t,s} = MTAPlaceField(Trial,[],state_set{s},1)
end
Trial.Bhv = MTABhv(Trial,ts_bhv_names{t});
for s = 1:length(state_set),
pfs{2,t,s} = MTAPlaceField(Trial,[],state_set{s},1)
end
end


f = [];
f(1) = figure;
f(2) = figure;

for unit = 1:length(pfs{1}.cluMap),
for a = 1:2,
set(0,'CurrentFigure',f(a));
for t = 1:length(trial_set),
for s = 1:length(state_set),
subplot(2,5,t+length(trial_set)*(s-1));
%subplot(5,2,t+5*(s-1));
pfs{a,t,s}.plot(unit)
end
end
title(num2str(unit));
end
waitforbuttonpress
end
