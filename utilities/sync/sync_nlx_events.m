function evts = sync_nlx_events(Trial,eventLabels)

Trial = MTATrial.validate(Trial);

evts = LoadEvents(fullfile(Trial.spath, [Trial.name '.all.evt']));

eClu = find(~cellfun(@isempty,regexp(evts.Labels,eventLabels)));

eTime = evts.time(evts.Clu==eClu);

