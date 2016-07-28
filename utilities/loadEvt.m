function evts = loadEvt(Trial,eventLabels)
MTAstartup('vr_exp');
Trial = MTATrial.validate('Ed10-20140820.rov.all');
eventLabels = {'teleport'};

evts = LoadEvents(fullfile(Trial.spath, [Trial.name '.all.evt']));

eClu = find(~cellfun(@isempty,regexp(evts.Labels,eventLabels)));

eTime = evts.time(evts.Clu==eClu);

