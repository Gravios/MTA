function set_unit_set(Spk,Session,name,units)
save(fullfile(Session.spath,[Session.filebase,'.unit_set.',name,'.mat']),'units');