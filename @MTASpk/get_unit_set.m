function units = get_unit_set(Spk,Session,name)
units = load(fullfile(Session.spath,[Session.filebase,'.unit_set.',name,'.mat']));