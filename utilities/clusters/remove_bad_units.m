function units = remove_bad_units(Trial,units)
ds = load(fullfile(Trial.spath,[Trial.name,'.bad_units.mat']));
units(ismember(units,ds.badUnits)) = [];