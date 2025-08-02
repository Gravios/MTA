function units = remove_bad_units(Trial,units)
units( ismember( units, Trial.spk.get_unit_set(Trial,'badcells'))) = [];