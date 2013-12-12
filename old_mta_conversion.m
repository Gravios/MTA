function old_mta_conversion(session_name)

path = '/data/homes/gravio/data/analysis/';

load(fullfile(path,session_name,[session_name 'cof.all.xyz.mat']));

save(fullfile(Session.spath.analysis,[Session.filebase '.pos.mat']),'data');

