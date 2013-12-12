function old_mta_conversion(session_name)

path = '/data/homes/gravio/data/analysis/';
xyzfilename =[session_name 'cof.all.xyz.mat'];
xyzfilepath= fullfile(path,session_name,xyzfilename);
load(xyzfilepath);

save(fullfile(Session.spath.analysis,[Session.filebase '.pos.mat']),'data');

new_path = ['/data/homes/gravio/data/' session_name];

Session = MTASession(session_name,[],1); % add the extra parameters

posfilename =[session_name 'cof.all.pos.mat'];
angfilename =[session_name 'cof.all.ang.mat'];

system(['mv ' fullfile(path,posfilename) ' ' fullfile(new_path,posfilename)])
system(['mv ' fullfile(path,angfilename) ' ' fullfile(new_path,angfilename)])