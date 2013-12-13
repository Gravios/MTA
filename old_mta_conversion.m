function old_mta_conversion(session_name)

session_name = 'jg05-20120310';

path = '/data/homes/gravio/data/analysis/';
xyzfilename =[session_name '.cof.all.xyz.mat'];
xyzfilepath= fullfile(path,session_name,xyzfilename);
load(xyzfilepath);
data = Session.xyz;

save(fullfile(Session.spath.analysis,[Session.filebase '.pos.mat']),'data');

clear('Session');

new_path = ['/data/homes/gravio/data/' session_name];

Session = MTASession(session_name,[],1,'0x0040'); % add the extra parameters


dsync(dsync~=0)=1;
tsync = ThreshCross(dsync,0.5,1)+Session.xyz.origin;
tsync(1) = tsync(1)+800;
Trial = MTATrial(Session,'all','cof',true,tsync./Session.xyz.sampleRate);

figure
plot(data(:,7,3))
hold on,
plot(Session.xyz.data(:,7,3),'c')
plot(dsync*10,'m')
hold on,plot([1:Trial.xyz.size]+Trial.xyz.origin-Session.xyz.origin,Trial.xyz.data(:,7,3),'m')





%posfilename =[session_name '.cof.all.pos.mat'];
angfilename =[session_name '.cof.all.ang.mat'];

%system(['mv ' fullfile(path,posfilename) ' ' fullfile(new_path,posfilename)])
system(['mv ' fullfile(path,angfilename) ' ' fullfile(new_path,angfilename)])