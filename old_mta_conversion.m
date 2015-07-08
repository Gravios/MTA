function old_mta_conversion(session_name)

session_name = 'er06-20130614';

path = '/data/homes/gravio/data/analysis/';
xyzfilename =[session_name '.cof.all.xyz.mat'];
xyzfilepath= fullfile(path,session_name,xyzfilename);
load(xyzfilepath);
data = Session.xyz;

save(fullfile(Session.spath.analysis,[Session.filebase '.pos.mat']),'data');

clear('Session');

new_path = ['/data/homes/gravio/data/' session_name];

Session = MTASession(session_name,[],1,'0x8000',199.997752); % add the extra parameters

allcof = [1,120732;122051,242449;243980,364905;609966,730366;...
          731912,852910;853412,974678;1219663,1340812;1341248,1461592];

allcof = (allcof+Session.xyz.origin)./Session.xyz.sampleRate;

%for er06-20130614
Trial = MTATrial(Session,'all-cof','cof',true,allcof./Session.xyz.sampleRate);

dsync(dsync~=0)=1;
tsync = ThreshCross(dsync,0.5,1)+Session.xyz.origin;
tsync(1) = tsync(1)+800;
Trial = MTATrial(Session,'all-coff','cof',true,tsync./Session.xyz.sampleRate);

xsync = Session.xyz.sync.copy;
xsync.resample(1);
xsync.data = xsync.data([1,2,3,6,7,8,11,12],:)
Trial = MTATrial(Session,'all-cof','cof',true,xsync.data);

figure
plot(data(:,7,3))
plot(Session.xyz.data(:,7,3),'m')

figure,
hold on,
plot([1:Session.xyz.size(1)]+Session.xyz.origin,Session.xyz.data(:,7,3))
%plot(dsync*10,'m')
hold on,plot([1:Trial.xyz.size]+Trial.xyz.origin-Session.xyz.origin,Trial.xyz.data(:,7,3),'m')

hold on,plot([1:Trial.xyz.size]+Trial.xyz.origin,Trial.xyz.data(:,7,3),'m')



%posfilename =[session_name '.cof.all.pos.mat'];
angfilename =[session_name '.cof.all.ang.mat'];

%system(['mv ' fullfile(path,posfilename) ' ' fullfile(new_path,posfilename)])
system(['mv ' fullfile(path,angfilename) ' ' fullfile(new_path,angfilename)])