
% ZY1219-20191107
 
sessionList = get_session_list('ZY1219');

sid = 1;
link_session(sessionList(sid).sessionName,sessionList(sid).dPaths);


filename = '/storage2/ziyan/data/raw/video/ZY1219/ZY1219-20191107/ZY1219-20191107.mp4'; 

VideoReader(filename);
AVPlay(filename);


build_sessions(sessionList);

Session  = MTASession.validate(sessionList);
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));            
lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),36,Par.nChannels,4)';
offset = find(lfp>2500,1,'first')./1250;