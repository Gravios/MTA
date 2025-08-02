function Session = sync_vicon_openephys(Session,TTLValue,viconSampleRate)

Par = LoadPar('/storage/javier/Raw_data/ephys/JM11/2016-11-25_18-34-06/processed/2016-11-25_18-34-06.xml');;
lfp = LoadBinary('/storage/javier/Raw_data/ephys/JM11/2016-11-25_18-34-06/processed/2016-11-25_18-34-06.lfp',...
                 145:164,Par.nChannels,4)';    

Session.spath = '/storage/javier/data/processed/xyz/JM11-20161126/';
Session.maze.name = 'cof';
Session.name = 'JM11-20161126';

[xyzData, markers, viconSampleRate] = concatenate_vicon_files(Session);            

xyzSampleRate = 199.997752;
(diff(ThreshCross(lfp(:,13),1e4,1),1,2)-4)./1250-cellfun(@length,xyzData(1:5))'./xyzSampleRate

Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
lfp = LoadBinary(fullfile(Session.spath, [Session.name '.lfp']),1,Par.nChannels,4)';    

