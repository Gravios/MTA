% *** req20190829
%     Tags: synchronization spectrometer ach
%     Status: retired
%     Type: utility
%     Author: Justin Graboski
%     Final_Forms: 
%     Project: 
%     Description: working on -
%                  synchronization of video to openephys signal
%                  synchronizaiton of ACh signal, processed from spectrometer data


link_session('ZY0519-20190820',...
             '/storage2/ziyan/data/processed/xyz/ZY0519/',...
             '/storage2/ziyan/data/processed/nlx/ZY0519/');


session = MTASession.validate('ZY0519-20190820.tfr.all');

xyz = session.load('xyz');



% SYNCHRONIZE the video
% ZY0519-20190820.mp4  ZY0519-20190820_videotime.txt.txt

subject = 'ZY0519';
session = MTASession.validate('ZY0519-20190820.tfr.all');

%session.path.data
%getenv('MTA_ROOT_PATH')

videoFile = '/storage2/ziyan/data/raw/video/ZY0519/ZY0519-20190820/ZY0519-20190820.mp4';
vs = VideoReader(videoFile);
vsSyncChannel = 40;
par = LoadPar(fullfile(session.spath,[session.name,'.xml']));

vsSyncData = LoadBinary(fullfile(session.spath,[session.name,'.dat']),40,par.nChannels)';
vsSyncPulses = ThreshCross(vsSyncData,2000,0);
vsTimeStamps = mean(vsSyncPulses,2)./session.sampleRate;

% DO they want the video in a matfile?

% Specify that reading should start at 0.5 seconds from the
% beginning.
vs.CurrentTime = 1/vs.FrameRate;

% Create an axes
hfig = figure();
sax = axes();
        
% Read video frames until available
while hasFrame(vs),
    frame = readFrame(vs);
    image(frame, 'Parent', sax);
    sax.Visible = 'off';
    pause(1/vs.FrameRate);
end


% CREATE MTAData object for spectral info
% IS this the final format for the processed objects

dsBld = load('/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/spect_446b.mat');
dsAch = load('/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/ZY0519-20190820_ACh_spect.mat');

% dsAch.AChAll:  { []    [1×1 struct] }
% dsAch.AChAll{2}
%           sigHbclean: [10995×845 double]
%                Hbfit: [10995×3 double]
%                 sign: [10995×845 double]
%                tnorm: [1×10995 double]
%            ACh_spect: [10995×845 double]
%                  ACh: [10995×1 double]
%           SampleRate: 10
%     IntergrationTime: 0.09
%    
% dsBld.spect_446b: {[]  [10851×845 double]}
%