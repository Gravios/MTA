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

Dpath.xyz = '/storage2/ziyan/data/processed/xyz/ZY0519/'
Dpath.nlx = '/storage2/ziyan/data/processed/nlx/ZY0519/'
Dpath.spect = '/storage2/ziyan/data/processed/spect/ZY0519/'

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

frame = readFrame(vs);
vsSyncData = LoadBinary(fullfile(session.spath,[session.name,'.dat']),40,par.nChannels)';
vsSyncPulses = ThreshCross(vsSyncData,2000,0);
vsTimeStamps = mean(vsSyncPulses,2)./session.sampleRate;
vsFrameCount = vs.Duration*vs.FrameRate;
% DO they want the video in a matfile?

% Specify that reading should start at 0.5 seconds from the
% beginning.


% Create an axes
hfig = figure();
sax = axes();




frame = zeros([vsFrameCount,vs.Height,vs.Width]);
vs.CurrentTime = 0;
% Read video frames until available
for f = 1:vsFrameCount,
    tf = readFrame(vs);
    frame(f,:,:) = tf(:,:,1);
end


% CREATE MTAData object for spectral info
% IS this the final format for the processed objects

dsBld = load('/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/spect_446b.mat');
dsAch = load('/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/ZY0519-20190820_ACh_spect.mat');

fbr = load('/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/ZY0519-20190820.fbr.hpc.SENR.mat');

/storage2/ziyan/data/raw/spect/ZY0519/ZY0519-20190831/



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


% $$$ % NEW                    OLD
% $$$ sigHbClean               sigHbclean: [10995×845 double]
% $$$ hbFit                         Hbfit: [10995×3 double]
% $$$ sigNorm                        sign: [10995×845 double]
% $$$ t                             tnorm: [1×10995 double]
% $$$ spectACh                  ACh_spect: [10995×845 double]
% $$$ signal                          ACh: [10995×1 double]
% $$$ sampleRate               SampleRate: 10
% $$$ intergrationTime   IntergrationTime: 0.09
% $$$ % NEW fields
% $$$                            location: 'hpc'
% $$$                          signalType:'NULL'




Dpaths = struct('xyz',  '/storage2/ziyan/data/processed/xyz/ZY0519/',         ...
                'ephys','/storage2/ziyan/data/processed/nlx/ZY0519/',         ...
                'spect','/storage2/ziyan/data/processed/spect/ZY0519/',       ...
                'video','/storage2/ziyan/data/processed/video/ZY0519/');

MTAstartup('RSZY2019');

sessionList = get_session_list('ZY0519');


link_session(sessionList.sessionName,sessionList.dPaths);

! ln -s ../xyz/ZY0519-20190820/tfr/ZY0519-20190820-Trial001.c3d.mat /storage/gravio/data/project/RSZY2019/ZY0519-20190820/tfr/ZY0519-20190820-Trial001.c3d.mat
cd('/storage/gravio/data/project/RSZY2019/ZY0519-20190820/tfr/');
system(['find . -type l -delete'])
cd('/storage/gravio/data/project/RSZY2019/xyz/ZY0519-20190820');
system(['find . -type l -exec ln -s ../../xyz/ZY0519-20190820/{} ',fullfile(sessionPath,'{}'),' \;']);
ls -lh /storage/gravio/data/project/RSZY2019/ZY0519-20190820/tfr/
system(['find ',fullfile('..',field,sessionName),' -type l -exec  ln -s {} ',fullfile(sessionPath,'{}'),' \;']);
system(['find . -type d -exec  \;'])
system(['find ./tfr -type l -exec {} \;'])

cd('/storage/gravio/data/project/RSZY2019/ZY0519-20190820/tfr/');
system(['find . -type l -delete'])
cd('/storage/gravio/data/project/RSZY2019/xyz/ZY0519-20190820');
ls -lh /storage/gravio/data/project/RSZY2019/ZY0519-20190820/tfr/

system(['find . -type l -exec '...
        'file={}\; '...
        'file=echo "$pid"|awk ''{print $2}''\; '...
        'ln -s ../../xyz/ZY0519-20190820/$file ',fullfile(sessionPath,'$file'),' \;']);



s = MTASession(sessionList.sessionName,...
               sessionList.mazeName,...
               false,...
               sessionList.TTLValue,...
               sessionList.dLoggers,...
               sessionList.xyzSampleRate);


% Make the fbr load and sync
list_files(s.name,'fbr')


s.fbr = MTADfbr([]);
ACh = s.load('fbr','hpc','SENR','signal');

xyz = s.load('xyz');
xyz.resample(ACh);
vxy = xyz.vel(1,[1,2]);
rhm = fet_rhm(s);

for f = list_files(s.name,'fbr')
    fbr = load(fullfile(s.spath(),f{1}));
    fbr.integrationTime = fbr.intergrationTime;
    fbr = rmfield(fbr,'intergrationTime');
    system(['rm ' fullfile(s.spath(),f{1})]);    
    save(fullfile(s.spath(),f{1}),'-struct','fbr');
end


figure,plot(log10(vxy.data),ACh.data,'.');

fbr = {};
fbrFileNames = list_files(s.name,'fbr');
for f = fbrFileNames,
    fbr{end+1} = load(fullfile(s.spath(),f{1}));
end



Trial = MTATrial(s.name,'tfr','fs',true,s.sync.sync);
Trial = MTATrial(s.name,'tfr','fs',true,[1,2030;2323,4004])


xyz = Trial.load('xyz');
Trial.fbr = MTADfbr([]);
ACh = Trial.load('fbr','hpc','SENR','signal');
lfp = Trial.load('lfp',22);

lfp = s.load('lfp',22);

[ys,fs,ts] = fet_spec(s,                              ... Trial
                      lfp,                            ... lfp 
                      'mtchglong',                    ... mode
                      true,                           ... whiten signal
                      [],                             ... sampleRate (ignore)
                      struct('nFFT'     , 2^11, 'Fs', s.lfp.sampleRate,... spectrum parameters
                             'WinLength', 2^10, 'nOverlap',2^10*.875,   ...
                             'FreqRange',[1,40]),                     ...
                      'overwrite', true);

                             
%ACh = s.load('fbr','hpc','SENR','signal');
%xyz = s.load('xyz');
%xyz.resample(ACh);
%vxy = xyz.vel(1,[1,2]);
ys.resample(ACh);

figure,plot(ys(:,10),ACh.data,'.')

figure,plot(log10(ys(:,15)),vxy.data,'.')

figure,plot(log10(rys(:,6)),ACh.data,'.')

