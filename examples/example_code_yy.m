%% run in shell
% mkdir /storage/weiwei/data/project
% mkdir /storage/weiwei/data/project/general
% mkdir /storage/weiwei/data/project/general/jg05-20120310
% cp /storage/gravio/data/project/general/jg05-20120310/jg05-20120310.*{ses,pos,trl,stc}*.mat /storage/weiwei/data/project/general/jg05-20120310/
% bash ./storage/share/matlab/MTA/configure

%% run in MATLAB
MTAstartup('general');

xyz_path = '/storage/gravio/data/processed/xyz/jg05';
nlx_path = '/storage/gravio/data/processed/nlx/';
states = {'theta','rear','walk','turn','pause','groom','sit'};
statesAUX = {'loc','hloc','lloc','hpause','lpause'};


% Link data from specified path to the MTA root directory of
% the current user.
if ~isempty(xyz_path) && ~isempty(nlx_path)
    linkSession('jg05-20120310',xyz_path,nlx_path);
end



%% Loading stuff

% This loads the trial
Trial = MTATrial.validate('jg05-20120310.cof.all');

% This loads the trajectory positions -> xyz.data [Time,Marker,cartDim];
xyz = Trial.load('xyz');

% list markers in the order observed in xyz.data's second dim
xyz.model.ml()

% You can filter data like this ...
xyz.filter('ButFilter',3,50,'low'); %see help gtwin for gaussian kernals

% get smoothed inter-marker angels -> ang.data [Time, Marker1, Marker2, sphDim]
% Marker segment angles are relative to the room coordinate system of Marker1.
ang = create(MTADang,Trial,filter(copy(xyz),'ButFilter',3,50,'low'));

figure,
    hist(ang(Trial.stc{'r'},'head_back','head_front',2),100)

%% Loading LFP

% for jg05 how to load the raw lfp of the hippocampus H64BUZ
lfp = Trial.load('lfp',1:64);

% for jg05 how to load the raw lfp of the hippocampus H32LIN
lfp = Trial.load('lfp' ,65:96);



%% spectrograms made easy
specArgs = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,20]);

[ys,fs,ts] = fet_spec(Trial,lfp,'mtcsdglong',true,[],specArgs,true,false);

figure();
    imagesc(ts,fs,log10(ys)');
    xlabel('Time');
    ylabel('Frequency');

    
    
%% Loading spikes 
spk = Trial.spk.copy();
spk.create(Trial,xyz.sampleRate,states{s},[],'deburst');
% or
spk.create(Trial,lfp.sampleRate,states{s},[],'deburst');



%% Behavior Segmentation stuff
% most behaviors are stored as periods in the stc field of a Trial
% Stc should be a MTAStateCollection object.

Stc = Trial.load('stc','msnn_ppsvd_raux');

% you can create a state matrix which fits your other data, such as xyz
stcm = stc2mat(Stc,xyz,states);
% or lfp
stcm = stc2mat(Stc,xyz,states);













