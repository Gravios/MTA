function [shift] = align_depth_probe_csd(Trial,meta,lfp,varargin)
% function [shift] = align_depth_probe_csd(Trial,lfp,varargin)
%  
% Temporary function primarly for the alignment of a linear probe in subjet jg05
% 
% Return the depth shift required to algin the csd profiles of a target and a reference trial 
% compute theta phase triggered csd average 
%
% NOTE : only good for shifts less than 150um.

% DEFARGS ------------------------------------------------------------------------------------------

defargs = struct(   'RefTrial' , 'jg05-20120310.cof.all',                                       ...
                 'refChannels' , [65:96],                                                       ...
             'refThetaChannel' , 65,                                                            ...
               'corThetaPhase' , 0.785398163397448,                                             ...  
             'channelInterval' , 1,                                                             ...
                'channelPitch' , 50,                                                            ...
                   'probeName' , 'LIN32');
[RefTrial,refChannels,refThetaChannel,corThetaPhase,channelInterval,channelPitch,probeName] =   ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% $$$ % TESTVARS 
% $$$ configure_default_args();
% $$$ MjgER2016_load_data();
% $$$ Trial = Trials{17};
% $$$ meta = sessionList(17);
% $$$ lfp = Trial.load('lfp',[65:96]);
% $$$ 
% $$$ RefTrial = 'jg05-20120310.cof.all';
% $$$ refChannels = [65:96];
% $$$ refThetaChannel = 65;
% $$$ corThetaPhase = 0.785398163397448;
% $$$ channelInterval = 1;
% $$$ channelPitch =  50;
% $$$ probeName = 'LIN32';

%Trial.subject.ephys.channelGroup.theta;

RefTrial = MTATrial.validate(RefTrial);
sampleRate = 250;
thetaBins = 32;
numEdgeInds = 4;

%pathFileCsd = fullfile(RefTrial.spath,[RefTrial.filebase,'.',mfilename(),'-',probeName,'.mat']);
pathFileCsd = fullfile(RefTrial.spath,...
                       [RefTrial.filebase,'.',mfilename(),'-',probeName,'.mat']);
% $$$ pathFileCsd = fullfile(RefTrial.spath,...
% $$$                        [RefTrial.filebase,'.','align_depth_probe_csd','-',probeName,'.mat']);
if ~exist(pathFileCsd,'file')
    
    thetaPhase = load_theta_phase(RefTrial, ...
                           sampleRate,      ...
                           refThetaChannel, ...
                           corThetaPhase);
    
    rcsd = csd(resample(RefTrial.load('lfp',refChannels), ...
                        thetaPhase),                      ...
               channelInterval,                           ...
               channelPitch);
    
    sind = logical(get(resample(cast([RefTrial.stc{'theta-groom-sit&gper'}],'TimeSeries'),rcsd),'data'));

    thetaPhaseEdgs = linspace(0,2*pi,thetaBins+1);
    thetaPhaseCntr = mean([thetaPhaseEdgs(1:end-1);thetaPhaseEdgs(2:end)]);
    thetaPhaseInds = discretize(thetaPhase.data,thetaPhaseEdgs);

    meanThetaCsdRef = zeros([size(rcsd,2),numel(thetaPhaseCntr)]);

% COMPUTE theta phase aligned average csd
    for pt = 1:numel(thetaPhaseCntr),
        ind = pt==thetaPhaseInds & sind;
        meanThetaCsdRef(:,pt) = mean(rcsd(ind,:));
    end
    
    save(pathFileCsd,'meanThetaCsdRef');
else
    load(pathFileCsd,'meanThetaCsdRef');
end


thetaPhase = load_theta_phase(Trial,                            ...
                              sampleRate,                       ...
                              meta.subject.channelGroup.theta,  ...
                              meta.subject.correction.thetaPhase);

tcsd = csd(resample(filter(copy(lfp),'ButFilter',4,[0.8,20],'bandpass'),                        ...
                    thetaPhase),                      ...
           channelInterval,                           ...
           channelPitch);

tind = logical(get(resample(cast([Trial.stc{'theta-groom-sit&gper'}],'TimeSeries'),tcsd),'data'));

thetaPhaseEdgs = linspace(0,2*pi,thetaBins+1);
thetaPhaseCntr = mean([thetaPhaseEdgs(1:end-1);thetaPhaseEdgs(2:end)]);
thetaPhaseInds = discretize(thetaPhase.data,thetaPhaseEdgs);

meanThetaCsdTrg = zeros([size(tcsd,2),numel(thetaPhaseCntr)]);

% COMPUTE theta phase aligned average csd
for pt = 1:numel(thetaPhaseCntr),
    ind = pt==thetaPhaseInds & tind;
    meanThetaCsdTrg(:,pt) = mean(tcsd(ind,:));
end

filtMeanThetaCsdTrg = imgaussfilt(meanThetaCsdTrg,1);
filtMeanThetaCsdRef = imgaussfilt(meanThetaCsdRef,1);

% FIND shift which minimizes the mean squared error
for shift = 1:2*numEdgeInds,
    mse(shift) = sqrt(sum(nonzeros(filtMeanThetaCsdTrg(shift:end-(2*numEdgeInds-shift),:)-filtMeanThetaCsdRef(numEdgeInds:end-numEdgeInds,:)).^2));
end
[~,shift] = min(mse);
shift = shift-numEdgeInds;
