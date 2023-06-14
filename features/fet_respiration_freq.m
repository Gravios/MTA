function [fet] = fet_respiration_freq(Trial,varargin)
% [ifreq] = fet_respiration_freq(Trial,varargin)
% 
% Computes a rostrocaudal motion feature of the head.
% 
% Varargin:
%           NAME : TYPE    : DEFAULT VAL          : Description
%     sampleRate : Numeric : Trial.xyz.sampleRate : computational sampling rate
%            tag : String  : ''                   : save file tag 
%      overwrite : Logical : false                : flag for overwriting saved data
%

% DEFARGS ------------------------------------------------------------------------------------------
Trial = MTATrial.validate(Trial);
defargs = struct('sampleRate',   250,                                                            ...
                        'tag',   '',                                                             ...
                  'overwrite',   false                                                           ...
);
[sampleRate,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash(struct('function',mfilename()));
end
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [Trial.filebase,'.fet_respiration_freq.',tag,'.mat'],...
              [],...
              sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'instantaneous respiration frequencey','fet_respiration_freq','f');

if overwrite || ~exist(fet)

% DEFVAR position
    xyz = preproc_xyz(Trial,'trb',sampleRate);

    [amins,ncpSampleRate] = find_respiration_troughs(Trial,Trial.meta.channelGroup.respiration);

% COMPUTE the instantaneous respiratory frequency in the target sampling rate
    rmins = round((amins-1)/ncpSampleRate.*xyz.sampleRate);
    fet.data(1:amins(1),1) = 0;
    for mm = 1:numel(rmins)-1,
        fet.data(rmins(mm):rmins(mm+1),1) = 1./((amins(mm+1)-amins(mm))./ncpSampleRate);
    end
    fet.data(amins(mm+1):end,1) = 0;
    
% SAVE respiration frequency feature to file
    fet.save();
else
% LOAD respiration frequency feature from file    
    fet.load();
end
