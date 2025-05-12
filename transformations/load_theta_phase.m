function phz = load_theta_phase(Trial, varargin)
%function phz = load_theta_phase(Trial,samplerate)
%
%
% DEFARGS ---------------------------------------------------------------------
defargs = struct('samplerate', 250,                                         ...
                 'channel',    []                                           ...
);
[samplerate, channel] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------
Trial.lfp.filename = [Trial.name,'.lfp'];        
if ~isempty(Trial.meta) & isempty(channel)
    lfp = Trial.load('lfp', Trial.meta.channelGroup.theta);
else
    lfp = Trial.load('lfp', channel);
end

phz = lfp.phase([5,13]);
phz.data = unwrap(phz.data);
phz.resample(samplerate);
if ~isempty(Trial.meta)
    phz.data = phz.data + Trial.meta.correction.thetaPhase;
end    
phz.data = mod(phz.data + 2*pi, 2*pi);
