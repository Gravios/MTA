function phz = load_thetarc_phase(Trial,varargin)
%function phz = load_theta_phase(Trial, varargin)
%
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(...
    'samplerate', 250, ...
    'channels', [] ...
);
[samplerate, channels] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

Trial.lfp.filename = [Trial.name,'.lfp'];        

if ~isempty(Trial.meta) & ~isempty(channels);
    lfp = Trial.load('lfp',Trial.meta.channelGroup.thetarc);
else
    lfp = Trial.load('lfp', channels);
end

lfp.data = diff(lfp.data,1,2);

phz = lfp.phase([5,13]);
phz.data = unwrap(phz.data);
phz.resample(samplerate);
if ~isempty(Trial.meta)
    phz.data = phz.data + Trial.meta.correction.thetaPhase;
end

phz.data = mod(phz.data + 2*pi, 2*pi);
