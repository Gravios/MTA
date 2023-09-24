function phz = load_thetarc_phase(Trial,sampleRate)
%function phz = load_theta_phase(Trial,sampleRate)
%
%    
Trial.lfp.filename = [Trial.name,'.lfp'];        
try,
    lfp = Trial.load('lfp',Trial.meta.channelGroup.thetarc);
catch
    lfp = Trial.load('lfp',Trial.meta.channelGroup.thetarc);
end
lfp.data = diff(lfp.data,1,2);

phz = lfp.phase([5,13]);
phz.data = unwrap(phz.data);
phz.resample(sampleRate);
phz.data = mod(phz.data+2*pi+Trial.meta.correction.thetaPhase,...
               2*pi);
