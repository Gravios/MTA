function phz = load_theta_phase(Trial,sampleRate,channel,phzCorrection)
%function phz = load_theta_phase(Trial,sampleRate,channel)
%
%    
Trial.lfp.filename = [Trial.name,'.lfp'];        
lfp = Trial.load('lfp',channel);
phz = lfp.phase([5,13]);    
phz.data = unwrap(phz.data);
phz.resample(sampleRate);
phz.data = mod(phz.data+2*pi+phzCorrection,2*pi);
