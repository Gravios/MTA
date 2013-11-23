function bhv_psd(Trial,varargin)
[trialName,mazeName] = DefaultArgs(varargin,{'all'});
if isa(Trial,'MTASession')&~isa(Trial,'MTATrial')
    Trial = MTATrial(Trial,{},trialName,[],[],mazeName);
else
    Trial = MTATrial(Trial,{},trialName,[],[],mazeName);
end



h = spectrum.mtm;

hpsd = psd(h,