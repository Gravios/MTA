function [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,halfSpkWindow)
%function [posterior,unitCount] = compute_bayesian_posterior_ufr_boxcar(Trial,index,sampleRate,ufr,ratemaps,halfSpkWindow)
% COMPUTE Posterior distribution base on ratemaps and unit firing rates
logRatemap = log(ratemaps);
priorRatemap = -sum(ratemaps,2,'omitnan')*ufr.spikeWindow;

% SELECT ratemap bins within mask
twind = index+round(halfSpkWindow.*sampleRate).*[-1,1];
cufr = sum(ufr(twind(1):twind(2),:));

gind = cufr>0.5/sampleRate;

if any(gind)
    unitCount = sum(gind);            
    % COMPUTE posterior
    posterior = exp(priorRatemap+logRatemap*(cufr+eps)');
    posterior = posterior./sum(posterior,'omitnan');
else
    unitCount = [];
    posterior = [];
end
