function time = d2t(X,sfreq,offset)
%function time = d2t(X,sfreq,offset)
time = (0:1:(X-1))/sfreq+0.5*offset/sfreq;
