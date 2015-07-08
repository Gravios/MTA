function out = exp_pfsOD(Trial)

if matlabpool('size')<12,matlabpool open 12,end
svar = [];states={};stateSize =[];velMean =[];velStd =[];
parfor i = 1:200,
    [svar(i),states{i},stateSize(i),velMean(i),velStd(i)] = pfs_overdispertion(Trial,'rnd1');
end
out.svar = svar;
out.states = states;
out.stateSize = stateSize;
out.velMean = velMean;
out.velStd = velStd;
matlabpool close force local