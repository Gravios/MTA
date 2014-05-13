function gen_unit_profile(Trial)

display=true;
stc_mode = 'auto_wbhr';
states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
numsts = numel(states);
Trial.stc.updateMode(stc_mode);Trial.stc.load;
Trial.load('nq');
Trial.load('xyz');

units = select_units(Trial,18,'pyr');



pfs = {};
for i=1:numsts,
pfs{i} = MTAAknnpfs(Trial,units,states{i},0,'numIter',1, ...
                     'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
end

if display
figure,
for u = pfs{1}.data.clu,
clims = [0, max(cellfun(@maxRate,pfs,repmat({u},1,5)))];
s = 1;
    for s = 1:numsts,
        subplotfit(s,numsts);pfs{s}.plot(u,[],[],clims);
        title([pfs{s}.parameters.states ' ' num2str(u)]);
    end
waitforbuttonpress
end
end
