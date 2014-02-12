MTAConfiguration('/gpfs01/sirota/bach/data/gravio/','absolute');
Trial = MTATrial('jg05-20120310');
Trial.xyz.load(Trial);
Trial.load('nq');

units = find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19);
states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
numsts = numel(states);

pfs={};
for i = 1:numsts,
    pfs{i}  =  MTAAknnpfs(Trial,units,states{i},1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end
 