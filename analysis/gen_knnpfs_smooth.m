%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
MTAstartup;
%Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');

states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
%states = {'hswalk&theta','lswalk&theta'};
nsts = size(states,2);


overwrite = true;
units = select_units(Trial,18,'pyr');

pfs ={};
for i = 1:size(states,2),
MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[30,30],'distThreshold',125,'nNearestNeighbors',110);
end
