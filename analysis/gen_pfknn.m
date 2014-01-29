% pfrt = MTAAknnpfs(Trial,[15:30],'rear&theta',1);
% pfht = MTAAknnpfs(Trial,[15:30],'hwalk&theta',1);
% pflt = MTAAknnpfs(Trial,[15:30],'lwalk&theta',1);
% pfwt = MTAAknnpfs(Trial,[15:30],'walk&theta',1);
% pfwt_ps = MTAAknnpfs(Trial,[15:30],'walk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5);
% pflt_ps = MTAAknnpfs(Trial,[15:30],'lwalk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5);
% pfht_ps = MTAAknnpfs(Trial,[15:30],'hwalk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5);
% pfrt_ps = MTAAknnpfs(Trial,[15:30],'rear&theta',1,'numIter',1000,'ufrShufBlockSize',0.5);
% 
% pfrt = MTAAknnpfs(Trial,[15:30],'rear&theta',1,'binDims',[20,20],'distThreshold',70);
% pfht = MTAAknnpfs(Trial,[15:30],'hwalk&theta',1,'binDims',[20,20],'distThreshold',70);
% pflt = MTAAknnpfs(Trial,[15:30],'lwalk&theta',1,'binDims',[20,20],'distThreshold',70);
% pfwt = MTAAknnpfs(Trial,[15:30],'walk&theta',1,'binDims',[20,20],'distThreshold',70);
% pfwt_ps = MTAAknnpfs(Trial,[15:30],'walk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
% pflt_ps = MTAAknnpfs(Trial,[15:30],'lwalk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
% pfht_ps = MTAAknnpfs(Trial,[15:30],'hwalk&theta',1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
% pfrt_ps = MTAAknnpfs(Trial,[15:30],'rear&theta',1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
% 

Trial = MTATrial('jg05-20120317');
states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
pfns = cell(1,numel(states));
parfor s = 1:numel(states)
pfns{s} = MTAAknnpfs(Trial,[],states{s},1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end
