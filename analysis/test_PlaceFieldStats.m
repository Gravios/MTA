
Trial = MTATrial('jg05-20120317','all');
Trial.load('nq');
units = find(Trial.nq.SpkWidthR>0.8&Trial.nq.eDist>18)';
pfs =     MTAAknnpfs(Trial,units,'walk&theta',0,'numIter',1000, ...
                'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);


if matlabpool('size')~=12,matlabpool('open',12);end
clear ts,clear tss
parfor i = 1:numel(units),
[ts(i),tss(i)] = PlaceFieldStats(Trial,pfs,units(i));
end
%if matlabpool('size')==12,matlabpool('close');end

pfstats = CatStruct(ts,[],2);
pfshuff = CatStruct(tss,[],2);