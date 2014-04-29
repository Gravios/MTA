MTAConfiguration('/gpfs01/sirota/bach/data/gravio/','absolute');
Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);
Trial.load('nq');
Trial.stc.updateMode('auto_wbhr');Trial.stc.load;

units = find(Trial.nq.SpkWidthR>0.3&Trial.nq.eDist>19);
states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
numsts = numel(states);


units = [18,21:32,36:39]';

pfsc={};
for i = 1:numsts,
    pfsc{i}  =  MTAAPfknncorm(Trial,units,states{i},1,'numIter',1,'ufrShufBlockSize',0,'binDims',[30,30],'distThreshold',200,'nNearestNeighbors',300);
end

pfsr={};
for i = 1:numsts,
    pfsr{i}  =  MTAAknnpfs(Trial,units,states{i},1,'numIter',1,'ufrShufBlockSize',0,'binDims',[30,30],'distThreshold',70);
end



un=38;
pfrt  =  MTAAknnpfs(Trial,un,'theta',0,'numIter',1, ...
                       'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
figure,
subplot(122),pfrt.plot(un,[],1);


pfst  =  MTAAPfknncorm(Trial,un,'theta',1,'numIter',1, ...
                       'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
%matlabpool open 11
parfor ind = 1:2500,
[rc(ind),pc(ind)] = corr(sq(distdw(1:nnn,:,ind)),(smywut(1:nnn,:,ind)),'Type','Spearman');
end
rc(repmat(pfknnmdw,[1,1,size(ufr,2)])>dthresh) = nan;
rc(pc>0.001) = nan;
%figure,
subplot(121),imagescnan(rot90(reshape(rc',50,50)'),[],[],1)


ind = 962,plot(sq(distdw(1:nnn,:,ind)),(smywut(1:nnn,:,ind)),'.')