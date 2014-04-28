

sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};

pftype = 'MTAAknnpf';



numsts = numel(states);
f = figure;
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.xyz.load(Trial);
    Trial.load('nq');


for unit = pfs{ses,1}.data.clu
subplot(121),pfs{ses,1}.plot(unit);subplot(122),pfst{ses,1}.plot(unit);
saveas(f,['/gpfs01/sirota/bach/homes/gravio/figures/testing/' ...
          Trial.filebase '.test_knnpfs_batchmatch-' num2str(unit) '.png']);
end
end

for p = 1:5,

sum(pfs{2,p}.data.rateMap(:,unit,1)==pfst{2,p}.data.rateMap(:,unitN))

units = pfst{1,1}.data.clu;
for unit = units,
for unitN = units,
t(unit==units,unitN==units) = sum(pfs{ses,p}.data.rateMap(:,unit==pfs{ses,p}.data.clu,1)~=pfst{ses,p}.data.rateMap(:,unitN==pfst{ses,p}.data.clu)&~isnan(pfst{ses,p}.data.rateMap(:,unitN==pfst{ses,p}.data.clu)));
end
end



pfs = MTAAknnpfs(Trial,1:5,'theta',0,'numIter',1001, ...
                       'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);


pfstuff(u==units,state) = PlaceFieldStats(Trial,pfs{ses,state},u);    
