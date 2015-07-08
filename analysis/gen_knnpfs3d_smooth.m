%MTAstartup('cin','cin');
%Trial = MTATrial('Ed10-20140812');
%MTAstartup;
Trial = MTATrial('er06-20130614','fly');
Trial = MTATrial('er06-20130614','all-cof');
%Session = MTASession('er06-20130614');
%Trial = MTATrial('jg05-20120317');
%Trial = MTATrial('jg05-20120310');
%Trial = MTATrial('jg05-20120309');
display = true;

%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'flight'};
states = {'u_wr','flight'};

wrper = Trial.stc{'r'}.copy;
wrper.data = [Trial.stc{'r'}.data(1:2:end,:);Trial.stc{'w'}.data(1:2:end,:)];
wrper.data = sort(wrper.data);
wrper.key = 'x';
wrper.label = 'u_wr';
wrper.updateFilename([Trial.filebase '.sst.' wrper.label '.' wrper.key '.mat']);
Trial.stc.states{end+1} = wrper;

'er06-20130614.cof.all-cof.sst.rear.r.mat'

states = {'u_wr'};
%states = {'hswalk&theta','lswalk&theta'};
nsts = size(states,2);


overwrite = true;
units = select_units(Trial,18,'pyr');
units = 1:100;


for i = 1:nsts,
MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20,20],'distThreshold',125,'nNearestNeighbors',150,'type','xyz');
end

for i = 1:nsts,
MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1,'ufrShufBlockSize',0,'binDims',[40,40,40],'distThreshold',125,'nNearestNeighbors',150,'type','xyz');
end


pfs = {};
for i = 1:nsts,
pfs{i} = MTAApfs(Trial,units,states{1},overwrite,'numIter',1,'binDims',[40,40,40],'type','xyz','SmoothingWeights',[1.8,1.8,1.8]);
end


if display,
    pfs ={};
    for i = 1:nsts,
        pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20,20],'distThreshold',125,'nNearestNeighbors',150,'type','xyz');
    end

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        sp = [];
% $$$         mrate = zeros([1,nsts]);
% $$$         for i = 1:nsts,
% $$$             mrate(i) = max(pfs{i}.data.rateMap(:,pfs{i}.data.clu==unit));
% $$$         end
for i=1:3

 
subplot2(3,2,i,1), pfs{i}.plot(unit,'xz',1); 
title([pfs{i}.session.trialName,':',pfs{i}.parameters.states,': ',num2str(unit)]);
subplot2(3,2,i,2), pfs{i}.plot(unit,'xy',1); 
end
% $$$         for i = 1:nsts,
% $$$             sp(i) = subplotfit(i,nsts);cla
% $$$             pfs{i}.plot(unit,'isosurface')%,true,[0,max(mrate)]);
% $$$             title([pfs{i}.parameters.states,': ',num2str(unit)]);
% $$$         end
% $$$         linkprop(sp, {'CameraPosition'}); 
        unit = figure_controls(hfig,unit,units);
    end

end
