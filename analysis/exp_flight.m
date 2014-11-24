

Trial = MTATrial('er06-20130614','all-cof');

display = true;
overwrite = true;
units = select_units(Trial,18,'pyr');
pfs = {};
states = {'u_wr','flight'};
nsts = size(states,2);

wrper = Trial.stc{'r'}.copy;
wrper.data = [Trial.stc{'r'}.data(1:2:end,:);Trial.stc{'w'}.data(1:2:end,:)];
wrper.data = sort(wrper.data);
wrper.key = 'x';
wrper.label = 'u_wr';
wrper.updateFilename([Trial.filebase '.sst.' wrper.label '.' wrper.key '.mat']);
Trial.stc.states{end+1} = wrper;

pfs{1} = MTAApfs(Trial,units,'u_wr',overwrite,'numIter',1,'binDims',[20,20,20],'type','xyz','SmoothingWeights',[1.8,1.8,1.8]);



Trial = MTATrial('er06-20130614','fly');
if isempty(Trial.stc{'t'}),Trial.stc.states{end+1} = theta(Trial);end


pfs{2} = MTAApfs(Trial,units,'flight&theta',overwrite,'numIter',1,'binDims',[20,20,20],'type','xyz','SmoothingWeights',[1.8,1.8,1.8]);







if display,

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        sp = [];
for i=1:nsts,

 
subplot2(nsts,2,i,1), pfs{i}.plot(unit,'xz',1,'isCircular',false); 
title([pfs{i}.session.trialName,':',pfs{i}.parameters.states,': ',num2str(unit)]);
subplot2(nsts,2,i,2), pfs{i}.plot(unit,'xy',1,'isCircular',false); 
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
