

Trial = MTATrial('er06-20130614','all-cof');

display = true;
%overwrite = true;
units = select_units(Trial,18,'pyr');
pfs = {};
states = {'theta','flight'};
nsts = size(states,2);

% $$$ wrper = Trial.stc{'r'}.copy;
% $$$ wrper.data = [Trial.stc{'r'}.data(1:2:end,:);Trial.stc{'w'}.data(1:2:end,:)];
% $$$ wrper.data = sort(wrper.data);
% $$$ wrper.key = 'x';
% $$$ wrper.label = 'u_wr';
% $$$ wrper.updateFilename([Trial.filebase '.sst.' wrper.label '.' wrper.key '.mat']);
% $$$ Trial.stc.states{end+1} = wrper;

Trial.maze.boundaries(end) = 450;
pfs{1} = MTAApfs(Trial,units,'theta',overwrite,'numIter',1,'binDims',[40,40,30],'type','xyz','SmoothingWeights',[1.5,1.5,1.2]);



Trial = MTATrial('er06-20130614','fly');
if isempty(Trial.stc{'t'}),Trial.stc.states{end+1} = theta(Trial);end
Trial.maze.boundaries(end) = 450;

pfs{2} = MTAApfs(Trial,units,'theta-flight',overwrite,'numIter',1,'binDims',[40,40,30],'type','xyz','SmoothingWeights',[1.5,1.5,1.2]);


pfs{3} = MTAApfs(Trial,units,'flight&theta',overwrite,'numIter',1,'binDims',[40,40,30],'type','xyz','SmoothingWeights',[1.5,1.5,1.2]);







if display,

    
    hfig = figure;
    unit = units(1);
set(hfig,'paperposition',get(hfig,'paperposition').*[1,1,0,0]+[0,0,9,12])
    while unit~=-1,
        sp = [];

for i=1:numel(pfs),
mrxz(i) =  max(max(max(pfs{i}.plot(unit,'xz',1,'isCircular',false))));
end
Nmrxz = max(mrxz);
        for i=1:numel(pfs),
subplot2(numel(pfs),2,i,1), pfs{i}.plot(unit,'xz',1,[0,mrxz],'isCircular',false); 
title([pfs{i}.session.trialName,':',pfs{i}.parameters.states,': ',num2str(unit)]);
subplot2(numel(pfs),2,i,2), pfs{i}.plot(unit,'xy',1,[0,mrxz],'isCircular',false); 
end
% $$$         for i = 1:nsts,
% $$$             sp(i) = subplotfit(i,nsts);cla
% $$$             pfs{i}.plot(unit,'isosurface')%,true,[0,max(mrate)]);
% $$$             title([pfs{i}.parameters.states,': ',num2str(unit)]);
% $$$         end
% $$$         linkprop(sp, {'CameraPosition'}); 
% $$$ 
    saveas(gcf,fullfile('/gpfs01/sirota/home/gravio/',...
                        'figures','SFN2014',...
                        ['pfsFly2_' Trial.filebase '-' num2str(unit) '.png']),'png');

        unit = figure_controls(hfig,unit,units);
    end

end
