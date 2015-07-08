
Trial = MTATrial('jg05-20120310')

%states = {'flight'};
% $$$ states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'theta','theta-ripples'};

nsts = size(states,2);


display = true;
mode = 'std';
niter = 10000;
overwrite = true;
units = select_units(Trial,18,'pyr');
%units = 1:21;


pfs = {};
switch mode
  case 'std'
    for i = 1:nsts,
        pfs{i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
    end

% $$$       case 'shuff'
% $$$         for i = 1:nsts,
% $$$             MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
% $$$         end
end




if display,

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for i = 1:nsts,
            subplot2(1,nsts,1,i);cla
            %subplotfit(i,nsts+1);cla
            pfs{i}.plot(unit,[],true);
            title([pfs{i}.session.trialName ' ' pfs{i}.parameters.states,': ',num2str(unit)]);
        end
        unit = figure_controls(hfig,unit,units);
    end

end

