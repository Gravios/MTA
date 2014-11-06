Session = 'Ed10-20140815';
tnames = {'lof','lon'};
maze = 'cof';
nt = numel(tnames);

Trial = MTATrial('Ed10-20140815');
Trial = QuickTrialSetup(s,'lof',[],4:7,false);
Trial = QuickTrialSetup(s,'lon',[],1:3,false);


%states = {'flight'};
states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
nsts = size(states,2);


display = true;
mode = 'std';
niter = 10000;
overwrite = true;
%units = select_units(Trial,18,'pyr');
units = 1:100;

Trial = MTATrial(Session,'all',maze);
stc = Trial.stc.copy;

pfs = {};
for t = 1:nt,
    Trial = MTATrial(Session,tnames{t},maze);
    Trial.stc = stc.copy;
    Trial.stc.load(Trial);
    switch mode
      case 'std'
        for i = 1:nsts,
            pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
        end

% $$$       case 'shuff'
% $$$         for i = 1:nsts,
% $$$             MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',niter,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',125,'nNearestNeighbors',150);
% $$$         end
    end
end



if display,

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for t = 1:nt,
        for i = 1:nsts,
            subplot2(nt,nsts,t,i);cla
            %subplotfit(i,nsts+1);cla
            pfs{t,i}.plot(unit,[],true);
            title([pfs{t,i}.session.trialName ' ' pfs{t,i}.parameters.states,': ',num2str(unit)]);
        end
        end
        unit = figure_controls(hfig,unit,units);
    end

end

