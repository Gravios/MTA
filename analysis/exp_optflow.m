

Trials = SessionList('Ed10VR_opticflow');
nt = numel(Trials);

states = {'theta'};
nsts = size(states,2);

units = 1:100;
overwrite = true;

pfs = {};
for t = 1:nt,
    Trial = MTATrial(Trials(t).name,...
                     Trials(t).trialName,...
                     Trials(t).mazeName);
    Trial.stc.states = {};
    Trial.stc.states{1} = theta(Trial);
    for i = 1:nsts,
        pfs{t,i} = MTAApfs(Trial,units,states{i},overwrite,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
    end
end



if display,

    hfig = figure;
    unit = units(1);
    while unit~=-1,
        for t = 1:nt,
        for i = 1:nsts,
            subplot2(nt,nsts,t,i);cla
            try,pfs{t,i}.plot(unit,[],true,[],false);
            title([pfs{t,i}.session.trialName ' ' pfs{t,i}.parameters.states,': ',num2str(unit)]);
        end
        end
        end
        %reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
        %          ['exp_optflow-' Trial.name] ,false,Trial.name);

        unit = figure_controls(hfig,unit,units);
    end

end
