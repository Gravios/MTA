
Session = Trial.name;
tnames = {'low','high'};
maze = 'sof';
nt = numel(tnames);


%states = {'flight'};
%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'theta'};%,'rear','walk','hang','lang','hswalk','lswalk'};
nsts = size(states,2);


% $$$ display = true;
% $$$ mode = 'std';
% $$$ niter = 10000;
% $$$ overwrite = true;
% $$$ %units = select_units(Trial,18,'pyr')v
;
units = 1:100;


pfs = {};
for t = 1:nt,
    Trial = MTATrial(Session,tnames{t},maze);
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
            %subplotfit(i,nsts+1);cla
            try,pfs{t,i}.plot(unit,[],true);
            title([pfs{t,i}.session.trialName ' ' pfs{t,i}.parameters.states,': ',num2str(unit)]);
        end
        end
        end
        reportfig('/gpfs01/sirota/home/gravio/figures/',hfig,...
                  ['exp_zshift-' Trial.name] ,false,Trial.name);

        unit = figure_controls(hfig,unit,units);
    end

end
