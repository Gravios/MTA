function gen_pfs(Trial,varargin)
[units,display,overwrite,tnames,maze,states,mode,niter] = DefaultArgs(varargin,{Trial.spk.map(:,1),false,false,{'all'},'cof',{'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'},'std',10000});


%stc = Trial.stc.copy; 
 
Session = Trial.name;
tnames = {'lof','lon'};
maze = 'cof';
nt = numel(tnames);
stc = Trial.stc.copy;

%states = {'flight'};
%states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta','hswalk&theta','lswalk&theta'};
states = {'theta','rear','walk','hang','lang','hswalk','lswalk'};
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

