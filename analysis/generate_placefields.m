function generate_placefields(Trial,varargin)

% TESTARGS ------------------------------------------------------------------------------------------
Trial = MTATrial.validate('er01-20110719.cof.all');
Trial.load('stc','msnn_ppsvd_raux');
units = [];
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',          ...
          'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};
report = true;
overwrite = false;
%---------------------------------------------------------------------------------------------------

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',               [],                                                      ...
                 'states',              {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
                                         'pause&theta','lpause&theta','hpause&theta',            ...
                                         'theta-groom-sit'},                                     ...
                 'report',              true,                                                    ...
                 'overwrite',           false                                                    ...                 
);
[units,states,report,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

if isempty(units),
    units = select_placefields(Trial);
end

nsts = numel(states);

pfs = {};
for sts = 1:nsts,
    da = get_default_args('MjgER2016','MTAApfs','struct');
    da.units = units;
    da.states = states{sts};
    da.overwrite = overwrite;
    da = struct2varargin(da);        
    pfs{sts} = MTAApfs(Trial,da{:});      
end


if report,
    OwnDir = '/storage/gravio/figures/parts/placefields';
    FigDir = ['placefields',Trial.filebase];
    create_directory(fullfile(OwnDir,FigDir));
    
    hfig = figure;
    hfig.Units = 'centimeters';
    hfig.Position([1,2]) = [1,20];
    hfig.Position([3,4]) = [nsts*4,6];
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

% END MAIN -----------------------------------------------------------------------------------------