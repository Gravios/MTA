function generate_placefields(Trial,varargin)

% TESTARGS ------------------------------------------------------------------------------------------
% $$$ Trial = MTATrial.validate('er01-20110719.cof.all');
% $$$ Trial.load('stc','msnn_ppsvd_raux');
% $$$ units = [];
% $$$ states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',          ...
% $$$           'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};
% $$$ report = true;
% $$$ overwrite = false;
%---------------------------------------------------------------------------------------------------

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',               [],                                                      ...
                 'states',              {{'loc&theta','lloc&theta','hloc&theta','rear&theta',    ...
                                         'pause&theta','lpause&theta','hpause&theta',            ...
                                         'theta-groom-sit'}},                                    ...
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

pfk = {};
for sts = 1:nsts,
    da = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
    da.units = units;
    da.states = states{sts};
    da.overwrite = overwrite;
    da = struct2varargin(da);        
    pfk{sts} = MTAAknnpfs_bs(Trial,da{:});      
end


if report,
    OwnDir = '/storage/gravio/figures/parts/placefields';
    FigDir = ['placefields-',Trial.filebase];
    create_directory(fullfile(OwnDir,FigDir));
    
    hfig = figure;
    hfig.Units = 'centimeters';
    hfig.Position([1,2]) = [1,4];
    hfig.Position([3,4]) = [nsts*4+4,12];
    hfig.PaperPositionMode = 'auto';
    unit = units(1);
    autoIncr =true;
    while unit~=-1,
        clf();
        pmrs = cell2mat(cf(@(p) p.maxRate(unit), pfs));
        pmr = max(pmrs);
        for i = 1:nsts,
            hax = axes();
            hax.Units = 'centimeters';
            hax.Position = [(i-1)*4+2,2,1.5,1.5];
            pfs{i}.plot(unit,'mean',false,pmr,true);
            title([pfs{i}.session.trialName ' ' pfs{i}.parameters.states,': ',num2str(unit)]);
            hax.XTickLabel = {};
            hax.YTickLabel = {};            
            hax.XTick = [];            
            hax.YTick = [];
            xlabel(['max rate: ',num2str(pmrs(i))]);
        end

        pmrs = cell2mat(cf(@(p) p.maxRate(unit), pfk));
        pmr = max(pmrs);
        for i = 1:nsts,
            hax = axes();
            hax.Units = 'centimeters';
            hax.Position = [(i-1)*4+2,8,1.5,1.5];
            pfk{i}.plot(unit,'mean',false,pmr,true);
            title([pfs{i}.session.trialName ' ' pfs{i}.parameters.states,': ',num2str(unit)]);
            hax.XTickLabel = {};
            hax.YTickLabel = {};            
            hax.XTick = [];            
            hax.YTick = [];
            xlabel(['max rate: ',num2str(pmrs(i))]);
        end
        
        FigName = [FigDir,'_unit',num2str(unit)];
        print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
        print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));

        unit = figure_controls(hfig,unit,units,autoIncr);
    end

end

% END MAIN -----------------------------------------------------------------------------------------