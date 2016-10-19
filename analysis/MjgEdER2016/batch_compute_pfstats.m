function ds = batch_compute_pfstats(varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',         'MjgEdER2016_bhv',                                       ...
                 'stcMode',             'NN0317R',                                               ...
                 'states',              {{'loc','lloc','hloc','rear','pause','lpause','hpause'}},...
                 'tag',                 '',                                                      ...
                 'overwrite',           false                                                    ...
);%-------------------------------------------------------------------------------------------------

[sessionList,stcMode,states,tag,overwrite] = DefaultArgs(varargin,defargs,'--struct');

% modify states
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);




% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag
%
%    sessionList
%    stcMode 
%    states 
if isempty(tag),
    tag = DataHash(struct('sessionList',sessionList,'stcMode',stcMode,'states',{states}));
end
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

slist = get_session_list(sessionList);
parp = [];
ds = {};
for s = 1:numel(slist)

    Trial = MTATrial.validate(slist(s));

    analysisFileName = fullfile(Trial.spath,[Trial.filebase,'_pfstats_',tag,'.mat']);            
    
    if ~exist(analysisFileName,'file') || overwrite,    
        
        % load labeled behavior
        try,
            Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
        catch err
            disp(err)
            Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
        end

        pft = pfs_2d_theta(Trial,[],[],overwrite);
        mrt = pft.maxRate;

        % Reduce clu list based on theta pfs max rate
        units = select_units(Trial,18);
        units = units(mrt(pft.data.clu(units))>1);

        % Compute place fields and subsampled estimate
        for sts = 1:numel(states),
            defargs = get_default_args_MjgEdER2016('MTAApfs','struct');
            defargs.units = units;
            defargs.states = states{sts};
            defargs = struct2varargin(defargs);        
            pfkbs{sts} = MTAApfs(Trial,defargs{:});      
        end

        % Parse place field features
        cluMap = pfkbs{1}.data.clu;

        pfkstats = {};

        tic
        for t = 1:numel(states),
            for u = 1:numel(units),
                [pfkstats{t}{u},pfkboots{t}{u},pfmstats{t}{u}] = PlaceFieldStats(Trial,pfkbs{t},units(u),false);
            end
        end
        toc
        
        
        % Fuck matlab ... seriously 
        pfkstats = [pfkstats{:}];
        pfkstats = [pfkstats{:}];

        pfkstats = reshape(pfkstats,numel(units),numel(states))';
        
        
        % Select the biggest baddest place field patch for all units
        peakPatchArea = [];
        peakPatchCOM  = [];
        peakPatchRate = [];

% $$$         for t = 1:numel(states),
% $$$             for k = 1:pfkbs{1}.parameters.numIter,
% $$$                 % Retrieve the patch center of mass from patch with the highest firing rate
% $$$                 pcom = ...
% $$$                     ... %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
% $$$                 arrayfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
% $$$                          pfkboots(t,:),...
% $$$                          repmat(k,[1,numel(units)]),...
% $$$                          'UniformOutput',false);
% $$$ 
% $$$                 pind = ~cellfun(@isempty,pcom);
% $$$                 peakPatchCOM(t,k,pind,:) = ...
% $$$                     cell2mat(cellfun(@(x) x(:,1),...
% $$$                                      pcom(pind),...
% $$$                                      'uniformoutput',false))';
% $$$ 
% $$$                 % Retrieve the max patch Firing rate
% $$$                 peakPatchRate(t,k,:) = ...
% $$$                     sq(arrayfun(@(x,y) max(x.patchPFR(1,y,:)),...
% $$$                                 pfkboots(t,:),...
% $$$                                 repmat(k,[1,numel(units)])));
% $$$ 
% $$$                 % Retrieve the patch area from patch with highest firing rate
% $$$                 parea = ...
% $$$                     arrayfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
% $$$                              pfkboots(t,:),...
% $$$                              repmat(k,[1,numel(units)]),...
% $$$                              'UniformOutput',false);
% $$$                 pind = ~cellfun(@isempty,parea);
% $$$                 peakPatchArea(t,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
% $$$                                                               parea(pind),...
% $$$                                                               'uniformoutput',false))');
% $$$             end
% $$$         end
% $$$ 

        session = Trial.filebase;
        save(analysisFileName,'-v7.3',...
             'session','stcMode','states','cluMap',...
             'pfkstats','pfkboots','pfmstats',...
             'peakPatchArea','peakPatchCOM','peakPatchRate');
    end

    ds(end+1) = {load(analysisFileName)};

end    




% $$$         warning('MTAStateCollection not found, running label_bhv_nn.m')
% $$$         labelBhv_NN(Trial,...
% $$$                     stcMode,...
% $$$                     'jg05-20120317.cof.all',...
% $$$                     'hl_3_jg_r',...
% $$$                     [],[],[],[],[],[],[],{'loc','rear','pause','groom','sit'});
