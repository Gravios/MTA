function ds = batch_compute_pfstats_bs(varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',         'MjgER2016',                                             ...
                 'stcMode',             'msnn_ppsvd_raux',                                       ...
                 'states',              {{'loc','lloc','hloc','rear','pause','lpause','hpause'}},...
                 'tag',                 '',                                                      ...
                 'overwrite',           false,                                                   ...
                 'verbose',             true,                                                    ...
                 'Trial',               []                                                       ...
);%-------------------------------------------------------------------------------------------------

[sessionList,stcMode,states,tag,overwrite,verbose,Trial] = DefaultArgs(varargin,defargs,'--struct');

% modify states
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states{end+1} = 'theta';


% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash(struct('sessionList',  sessionList,...
                          'stcMode',      stcMode,...
                          'states',       {states}));
end
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

slist = get_session_list(sessionList);

% NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO !!!!!!!!!!!!!!!
if ~isempty(Trial),
    try,
        ds = load(fullfile(Trial.spath,[Trial.filebase,'_pfstatsBS_',tag,'.mat']));
        return;
    catch err,
        disp(err);
        slistInds = find(cell2mat(af(@(l,t) all([strcmp(l.sessionName,t{1}.name),...
                                                 strcmp(l.mazeName,t{1}.maze.name),...
                                                 strcmp(l.trialName,t{1}.trialName)]),...
                            slist,repmat({Trial},[1,numel(slist)]))));
    end
else
    slistInds = 1:numel(slist);
end


parp = [];
ds = {};
for s = slistInds

% LOAD trial meta data
    if isa(slist,'MTASession'), Trial = slist; else Trial = MTATrial.validate(slist(s)); end
    
% BUILD output file name 
    analysisFileName = fullfile(Trial.spath,[Trial.filebase,'_pfstatsBS_',tag,'.mat']);            
    
    if ~exist(analysisFileName,'file') || overwrite,    

% DISPLAY processing status
        if verbose, fprintf(['\nProcessing trial: %s\n','Output to: %s\n'],Trial.filebase,analysisFileName);end
        
% LOAD labeled behavior
        try,  Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
        catch err, disp(err),
            Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
        end

% REDUCE clu list based on theta pfs max rate        
        pft = pfs_2d_theta(Trial,[],[],overwrite);
        mrt = pft.maxRate;
        units = select_units(Trial,18);
        units = units(mrt(pft.data.clu(units))>1);

% COMPUTE shuffeled place fields in theta state
        defargs = get_default_args('MjgER2016','MTAAknnpfs','struct');
        defargs.units = units;
        defargs.states = 'theta-groom-sit';
        defargs.overwrite = false;
        defargs = struct2varargin(defargs);        
        pf = MTAAknnpfs(Trial,defargs{:});      

        
% COMPUTE place fields and subsampled estimate
        if verbose, fprintf('\nProcessing placefields...\n'); end
        for sts = 1:numel(states),
            if verbose, fprintf('process state: %s...\n',states{sts}); end            
            defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
            defargs.units = units;
            defargs.states = states{sts};
            defargs = struct2varargin(defargs);        
            pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
        end
        
% DISPLAY processing status        
        if verbose, fprintf('\nParse place field features...\n'); end

% PARSE place field features
% COMPILE placefield statistics
        try delete(gcp('nocreate')); end
        parp = parpool(7);
        cluMap = pfkbs{1}.data.clu;
        pfkstats = {};        pfkboots = {};        pfmstats = {};
        tic
        parfor t = 1:numel(states),
        %for t = 1:numel(states),            
            for u = 1:numel(units),
                [pfkstats{t}{u},pfkboots{t}{u},pfmstats{t}{u}] = PlaceFieldStats(Trial,pfkbs{t},units(u),false);
            end
        end
        toc
        try,delete(gcp('nocreate'));end            
        
% FUCK matlab
        % Fuck matlab ... seriously 
        pfkstats = [pfkstats{:}];
        pfkstats = [pfkstats{:}];
        pfmstats = [pfmstats{:}];
        pfmstats = [pfmstats{:}];
        pfkboots = [pfkboots{:}];
        pfkboots = [pfkboots{:}];
        pfkstats = reshape(pfkstats,numel(units),numel(states))';
        pfmstats = reshape(pfmstats,numel(units),numel(states))';
        pfkboots = reshape(pfkboots,numel(units),numel(states))';

% DISPLAY processing status                
        if verbose, fprintf('Extract greatest place feature...\n'); end
% SELECT the biggest baddest place field patch for all units
        peakPatchArea = [];
        peakPatchCOM  = [];
        peakPatchRate = [];
        for t = 1:numel(states),
            for k = 1:pfkbs{1}.parameters.numIter,
% RETRIEVE the patch center of mass from patch with the highest firing rate
                pcom = ...
                    ... %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                arrayfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                         pfkboots(t,:),...
                         repmat(k,[1,numel(units)]),...
                         'UniformOutput',false);

                pind = ~cellfun(@isempty,pcom);
                peakPatchCOM(t,k,pind,:) = ...
                    cell2mat(cellfun(@(x) x(:,1),...
                                     pcom(pind),...
                                     'uniformoutput',false))';

% RETRIEVE the max patch Firing rate
                peakPatchRate(t,k,:) = ...
                    sq(arrayfun(@(x,y) max(x.patchPFR(1,y,:)),...
                                pfkboots(t,:),...
                                repmat(k,[1,numel(units)])));

% RETRIEVE the patch area from patch with highest firing rate
                parea = ...
                    arrayfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                             pfkboots(t,:),...
                             repmat(k,[1,numel(units)]),...
                             'UniformOutput',false);
                pind = ~cellfun(@isempty,parea);
                peakPatchArea(t,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                              parea(pind),...
                                                              'uniformoutput',false))');
            end
        end


% SAVE pfstats
        if verbose, fprintf('\nSaving pfstats.\n'); end
        session = Trial.filebase;
        save(analysisFileName,'-v7.3',...
             'session','stcMode','states','cluMap',...
             'pfkstats','pfkboots','pfmstats',...
             'peakPatchArea','peakPatchCOM','peakPatchRate');
    end % IF ~exist(analysisFileName,'file') || overwrite,    

% DISPLAY processing status                    
    if verbose, fprintf('\nLoading %s...\n\n',analysisFileName); end
% LOAD pfstats
    ds(s) = {load(analysisFileName)};
end    

% END MAIN -------------------------------------------------------------------------------------------

