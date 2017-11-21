function ds = compute_pfstats_bs(Trial,varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',         'MjgER2016',                                             ...
                 'stcMode',             'msnn_ppsvd_raux',                                       ...
                 'states',              {{'loc&theta','lloc&theta','hloc&theta','rear&theta',    ...
                                          'pause&theta','lpause&theta','hpause&theta',           ...
                                          'theta-groom-sit'}},                                   ...
                 'tag',                 '',                                                      ...
                 'overwrite',           false,                                                   ...
                 'verbose',             true                                                     ...
);
[sessionList,stcMode,states,tag,overwrite,verbose] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash(struct('sessionList',  sessionList,...
                          'stcMode',      stcMode,...
                          'states',       {states}));
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------


% BUILD output file name 
analysisFileName = fullfile(Trial.spath,[Trial.filebase,'_pfstatsBS_',tag,'.mat']);            


if ~exist(analysisFileName,'file') || overwrite,    

% DISPLAY processing status
    if verbose, fprintf(['\nProcessing trial: %s\n','Output to: %s\n'],Trial.filebase,analysisFileName);end
    
% LOAD labeled behavior
    Trial.load('stc',stcMode);

% REDUCE clu list based on theta pfs max rate        
    pft = pfs_2d_theta(Trial,[],[],true); % overwrite
    mrt = pft.maxRate;
    units = select_units(Trial,18);
    units = units(mrt(pft.data.clu(units))>1);

% COMPUTE shuffeled place fields in theta state
    defargs = get_default_args('MjgER2016','MTAAknnpfs','struct');
    defargs.units = units;
    defargs.states = 'theta-groom-sit';
    defargs = struct2varargin(defargs);        
    pf = MTAAknnpfs(Trial,defargs{:});      
    
% COMPUTE place fields and subsampled estimate
    if verbose, fprintf('\nProcessing placefields...\n'); end
    for sts = 1:numel(states),
        if verbose, fprintf('process state: %s...\n',states{sts}); end            
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = units;
        %defargs.overwrite = true;
        defargs.states = states{sts};
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end
    
% DISPLAY processing status        
    if verbose, fprintf('\nParse place field features...\n'); end

% PARSE place field features
% COMPILE placefield statistics
    cluMap = units;
    pfkstats = {};        pfkboots = {};        pfmstats = {};
    if 0,%isempty(gcp('nocreate')),
        parp = parpool(8);
        parfor t = 1:numel(states),
            for u = 1:numel(units),
                [pfkstats{t}{u},pfkboots{t}{u},pfmstats{t}{u}] = ...
                    PlaceFieldStats(Trial,pfkbs{t},units(u),false);
            end % units loop
        end % states loop
        delete(gcp('nocreate'));
    else
        for t = 1:numel(states),
            for u = 1:numel(units),
                [pfkstats{t}{u},pfkboots{t}{u},pfmstats{t}{u}] = ...
                    PlaceFieldStats(Trial,pfkbs{t},units(u),false);
            end % units loop
        end % states loop
    end
        
    
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
        end % iteration loop
    end % states loop


% SAVE pfstats
    if verbose, fprintf('\nSaving pfstats.\n'); end
    session = Trial.filebase;
    save(analysisFileName,'-v7.3',...
         'session','stcMode','states','cluMap',...
         'pfkstats','pfkboots','pfmstats',...
         'peakPatchArea','peakPatchCOM','peakPatchRate');
end % if savefile exists

% DISPLAY processing status                    
if verbose, fprintf('\nLoading %s...\n\n',analysisFileName); end
% LOAD pfstats
ds = load(analysisFileName);


% END MAIN -------------------------------------------------------------------------------------------

