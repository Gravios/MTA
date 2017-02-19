function ds = batch_compute_pfstats_bs(varargin)


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',         'MjgEdER2016_bhv',                                       ...
                 'stcMode',             'NN0317R',                                               ...
                 'states',              {{'loc','rear','pause','lloc','hloc','lpause','hpause'}},...
                 'tag',                 '',                                                      ...
                 'overwrite',           false,                                                   ...
                 'verbose',             true                                                     ...
);%-------------------------------------------------------------------------------------------------

[sessionList,stcMode,states,tag,overwrite,verbose] = DefaultArgs(varargin,defargs,'--struct');

% modify states
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);




% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create 
%sessionList
%stcMode 
%states 
if isempty(tag),
    tag = DataHash(struct('sessionList','MjgEdER2016_bhv','stcMode',stcMode,'states',{states}));
end
%---------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

slist = get_session_list(sessionList);
parp = [];
ds = {};
for s = 1:numel(slist)

    % Load trial meta data
    if isa(slist,'MTASession'),
        Trial = slist;
    else
        Trial = MTATrial.validate(slist(s));
    end
    
    % Build output file name 
    analysisFileName = fullfile(Trial.spath,[Trial.filebase,'_pfstatsBS_',tag,'.mat']);            
    
    if ~exist(analysisFileName,'file') || overwrite,    

        % Display processing status
        if verbose,
            fprintf(['\nProcessing trial: %s\n',...
                     'Output to: %s\n'],Trial.filebase,analysisFileName);
        end
        
        
        
        % load labeled behavior
        try,
            Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
        catch err
            disp(err)
            Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
        end



        % Reduce clu list based on theta pfs max rate        
        pft = pfs_2d_theta(Trial,[],[],overwrite);
        mrt = pft.maxRate;
        units = select_units(Trial,18);
        units = units(mrt(pft.data.clu(units))>1);


        if verbose, fprintf('\nProcessing placefields...\n'); end
        % Compute place fields and subsampled estimate
        for sts = 1:numel(states),
            if verbose, fprintf('process state: %s...\n',states{sts}); end            
            defargs = get_default_args('MjgEdER2016','MTAAknnpfs_bs','struct');
            defargs.units = units;
            defargs.states = states{sts};
            defargs = struct2varargin(defargs);        
            pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
        end

        
        
        if verbose, fprintf('\nParse place field features...\n'); end
        % Parse place field features
        try delete(gcp('nocreate')); end
        parp = parpool(7);

        cluMap = pfkbs{1}.data.clu;

        pfkstats = {};
        pfkboots = {};
        pfmstats = {};

        tic
        parfor t = 1:numel(states),
            for u = 1:numel(units),
                [pfkstats{t}{u},pfkboots{t}{u},pfmstats{t}{u}] = PlaceFieldStats(Trial,pfkbs{t},units(u),false);
            end
        end
        toc
        try,delete(gcp('nocreate'));end            

        
              
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

        
        if verbose, fprintf('Extract greatest place feature...\n'); end
        % Select the biggest baddest place field patch for all units
        peakPatchArea = [];
        peakPatchCOM  = [];
        peakPatchRate = [];
        for t = 1:numel(states),
            for k = 1:pfkbs{1}.parameters.numIter,
                % Retrieve the patch center of mass from patch with the highest firing rate
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

                % Retrieve the max patch Firing rate
                peakPatchRate(t,k,:) = ...
                    sq(arrayfun(@(x,y) max(x.patchPFR(1,y,:)),...
                                pfkboots(t,:),...
                                repmat(k,[1,numel(units)])));

                % Retrieve the patch area from patch with highest firing rate
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


        % Save pfstats
        if verbose, fprintf('\nSaving pfstats.\n'); end
        session = Trial.filebase;
        save(analysisFileName,'-v7.3',...
             'session','stcMode','states','cluMap',...
             'pfkstats','pfkboots','pfmstats',...
             'peakPatchArea','peakPatchCOM','peakPatchRate');
    end

    if verbose, fprintf('\nLoading %s...\n\n',analysisFileName); end
    ds(s) = {load(analysisFileName)};

end    

% END MAIN -------------------------------------------------------------------------------------------


% $$$         warning('MTAStateCollection not found, running label_bhv_nn.m')
% $$$ labelBhv_NN(Trial,...
% $$$             stcMode,...
% $$$             'jg05-20120317.cof.all',...
% $$$             'hl_3_jg_r',...
% $$$             [],[],[],[],[],[],[],{'loc','rear','pause','groom','sit'});
% $$$ 
% $$$ label_aux_bhv_reduced(Trial,stcMode,'overwrite',overwrite);        
