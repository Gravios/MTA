


% label behavior
Trial = MTATrial.validate(sessionList(s));

try,  Trial.load('stc',stcMode); 
catch err
    warning('MTAStateCollection not found, running label_bhv_nn.m')
    labelBhv_NN(Trial,...
                stcMode,...
                'jg05-20120317.cof.all',...
                'hl_3_jg_r',...
                [],[],[],[],[],[],[],{'loc','rear','pause','groom','sit'});
    Trial.load('stc',stcMode);
end

label_aux_bhv_reduced(Trial,stcMode,'overwrite',overwrite);        

pft = pfs_2d_theta(Trial,[],[],overwrite);
mrt = pft.maxRate;

% Reduce clu list based on theta pfs max rate
units = pft.data.clu(mrt>1);

% Compute place fields and subsampled estimate
for sts = 1:numel(states),
    defargs = get_default_args_MjgEdER2016('MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end

% Parse place field features

analysisFileName = fullfile(Trial.spath,[Trial.filebase,tag,'_pfStats.mat']);        

if 1%~exist(analysisFileName,'file') || overwrite,
    
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
    
    
    pfkstats = reshape([pfkstats{:}],numel(units),numel(states))';
    pfmstats = reshape([pfkstats{:}],numel(units),numel(states))';
    pfkboots = reshape([pfkboots{:}],numel(units),numel(states))';


% $$$             
% $$$             pfkstats = reshape(pfkstats,numel(units),numel(states))';
% $$$             pfmstats = reshape(pfkstats,numel(units),numel(states))';
% $$$             pfkboots = reshape(pfkboots,numel(units),numel(states))';
    
    peakPatchArea = [];
    peakPatchCOM  = [];
    peakPatchRate = [];
    % Select the biggest baddest place field patch for all units
    for t = 1:numel(states),
        for k = 1:pfkbs{1}.parameters.numIter,
            % Retrieve the patch center of mass from patch with the highest firing rate
            pcom = ...
                ... %cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
            cellfun(@(x,y) sq(x.patchCOM(1,y,find(max(x.patchArea(1,y,:).*x.patchPFR(1,y,:))==x.patchArea(1,y,:).*x.patchPFR(1,y,:)),:)),...
                    pfkboots(t,:),...
                    repmat({k},[1,numel(units)]),...
                    'UniformOutput',false);

            pind = ~cellfun(@isempty,pcom);
            peakPatchCOM(t,k,pind,:) = ...
                cell2mat(cellfun(@(x) x(:,1),...
                                 pcom(pind),...
                                 'uniformoutput',false))';

            
            % Retrieve the max patch Firing rate
            peakPatchRate(t,k,:) = ...
                sq(cellfun(@(x,y) max(x.patchPFR(1,y,:)),...
                           pfkboots(t,:),...
                           repmat({k},[1,numel(units)])));

            
            % Retrieve the patch area from patch with highest firing rate
            parea = ...
                cellfun(@(x,y) sq(x.patchArea(1,y,find(max(x.patchPFR(1,y,:))==x.patchPFR(1,y,:)),:)),...
                        pfkboots(t,:),...
                        repmat({k},[1,numel(units)]),...
                        'UniformOutput',false);
            pind = ~cellfun(@isempty,parea);
            peakPatchArea(t,k,pind) = sq(cell2mat(cellfun(@(x) x(1),...
                                                          parea(pind),...
                                                          'uniformoutput',false))');
        end
    end



    session = Trial.filebase;
    save(analysisFileName,'-v7.3',...
         'session','stcMode','states','cluMap',...
         'pfkstats','pfkboots','pfmstats',...
         'peakPatchArea','peakPatchCOM','peakPatchRate');
else
    ds(end+1) = {load(analysisFileName)};
end
