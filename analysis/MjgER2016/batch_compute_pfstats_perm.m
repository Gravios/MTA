
%subjects = {'jg04','jg05','ER06'};
%subjects = {'ER06_BHV','Ed10_BHV'};
%subjects = {'Ed10};
%subjects = {'jg05'};
slist = get_session_list('MjgEdER2016_bhv');
tag = '';
tag = 'bhvtheta';
if ~isempty(tag),tag = ['_',tag];end

stcMode = 'msnn_ppsvd';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
overwrite = false;
parp = [];

ds = {};


for s = 1:numel(slist)

    % label behavior
    Trial = MTATrial.validate(slist(s));

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
    units = select_units(Trial,18);
    units = units(mrt(units)>1);

    % Compute place fields and subsampled estimate
    pfkbs = {};
    for sts = 1:numel(states),
        for sti = sts+1:numel(states),
            defargs = get_default_args('MjgER2016','MTAAknnpfs_perm','struct');
            defargs.units = units;
            defargs.states = [states(sts),states(sti)];
            defargs = struct2varargin(defargs);        
            pfkbs{end+1} = MTAAknnpfs_perm(Trial,defargs{:});      
        end
    end

% PARSE place field features

    analysisFileName = fullfile(Trial.spath,[Trial.filebase,tag,'_pfStats.mat']);        
    
    if 1,%~exist(analysisFileName,'file') || overwrite,
        
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

        
        
% FUCK matlab ... seriously 
        pfkstats = [pfkstats{:}];
        pfkstats = [pfkstats{:}];

        pfmstats = [pfmstats{:}];
        pfmstats = [pfmstats{:}];

        pfkboots = [pfkboots{:}];
        pfkboots = [pfkboots{:}];

        pfkstats = reshape(pfkstats,numel(units),numel(states))';
        pfmstats = reshape(pfmstats,numel(units),numel(states))';
        pfkboots = reshape(pfkboots,numel(units),numel(states))';
        
        
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

        


        session = Trial.filebase;
        save(analysisFileName,'-v7.3',...
             'session','stcMode','states','cluMap',...
             'pfkstats','pfkboots','pfmstats',...
             'peakPatchArea','peakPatchCOM','peakPatchRate');
        
    end

    ds(end+1) = {load(analysisFileName)};

end    


