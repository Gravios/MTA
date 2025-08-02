
sessionList = get_session_list('MjgER2016');
stcMode     = 'msnn_ppsvd_raux';
states      ={'loc&theta','lloc&theta','hloc&theta','rear&theta',...
              'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};
overwrite = true;

for t = 1:numel(sessionList),
    Trial = MTATrial.validate(sessionList(t));
    disp(['processing trial: ' Trial.filebase]);
    
    try,  Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
    catch err, disp(err),
        Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
    end

    pfstats = compute_pfstats_bs(Trial);

% COMPLETE 
    for sts = 1:numel(states)-1
        defargs = get_default_args('MjgER2016','MTAAknnpfs','struct');
        defargs.units = pfstats.cluMap;
        defargs.overwrite = overwrite;
        defargs.states = states{sts};
        defargs.nNodes = 12;
        defargs = struct2varargin(defargs);        
        pftbs = MTAAknnpfs(Trial,defargs{:});      
    end
end

    unit = 15;
    mmap = pftbs.plot(unit,'mean');
    smap = pftbs.plot(unit,'std');
    threshold = nanmean(mmap(:))+nanmean(smap(:))*3;

    
% COMPLETE 
    for sts = 1:numel(states),
        fprintf('processing state: %s...\n',states{sts});
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = pfstats.cluMap;
        defargs.overwrite = overwrite;
        defargs.states = states{sts};
        defargs.nNodes = 12;
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end


figure();
for unit = pfstats.cluMap,
clf();
for sts = 1:8,
    rateMap = 1./nansum((repmat(max(pftbs.data.rateMap(:,pftbs.data.clu==unit,2:end)),...
                                [size(pftbs.data.rateMap,1),1,1])...
                         -repmat(nanmean(pfkbs{sts}.data.rateMap(:,pfkbs{sts}.data.clu==unit,:),3),...
                                 [1,1,pfkbs{sts}.parameters.numIter-1]))<0,3)';
subplot2(2,8,1,sts);pfkbs{sts}.plot(unit); title(num2str(unit));
subplot2(2,8,2,sts);imagesc(reshape(rateMap,[50,50])),axis('xy');
end
waitforbuttonpress();
end


% RUNNING  
    statePairs = {{'loc&theta','rear&theta'}};%,{'hloc','lloc'},{'hloc','lloc'}};
    for sts = 1:numel(statePairs)
        %fprintf('processing state: %s...\n',states{sts});
        defargs = get_default_args('MjgER2016','MTAAknnpfs_perm','struct');
        defargs.units = pfstats.cluMap;
        defargs.overwrite = true;
        defargs.states = statePairs{sts}
        defargs.nNodes = 12;
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_perm(Trial,defargs{:});      
    end

    
end
    
    


