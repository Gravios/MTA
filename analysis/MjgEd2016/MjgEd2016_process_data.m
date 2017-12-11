
sessionList = get_session_list('MjgER2016');
stcMode     = 'msnn_ppsvd_raux';
states      ={'loc&theta','lloc&theta','hloc&theta','rear&theta',...
              'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};

for t = 1:numel(sessionList),
for t = 5:10
    
    Trial = MTATrial.validate(sessionList(t));
    disp(['processing trial: ' Trial.filebase]);
    
    try,  Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
    catch err, disp(err),
        Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
    end

    pfstats = compute_pfstats_bs(Trial);
    
    for sts = 1:numel(states)-1,
        fprintf('processing state: %s...\n',states{sts});
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = pfstats.cluMap;
        defargs.overwrite = true;
        defargs.states = states{sts};
        defargs.nNodes = 12;
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end
end
    
    
