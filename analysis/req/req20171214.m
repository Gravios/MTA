% Place fields in non theta periods minus ripples 


Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial = MTATrial.validate('er01-20110719.cof.all');
Trial.load('stc','msnn_ppsvd_raux');

%label_ripples(Trial);
pft = pfs_2d_theta(Trial,'overwrite',true);

units = select_placefields(Trial);

% LOAD placefields MTAApfs LIA
defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'gper-theta-spw';
defargs = struct2varargin(defargs);        
pfn = MTAApfs(Trial,defargs{:});      

% LOAD placefields MTAApfs THETA
defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'theta-groom-sit';
defargs.overwrite = false;
defargs = struct2varargin(defargs);        
pft = MTAApfs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = Trial.spk.map(:,1);
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs = struct2varargin(defargs);        
pfn = MTAApfs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.SmoothingWeights = [2.5,2.5];
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs = struct2varargin(defargs);        
pfm = MTAApfs(Trial,defargs{:});      



figure();
for unit = units,%pft.data.clu,
    clf();
    subplot(131);
    plot(pft,unit);
    subplot(132);
    %plot(pfk,unit,'mean');
% $$$     subplot(133);
% $$$     plot(pfm,unit);
    title(['unit: ',num2str(unit),' si:',num2str(pft.data.si(pft.data.clu==unit)),' spar:',num2str(pft.data.spar(pft.data.clu==unit))]);    
    waitforbuttonpress();
end




  


figure();
for unit = gunits,%pft.data.clu,
    clf();
    plot(pft,unit);
        title(['unit: ',num2str(unit),' si:',num2str(pft.data.si(pft.data.clu==unit)),' spar:',num2str(pft.data.spar(pft.data.clu==unit))]);
    waitforbuttonpress();
end




% LOAD placefields MTAAknnpfs_bs 
  
defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
defargs.units = units;
defargs.states = 'theta-groom-sit';
defargs = struct2varargin(defargs);        
pfkt = MTAAknnpfs_bs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
defargs.units = units;
defargs.states = 'loc&theta';
defargs = struct2varargin(defargs);        
pfkl = MTAAknnpfs_bs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
defargs.units = units;
defargs.states = 'pause&theta';
defargs = struct2varargin(defargs);        
pfkp = MTAAknnpfs_bs(Trial,defargs{:});      



% LOAD placefields MTAApfs 

overwrite = true;
defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'theta-groom-sit';
defargs.numIter = 1000;
defargs.halfsample = true;
defargs.overwrite = overwrite;
defargs = struct2varargin(defargs);        
pft = MTAApfs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'loc&theta';
defargs.numIter = 1000;
defargs.halfsample = true;
defargs.overwrite = overwrite;
defargs = struct2varargin(defargs);        
pfh = MTAApfs(Trial,defargs{:});      

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units = units;
defargs.states = 'pause&theta';
defargs.numIter = 1000;
defargs.halfsample = true;
defargs.overwrite = overwrite;
defargs = struct2varargin(defargs);        
pfl = MTAApfs(Trial,defargs{:});      



figure();
for unit = units,%pft.data.clu,
    clf();
    subplot2(3,6,1,1);    plot(pft,unit,1);
    subplot2(3,6,1,2);    plot(pft,unit,'mean');
    subplot2(3,6,1,3);    plot(pft,unit,'std');

    subplot2(3,6,2,1);    plot(pfh,unit,1);   
    subplot2(3,6,2,2);    plot(pfh,unit,'mean');
    subplot2(3,6,2,3);    plot(pfh,unit,'std');

    subplot2(3,6,3,1);    plot(pfl,unit,1);   
    subplot2(3,6,3,2);    plot(pfl,unit,'mean');
    subplot2(3,6,3,3);    plot(pfl,unit,'std');

    subplot2(3,6,1,4);    plot(pfkt,unit,1);
    subplot2(3,6,1,5);    plot(pfkt,unit,'mean');
    subplot2(3,6,1,6);    plot(pfkt,unit,'std');

    subplot2(3,6,2,4);    plot(pfkl,unit,1);   
    subplot2(3,6,2,5);    plot(pfkl,unit,'mean');
    subplot2(3,6,2,6);    plot(pfkl,unit,'std');

    subplot2(3,6,3,4);    plot(pfkp,unit,1);   
    subplot2(3,6,3,5);    plot(pfkp,unit,'mean');
    subplot2(3,6,3,6);    plot(pfkp,unit,'std');

    
    
    
        title(['unit: ',num2str(unit),' si:',num2str(pft.data.si(pft.data.clu==unit)),' spar:',num2str(pft.data.spar(pft.data.clu==unit))]);
    waitforbuttonpress();
end




% LOOP through sessions
sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta','theta-groom-sit'};

for t = 1:numel(Trials),
    clf();
    Trial = Trials{t};    
    pfs_2d_theta(Trial,'overwrite',true);
    units = select_placefields(Trial);
    if isempty(units), continue, end;
    
    for s = 1:numel(states),
        defargs = get_default_args('MjgER2016','MTAApfs','struct');
        defargs.units = units;
        defargs.states = states{s};
        defargs.overwrite = true;
        defargs = struct2varargin(defargs);
        MTAApfs(Trial,defargs{:});   
    end
    
    MjgER2016_drzfields(Trial,units,true);
    
end


