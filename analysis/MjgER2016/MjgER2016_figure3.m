

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

FigDir = '/storage/gravio/figures/placefields';
mkdir(FigDir);

 
% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',    ...
          'pause&theta','lpause&theta','hpause&theta',           ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear',    ...
          'pause','lpause','hpause',           ...
          'theta-groom-sit'};

overwrite = false;
testRequested = false;

% pfstats = cf(@(t)  compute_pfstats_bs(t),       Trials);
% $$$ cf(@(t)  MjgER2016_drzfields(t,true), Trials);

% FOR each Trial -------------------------------------------------------------

for t = 1:20,
    Trial = Trials{t};
    mkdir(fullfile(FigDir,Trial.filebase));
    stc = Trial.stc.copy();
    %spk = Trial.spk.copy();
    %spk.load_spk(Trial);
    xyz = Trial.load('xyz');
    %pch = fet_HB_pitch(Trial);
    %map_to_reference_session(pch,Trial,pitchReferenceTrial);    
    
    pfstats = compute_pfstats_bs(Trial,'overwrite',overwrite);
    units = pfstats.cluMap;
    pft = pfs_2d_theta(Trial);
    
    
% LOAD placefields and subsampled estimate
    for sts = 1:numel(states),        
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = pfstats.cluMap;
        %defargs.overwrite = true;
        defargs.states = states{sts};
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end

    unit = 150;
    s = 1;
    o = 1;
    rmap = plot(pfkbs{o},unit);
    patchInd = pfstats.pfmstats(s,unit==pfstats.cluMap).patchRateInd(1,1,:,:);
    oPatchRates = rmap(patchInd(1,:),patchInd(2,:));
    
    

    
    figure();
    if testRequested,
        for u = 1:numel(units),
            subplot(121);
            plot(pft,units(u));
            subplot(122);
            plot(pfkbs{8},units(u));
            waitforbuttonpress();
        end
    end
    
    
% LOAD DRZ fields
    dfs = cell([1,3]);
    [dfs{:}] = MjgER2016_drzfields(Trial,true);
    dfst = {'pitch','height','rhm'};

% COMPUTE place fields and subsampled estimate
% $$$ for sts = 1:numel(states),
% $$$         [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5);
% $$$     end
% $$$     sper{end}{1}(sper{end}{1}>(size(xyz,1)./xyz.sampleRate.*1250)) = [];
% $$$     sper{end}{2}(sper{end}{2}>(size(xyz,1)./xyz.sampleRate.*1250)) = [];    


    %[accg,tbins] = autoccg(Trial);

    % PLACEFIELD figure now located in placefield_summary.m
end