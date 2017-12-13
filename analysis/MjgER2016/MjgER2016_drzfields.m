function [pfs_dp,pfs_dh,pfs_dm] =  MjgER2016_drzfields(Trial,varargin)
%  Status: active
%  Type: Analysis
%  Final_Forms: MjgER2016_drzfields.m
%  Description: directional rate zone (DRZ) versus various features
%  Bugs: NA
%  Original_Form: req20171101.m

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                          [],                                           ...
                 'overwrite',                      false,                                        ...
                 'pfstats',                        []                                            ...
);
[units,overwrite,pfstats] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

testRequested = false;

pitchReferenceTrial = 'Ed05-20140529.ont.all';

states  = {'loc&theta','lloc&theta','hloc&theta','rear&theta','pause&theta','lpause&theta','hpause&theta'};
states{end+1} = 'theta-groom-sit';

% LOAD behavioral state collection
stcMode = 'msnn_ppsvd_raux';    
Trial.load('stc',stcMode);


% LOAD placefield statistics 
if isempty(pfstats),  pfstats = compute_pfstats_bs(Trial);    end
if isempty(units),    units   = pfstats.clu;                  end


% LOAD theta state placefields
% $$$ defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
% $$$ defargs.units = pfstats.cluMap;
% $$$ defargs.states = 'theta-groom-sit';
% $$$ defargs = struct2varargin(defargs);        
% $$$ pft = MTAAknnpfs_bs(Trial,defargs{:});      
% $$$ [mrt,mrp] = pft.maxRate([],'mean');
pft = pfs_2d_theta(Trial);
[mrt,mrp] = pft.maxRate();




units = pfstats.cluMap;

if overwrite,
    xyz = Trial.load('xyz');
    [rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
    % DIAGNOSTIC_FIG figure,imagesc(ts,fs,log10(rhm.data)');axis('xy');colormap('jet');
    rhmp = rhm.copy();
    rhmp.data = median(log10(rhm(:,5<fs&fs<12)),2);
    rhmp.resample(xyz);
    % DIAGNOSTIC_FIG figure,plot(rhmp.data)

    % COMPUTE direction rate zones
    drz = compute_drz(Trial,pft,units);
end

% LOAD pitches 
% MAP pitches to reference trial
pch = fet_HB_pitch(Trial);
map_to_reference_session(pch,Trial,pitchReferenceTrial);    


%% PLACEFIELD DRZ X HEAD PITCH --------------------------------------------------------------
if overwrite,
    for u = 1:numel(units);
        defargs = get_default_args('MjgER2016','MTAApfs','struct');
        defargs.tag       = 'drzXpitch';
        defargs.units     = units(u);
        defargs.states    = 'theta-groom-sit';
        defargs.overwrite = overwrite;    
        defargs.binDims   = [0.1,0.1];
        defargs.boundaryLimits = [-1,1;-pi/2,pi/2];    
        defargs.xyzp      = MTADxyz('data',[drz(:,u),pch(:,3)],...
                                    'sampleRate',xyz.sampleRate);
        defargs = struct2varargin(defargs);        
        pfs_dp = MTAApfs(Trial,defargs{:});      
    end
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.tag       = 'drzXpitch';
defargs.units     = units;
defargs.states    = 'theta-groom-sit';
defargs.binDims   = [0.1,0.1];
defargs.boundaryLimits = [-1,1;-pi/2,pi/2];
defargs = struct2varargin(defargs);        
pfs_dp = MTAApfs(Trial,defargs{:});      

if testRequested,
    figure();
    for u = 1:numel(units),
        subplot(121);  hold('on'); 
        plot(pft,units(u)); 
        plot(mrp(units(u)==pft.data.clu,1),mrp(units(u)==pft.data.clu,2),'*m');
        subplot(122);  plot(pfs_dp,units(u),'isCircular',false);
        colorbar();
        waitforbuttonpress();
    end
end


%% PLACEFIELD DRZ X HEAD HEIGHT -------------------------------------------------------------
if overwrite,
    for u = 1:numel(units);
        defargs = get_default_args('MjgER2016','MTAApfs','struct');
        defargs.tag       = 'drzXheight';    
        defargs.units     = units(u);
        defargs.states    = 'theta-groom-sit';
        defargs.overwrite = overwrite;    
        defargs.binDims   = [0.1,20];
        defargs.boundaryLimits = [-1,1;0,300];    
        defargs.xyzp      = MTADxyz('data',[drz(:,u),xyz(:,'head_front',3)],...
                                    'sampleRate',xyz.sampleRate);
        defargs = struct2varargin(defargs);        
        pfs_dh = MTAApfs(Trial,defargs{:});      
    end
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.tag       = 'drzXheight';    
defargs.units     = units;
defargs.states    = 'theta-groom-sit';
defargs.binDims   = [0.1,20];
defargs.boundaryLimits = [-1,1;0,350];
defargs = struct2varargin(defargs);        
pfs_dh = MTAApfs(Trial,defargs{:});      

if testRequested,
    figure();
    for u = 1:numel(units),
        subplot(121);  hold('on'); 
        plot(pft,units(u)); 
        plot(mrp(units(u)==pft.data.clu,1),mrp(units(u)==pft.data.clu,2),'*m');
        subplot(122);  plot(pfs_dh,units(u),'isCircular',false);
        colorbar();
        waitforbuttonpress();
    end
end

%% PLACEFIELD DRZ X RHMP HEIGHT -------------------------------------------------------------
if overwrite,
    for u = 1:numel(units);
        defargs = get_default_args('MjgER2016','MTAApfs','struct');
        defargs.tag       = 'drzXrhmp';        
        defargs.units     = units(u);
        defargs.states    = 'theta-groom-sit';
        defargs.binDims   = [0.1,0.15];
        defargs.overwrite = overwrite;
        defargs.boundaryLimits = [-1,1;-8.5,-3];    
        defargs.xyzp      = MTADxyz('data',[drz(:,u),rhmp.data],...
                                    'sampleRate',xyz.sampleRate);
        defargs = struct2varargin(defargs);        
        pfs_dm = MTAApfs(Trial,defargs{:});      
    end
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.tag       = 'drzXrhmp';        
defargs.units     = units;
defargs.states    = 'theta-groom-sit';
defargs.binDims   = [0.1,0.15];
defargs.boundaryLimits = [-1,1;-8.5,-3];    
defargs = struct2varargin(defargs);        
pfs_dm = MTAApfs(Trial,defargs{:});      

if testRequested,
    figure();
    for u = 1:numel(units),
        subplot(121);  hold('on'); 
        plot(pft,units(u)); 
        plot(mrp(units(u)==pft.data.clu,1),mrp(units(u)==pft.data.clu,2),'*m');
        subplot(122);  plot(pfs_dm,units(u),'isCircular',false);
        colorbar();
        waitforbuttonpress();
    end
end


% $$$ figure();
% $$$ for unit = units,
% $$$     clf();
% $$$     subplot(141); % placefield theta
% $$$     hold('on');
% $$$     pft.plot(unit);
% $$$     plot(mrp(pft.data.clu==unit,1),mrp(pft.data.clu==unit,2),'*m');
% $$$     plot(pfstats{1}.peakPatchCOM(8,1,unit==units,2),pfstats{1}.peakPatchCOM(8,1,unit==units,1),'*g')
% $$$     subplot(142); % drz X pitch 
% $$$     pfs_dp.plot(unit,'isCircular',false);
% $$$     colorbar();
% $$$     subplot(143); % drz X pitch 
% $$$     pfs_dm.plot(unit,'isCircular',false);
% $$$     colorbar();
% $$$     subplot(144); % drz X height
% $$$     pfs_dh.plot(unit,'isCircular',false);
% $$$     colorbar();    
% $$$     waitforbuttonpress();
% $$$ end