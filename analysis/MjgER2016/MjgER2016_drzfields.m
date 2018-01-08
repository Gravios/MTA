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


pitchReferenceTrial = 'Ed05-20140529.ont.all';

states  = {'loc&theta','lloc&theta','hloc&theta','rear&theta','pause&theta','lpause&theta','hpause&theta'};
states{end+1} = 'theta-groom-sit';

% LOAD behavioral state collection
stcMode = 'msnn_ppsvd_raux';    
Trial.load('stc',stcMode);


% LOAD placefield statistics 
% $$$ if isempty(pfstats),  pfstats = compute_pfstats_bs(Trial);    end
%if isempty(units),    units   = pfstats.clu;                  end
if isempty(units),    units   = select_placefields(Trial);                  end

% LOAD theta state placefields
% $$$ defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
% $$$ defargs.units = pfstats.cluMap;
% $$$ defargs.states = 'theta-groom-sit';
% $$$ defargs = struct2varargin(defargs);        
% $$$ pft = MTAAknnpfs_bs(Trial,defargs{:});      
% $$$ [mrt,mrp] = pft.maxRate([],'mean');
pft = pfs_2d_theta(Trial);
[mrt,mrp] = pft.maxRate();


if overwrite,
    xyz = Trial.load('xyz');
% COMPUTE rhythmic head motion spectra
% COMPUTE rhythmic head motion power within 5-12 Hz frequency band
    [rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
    rhmp = rhm.copy();
    rhmp.data = median(log10(rhm(:,5<fs&fs<12)),2);
    rhmp.resample(xyz);
    % DIAGNOSTIC_FIG figure();  imagesc(ts,fs,log10(rhm.data)');axis('xy');colormap('jet');    
    % DIAGNOSTIC_FIG figure();  plot(rhmp.data);

% COMPUTE direction rate zones
    drz = compute_drz(Trial,units,pft);
end



% COMPUTE place fields
pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units      = units;
pargs.states     = 'theta-groom-sit';
pargs.overwrite  = overwrite;
pargs.numIter    = 1001;
pargs.halfsample = true;

% DRZ X HEAD PITCH
pargs.tag            = 'DRZxPITCH';
pargs.boundaryLimits = [-1,1;-pi/2-0.3,pi/2+0.3];
pargs.binDims        = [0.1,0.1];
if overwrite,
% LOAD pitches 
% MAP pitches to reference trial
    pargs.xyzp = MTADxyz('data',[drz(:,1),pch(:,3)],'sampleRate',xyz.sampleRate);
    pch = fet_HB_pitch(Trial);
    map_to_reference_session(pch,Trial,pitchReferenceTrial);
    for u = 1:numel(units);
        pargs.units      = units(u);        
        pargs.xyzp.data = [drz(:,u),pch(:,3)];
        pfsArgs = struct2varargin(pargs);
        pfs_dp = MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units      = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfs_dp = MTAApfs(Trial,pfsArgs{:});



% DRZ X HEAD HEIGHT
pargs.tag              = 'DRZxHEIGHT';    
pargs.boundaryLimits   = [-1,1;-10,350];
pargs.binDims          = [0.1,20];
if overwrite,
    pargs.xyzp = MTADxyz('data',[drz(:,1),xyz(:,'head_front',3)],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units      = units(u);
        pargs.xyzp.data = [drz(:,u),xyz(:,'head_front',3)];
        pfsArgs = struct2varargin(pargs);
        pfs_dh = MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units      = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfs_dh = MTAApfs(Trial,pfsArgs{:});



% DRZ X RHMP HEIGHT
pargs.tag              = 'DRZxRHMP';        
pargs.boundaryLimits   = [-1,1;-9,-2];
pargs.binDims          = [0.1,0.15];
if overwrite,
    pargs.xyzp = MTADxyz('data',[drz(:,1),rhmp.data],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units      = units(u);
        pargs.xyzp.data = [drz(:,u),rhmp.data];
        pfsArgs = struct2varargin(pargs);
        pfs_dm = MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units      = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfs_dm = MTAApfs(Trial,pfsArgs{:});

