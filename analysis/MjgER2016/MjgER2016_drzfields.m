function [pfd] =  MjgER2016_drzfields(Trial,varargin)
%  Status: active
%  Type: Analysis
%  Final_Forms: MjgER2016_drzfields.m
%  Description: directional rate zone (DRZ) versus various features
%  Bugs: NA
%  Original_Form: req20171101.m

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('units',                          [],                                           ...
                 'overwrite',                      false,                                        ...
                 'pfstats',                        [],                                           ...
                 'marker',                         'nose'                                        ...
);
[units,overwrite,pfstats,marker] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


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
pft = pfs_2d_theta(Trial);
[mrt,mrp] = pft.maxRate();


if overwrite,
    xyz = preproc_xyz(Trial,'trb');
% COMPUTE rhythmic head motion spectra
% COMPUTE rhythmic head motion power within 5-12 Hz frequency band
    [rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
    rhmp = rhm.copy();
    rhmp.data = median(log10(rhm(:,5<fs&fs<12)),2);
    rhmp.resample(xyz);
    % DIAGNOSTIC_FIG figure();  imagesc(ts,fs,log10(rhm.data)');axis('xy');colormap('jet');    
    % DIAGNOSTIC_FIG figure();  plot(rhmp.data);

% COMPUTE direction rate zones
    drz = compute_drz(Trial,units,pft,pfstats,[],marker);
% COMPUTE head and body pitch
    pch = fet_HB_pitchB(Trial,[],[],[],pitchReferenceTrial);    
end

pfd = {};

% $$$ % COMPUTE place fields
% $$$ pargs = get_default_args('MjgER2016','MTAApfs','struct');
% $$$ pargs.units      = units;
% $$$ pargs.states     = 'theta-groom-sit';
% $$$ pargs.numIter    = 1001;
% $$$ pargs.halfsample = true;
% $$$ 
% $$$ % DRZ X HEAD PITCH
% $$$ pargs.tag            = 'DRZxPITCH';
% $$$ pargs.boundaryLimits = [-1,1;-pi/2-0.3,pi/2+0.3];
% $$$ pargs.binDims        = [0.1,0.1];
% $$$ if overwrite,
% $$$ % LOAD pitches 
% $$$ % MAP pitches to reference trial
% $$$     pch = fet_HB_pitch(Trial);
% $$$     map_to_reference_session(pch,Trial,pitchReferenceTrial);    
% $$$     pargs.overwrite  = overwrite;    
% $$$     pargs.xyzp = MTADxyz('data',[drz(:,1),pch(:,3)],'sampleRate',xyz.sampleRate);
% $$$     for u = 1:numel(units);
% $$$         pargs.units      = units(u);        
% $$$         pargs.xyzp.data = [drz(:,u),pch(:,3)];
% $$$         pfsArgs = struct2varargin(pargs);
% $$$         pfs_dp = MTAApfs(Trial,pfsArgs{:});
% $$$     end
% $$$ end
% $$$ pargs.units     = units;
% $$$ pargs.overwrite = false;
% $$$ pfsArgs = struct2varargin(pargs);
% $$$ pfd{end} = MTAApfs(Trial,pfsArgs{:});



% DRZ X HEAD HEIGHT
pargs.tag              = 'DRZxHEIGHT_v2';%'DRZxHEIGHT';
pargs.boundaryLimits   = [-1,1;-10,350];
pargs.binDims          = [0.05,20];
pargs.numIter = 1001;
pargs.halfsample = true;
pargs.SmoothingWeights = [2,2]; %[1.5,1.5]
if overwrite,
    pargs.overwrite = true;
    pargs.xyzp = MTADxyz('data',[drz(:,1),xyz(:,'hcom',3)],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units      = units(u);
        pargs.xyzp.data = [drz(:,u),xyz(:,'head_front',3)];
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{end+1} = MTAApfs(Trial,pfsArgs{:});


% DRZ X RHMP HEIGHT
pargs.tag              = 'DRZxRHMP_v2';%'DRZxRHMP';
pargs.boundaryLimits   = [-1,1;-9,-2];
pargs.binDims          = [0.05,0.1];
pargs.numIter = 1001;
pargs.halfsample = true;
pargs.SmoothingWeights = [2,2]; %[1.5,1.5]
if overwrite,
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',[drz(:,1),rhmp.data],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units      = units(u);
        pargs.xyzp.data = [drz(:,u),rhmp.data];
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{end+1} = MTAApfs(Trial,pfsArgs{:});


pargs.tag            = 'DRZxbp_v2';%'DRZxbp';
pargs.boundaryLimits = [-1,1;-pi/2-0.2,pi/2+0.2];
pargs.binDims        = [0.05,0.05];%[0.1,0.1]
pargs.numIter        = 1001;
pargs.halfsample     = true;
pargs.SmoothingWeights = [2,2]; %[1.5,1.5]
if overwrite,
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',[drz(:,1),pch(:,1)],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units      = units(u);
        pargs.xyzp.data = [drz(:,u),pch(:,1)];
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{end+1} = MTAApfs(Trial,pfsArgs{:});


pargs.tag            = 'DRZxhbp_v2';%'DRZxhbp';
pargs.boundaryLimits = [-1,1;-pi/2-0.2,pi/2+0.2];
pargs.binDims        = [0.05,0.05];
pargs.numIter        = 1001;
pargs.halfsample     = true;
pargs.SmoothingWeights = [2,2];
if overwrite,
    pargs.overwrite = true;    
    pargs.xyzp = MTADxyz('data',[drz(:,1),pch(:,2)],'sampleRate',xyz.sampleRate);
    for u = 1:numel(units);
        pargs.units     = units(u);
        pargs.xyzp.data = [drz(:,u),pch(:,2)];
        pfsArgs = struct2varargin(pargs);
        MTAApfs(Trial,pfsArgs{:});
    end
end
pargs.units     = units;
pargs.overwrite = false;
pfsArgs = struct2varargin(pargs);
pfd{end+1} = MTAApfs(Trial,pfsArgs{:});
