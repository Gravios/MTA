% req20171101 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: MjgER2016_drzfields.m
%  Description: directional rate zone (DRZ) versus various features
%  Bugs: NA

sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);

states  = {'loc','lloc','hloc','rear','pause','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states{end+1} = 'theta-groom-sit';


t = 1;
Trial = MTATrial.validate(sessionList(t));
stc = label_bhv_reduced(Trial.load('stc','msnn_ppsvd'),Trial);
stcMode = 'msnn_ppsvd_raux';    
Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);


% LOAD theta state placefields
pft = pfs_2d_theta(Trial);

pfstats = batch_compute_pfstats_bs(sessionListName);

units = pfstats{t}.cluMap;
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong');
% DIAGNOSTIC_FIG figure,imagesc(ts,fs,log10(rhm.data)');axis('xy');colormap('jet');
rhmp = rhm.copy();
rhmp.data = median(log10(rhm(:,5<fs&fs<12)),2)
rhmp.resample(xyz);
% DIAGNOSTIC_FIG figure,plot(rhmp.data)

% Get the mean firing rate for each xy position along trajectory 
wpmr{t} = zeros(xyz.size(1),numel(units));
[~,indx] = min(abs( repmat(pft.adata.bins{1}',xyz.size(1),1)...
                    -repmat(xyz(:,Trial.trackingMarker,1),1,numel(pft.adata.bins{1}))),...
               [],2);
[~,indy] = min(abs( repmat(pft.adata.bins{2}',xyz.size(1),1)...
                    -repmat(xyz(:,Trial.trackingMarker,2),1,numel(pft.adata.bins{2}))),...
               [],2);
rateMapIndex = sub2ind(pft.adata.binSizes',indx,indy);
for unit = units,
    rateMap = pft.plot(unit,'mean');      %  for MTAApfs
                                             %rateMap = rot90(rot90(rateMap)'); % for MTAAknnpfs
    wpmr{t}(:,unit==units) = rateMap(rateMapIndex);
end


% Get the rat's heading 
pfds = [];
pfdd = [];

for unit = units
    pfhxy = xyz(:,{'head_back','head_front'},:);
    pfhxy = cat(2,pfhxy,permute(repmat([fliplr(sq(pfstats{t}.peakPatchCOM(8,1,unit==units,:))'),0],...
                                       [size(xyz,1),1]),...
                                [1,3,2]));
    pfhxy = MTADxyz([],[],pfhxy,xyz.sampleRate);
    
    cor = cell(1,3);
    [cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
    cor = cell2mat(cor);
    
    por = cell(1,3);
    [por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
    por = cell2mat(por);
    
    pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));
    pfdd(:,unit==units) = por(:,3);
    
end

pfdist{t} = pfdd;
pfd{t} = zeros(size(pfds));
pfd{t}(abs(pfds)<=pi/2)=-1;
pfd{t}(abs(pfds)>pi/2)=1;

% Calculate DRZ 
DRZ{t} = pfd{t}.*(1-wpmr{t}./mean(pfstats{t}.peakPatchRate(8,:,1)));


%% PLACEFIELD DRZ X HEAD PITCH --------------------------------------------------------------
for u = 1:numel(units);
    defargs = get_default_args('MjgER2016','MTAApfs','struct');
    defargs.units     = units(u);
    defargs.states    = 'theta';
    defargs.binDims   = [0.1,0.1];
    defargs.boundaryLimits = [-1,1;-pi/2,pi/2];    
    defargs.xyzp      = MTADxyz('data',[DRZ{t}(:,u),ang(:,'head_back','head_front',2)],...
                                'sampleRate',xyz.sampleRate);
    defargs = struct2varargin(defargs);        
    pfs_dp = MTAApfs(Trial,defargs{:});      
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units     = units;
defargs.states    = 'theta';
defargs.binDims   = [0.1,0.1];
defargs.boundaryLimits = [-1,1;-pi/2,pi/2];
defargs = struct2varargin(defargs);        
pfs_dp = MTAApfs(Trial,defargs{:});      



%% PLACEFIELD DRZ X HEAD HEIGHT -------------------------------------------------------------
for u = 1:numel(units);
    defargs = get_default_args('MjgER2016','MTAApfs','struct');
    defargs.units     = units(u);
    defargs.states    = 'theta';
    defargs.binDims   = [0.1,20];
    defargs.boundaryLimits = [-1,1;0,300];    
    defargs.xyzp      = MTADxyz('data',[DRZ{t}(:,u),xyz(:,'head_front',3)],...
                                'sampleRate',xyz.sampleRate);
    defargs = struct2varargin(defargs);        
    pfs_dh = MTAApfs(Trial,defargs{:});      
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units     = units;
defargs.states    = 'theta';
defargs.binDims   = [0.1,20];
defargs.boundaryLimits = [-1,1;0,300];
defargs = struct2varargin(defargs);        
pfs_dh = MTAApfs(Trial,defargs{:});      


%% PLACEFIELD DRZ X RHMP HEIGHT -------------------------------------------------------------
for u = 1:numel(units);
    defargs = get_default_args('MjgER2016','MTAApfs','struct');
    defargs.units     = units(u);
    defargs.states    = 'theta';
    defargs.binDims   = [0.1,0.2];
    defargs.boundaryLimits = [-1,1;-8.5,-4.5];    
    defargs.xyzp      = MTADxyz('data',[DRZ{t}(:,u),rhmp.data],...
                                'sampleRate',xyz.sampleRate);
    defargs = struct2varargin(defargs);        
    pfs_dm = MTAApfs(Trial,defargs{:});      
end

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.units     = units;
defargs.states    = 'theta';
defargs.binDims   = [0.1,20];
defargs.boundaryLimits = [-1,1;0,300];
defargs = struct2varargin(defargs);        
pfs_dm = MTAApfs(Trial,defargs{:});      



unit = units(2);
figure();
for unit = units,
    clf();
    subplot(131); % placefield theta
    pft.plot(unit);
    subplot(132); % drz X pitch 
    pfs_dp.plot(unit,'isCircular',false);
    colorbar();
    subplot(133); % drz X height
    pfs_dh.plot(unit,'isCircular',false);
    colorbar();    
    waitforbuttonpress();
end