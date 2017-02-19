% Update and test the MTAAknnpfs class for the original placefield
% thresholding for significant patch detection.

Trial = MTATrial.validate('Ed10-20140817.cof.all');
Trial.load('stc','NN0317R');

% SELECT units based on Isolation Distance and Theta place field activity
pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;
units = select_units(Trial,18);
units = units(mrt(pft.data.clu(units))>1);

defargs = get_default_args('MjgEdER2016','MTAAknnpfs','struct');
defargs.units = units;
defargs.states = 'loc&theta';
defargs.overwrite = false;
defargs = struct2varargin(defargs);        
pfs = MTAAknnpfs(Trial,defargs{:});      

figure
unit = 4;
clf
subplot(131);
pfs.plot(unit);
colorbar;
subplot(132);
pfs.plot(unit,'mean');
colorbar;
subplot(133);
pfs.plot(unit,'std');
colorbar;


mmap = pfs.plot(unit,'mean');
smap = pfs.plot(unit,'std');
nanmean(mmap(:))+nanstd(smap(:))*5;