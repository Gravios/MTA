% Update and test the MTAAknnpfs class for the original placefield
% thresholding for significant patch detection.

Trial = MTATrial.validate('Ed10-20140817.cof.all');
Trial.load('stc','NN0317R');

defargs = get_default_args('MjgEdER2016','MTAAknnpfs','struct');
defargs.units = units;
defargs.states = 'loc&theta';
defargs = struct2varargin(defargs);        
pfs = MTAAknnpfs(Trial,defargs{:});      
