% req20180127 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Comparison of trb hcom and nose prediction for placefields
%  Bugs: NA




Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial.load('stc','msnn_ppsvd_raux');

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.binDims = [40,40];
defargs.SmoothingWeights = [1.25,1.25];
defargs.numIter = 1;
defargs.units = [];
defargs.tag = 'nosexl';
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs.halfsample = false;
defargs.trackingMarker = 'nose';

pfsArgs = struct2varargin(defargs);
pft = MTAApfs(Trial,pfsArgs{:});

defargs.tag = 'hcomxl';
defargs.trackingMarker = 'hcom';
pfsArgs = struct2varargin(defargs);
pfc = MTAApfs(Trial,pfsArgs{:});

defargs = get_default_args('MjgER2016','MTAApfs','struct');
defargs.numIter = 1;
defargs.units = [];
defargs.tag = '';
defargs.states = 'theta-groom-sit';
defargs.overwrite = true;
defargs.halfsample = false;
defargs.trackingMarker = 'hcom';
pfsArgs = struct2varargin(defargs);
pfth = MTAApfs(Trial,pfsArgs{:});


figure();
for u = pft.data.clu,
    clf();
    subplot(131);
    plot(pft,u);
    title(num2str(u));
    colorbar();
    subplot(132);    
    plot(pfc,u);
    colorbar();
    subplot(133);    
    plot(pfth,u);
    colorbar();
    waitforbuttonpress();
end

