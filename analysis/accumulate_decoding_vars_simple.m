function dc = accumulate_decoding_vars_simple(Trial,units,varargin)

global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(     'sampleRate', 250,                                                        ...
                 'halfSpikeWindow', 0.022,                                                      ...
                       'ufrWindow', 0.008,                                                       ...
                'smoothingWeights', [800.^2,800.^2]                                             ...
);
[sampleRate,halfSpikeWindow,ufrWindow,smoothingWeights] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

phzBins = linspace(0,2*pi,25);
states = {'theta','vel','gper'};
xyz = preproc_xyz(Trial,'trb',sampleRate);

switch Trial.maze.name
    case 'sof'
% PLACEFIELDS 
%pfs = compute_xyhb_ratemaps( Trial, units);
% compute_ratemaps ---------------------------------------------------------------------------------
global AP
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'vel&gper&theta',             ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-1000,1000;-1000,1000],      ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------
    case 'lin'
% PLACEFIELDS 
%pfs = compute_xyhb_ratemaps( Trial, units);
% compute_ratemaps ---------------------------------------------------------------------------------
global AP
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'vel&gper&theta',             ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-1600,1600;-300,300],      ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------
end
pfs = compute_ratemaps( Trial, units,'overwrite',false);
%ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
%mask = ds.mask;
switch Trial.maze.name
  case 'lin'
    mask = create_tensor_mask(pfs.adata.bins,struct('shape','line',    'edgeLength',2800));
  case 'sof'
    mask = create_tensor_mask(pfs.adata.bins,struct('shape','square',  'edgeLength',1600));
  case 'cof'    
    mask = create_tensor_mask(pfs.adata.bins,struct('shape','circular','radius',     440)); 
  otherwise
    error('MTA:analysis:accumulate_decoding_vars_simple: Maze Type Not Found');
end


% UNITS 
ufr = Trial.load('ufr', ...
                 xyz,...
                 Trial.load('spk', sampleRate, '', units, ''),...
                 units,...
                 ufrWindow,...
                 'boxcar',...
                 true);

% DECODE 
dc = decode_ufr_boxcar(Trial,          ...
                       units,          ...
                       sampleRate,     ...
                       ufr,            ...
                       pfs,            ...
                       mask,           ...
                       halfSpikeWindow,...
                       smoothingWeights,...
                       'overwrite',false);
% $$$ load('/storage/gravio/data/project/general/ec14-00271287/ec14-00271287.sof.all.decode_ufr_boxcar.dfe8952bbcedcdb67cc5bc38040c1dee.mat','dc');

dc.hRot = -Trial.meta.correction.headYaw;

% THETA PHASE at xyz sampling rate
dc.phz = load_theta_phase(Trial, xyz);
dc.phz = dc.phz(dc.ind,1);
dc.iphz = discretize(dc.phz,phzBins);

dc.pfs = pfs;

dc.xyz  = filter(copy(xyz),'ButFilter',4,20,'low');;
dc.xyz.data  = xyz.data(dc.ind,:,:);

% SPEED head and body
dc.vxy  = filter(copy(xyz),'ButFilter',4,1.5,'low');;
dc.vxy  = dc.vxy.vel('hcom',[1,2]);
dc.vxy.data  = dc.vxy(dc.ind,1);
dc.lvxy = copy(dc.vxy);
dc.lvxy.data(dc.lvxy.data<=0) = 0.001;
dc.lvxy.data = log10(dc.lvxy.data);

% HEAD VECTOR
dc.hvec = dc.xyz(:,'nose',[1,2])-dc.xyz(:,'hcom',[1,2]);
dc.hvec = sq(bsxfun(@rdivide,dc.hvec,sqrt(sum(dc.hvec.^2,3))));
dc.hvec = cat(3,dc.hvec,sq(dc.hvec)*[0,-1;1,0]);
dc.hvec = multiprod(dc.hvec,                           ...
                    [cos(dc.hRot),-sin(dc.hRot); ...
                    sin(dc.hRot), cos(dc.hRot)],...
                    [2,3],                                   ...
                    [1,2]);

% TRAJECTORY VECTOR
dc.tvec = circshift(xyz(:,'hcom',[1,2]),-round(250*0.1)) - ...
          circshift(xyz(:,'hcom',[1,2]),round(250*0.1));
dc.tvec = sq(bsxfun(@rdivide,dc.tvec,sqrt(sum(dc.tvec.^2,3))));
dc.tvec = cat(3,dc.tvec,sq(dc.tvec)*[0,-1;1,0]);
dc.tvec = dc.tvec(dc.ind,:,:);

% POSTURE pitch
% $$$ dc.fet  = fet_HB_pitchB(Trial,sampleRate);
% $$$ dc.fet.data  = dc.fet.data(dc.ind,:);

% DISTANCE to maze center 
dc.hdist = filter(copy(xyz),'ButFilter',3,5,'low');
xycoor =    xyz(:,'hcom',[1,2]);
[~,dc.hdist.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
dc.hdist.data = dc.hdist.data(dc.ind);


% COM 
dc.ecom = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.com(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.hvec(:,:,:),2,[2,3]));
% SAX 
dc.esax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.sax(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.hvec(:,:,:),2,[2,3]));
% MAX 
dc.emax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.max(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.hvec(:,:,:),2,[2,3]));
% LOM 
dc.elom = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.lom(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.hvec(:,:,:),2,[2,3]));

% LAX 
dc.elax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.lax(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.hvec(:,:,:),2,[2,3]));



% COM 
dc.tcom = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.com(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.tvec(:,:,:),2,[2,3]));
% SAX 
dc.tsax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.sax(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.tvec(:,:,:),2,[2,3]));
% MAX 
dc.tmax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.max(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.tvec(:,:,:),2,[2,3]));
% LOM 
dc.tlom = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.lom(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.tvec(:,:,:),2,[2,3]));

% LAX 
dc.tlax = sq(multiprod(permute(bsxfun(@minus,...
                                      dc.lax(:,[1,2],:),...
                                      sq(dc.xyz(:,'hcom',[1,2]))),...
                               [1,2,4,3]),...
                       dc.tvec(:,:,:),2,[2,3]));


% ANGULAR velocity head
% $$$ dc.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
% $$$ xycoor = cat(2,...
% $$$              dc.hvang(:,'spine_upper',[1,2])-dc.hvang(:,'bcom',[1,2]),...
% $$$              dc.hvang(:,'nose',[1,2])-dc.hvang(:,'hcom',[1,2]));
% $$$ dc.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% $$$ % Positive: CCW (Left)     Negative: CW (Right)
% $$$ dc.hvang.data = circ_dist(circshift(dc.hvang.data(:,2),-10),...
% $$$                           circshift(dc.hvang.data(:,2),+10));
% $$$ dc.hvang.data = dc.hvang.data(dc.ind);


% COMPUTE corrected hbang
% $$$ dc.hbang = filter(copy(xyz),'ButFilter',3,14,'low');    
% $$$ xycoor = cat(2,...
% $$$              dc.hbang(:,'spine_upper',[1,2])-dc.hbang(:,'bcom',[1,2]),...
% $$$              dc.hbang(:,'nose',[1,2])-dc.hbang(:,'hcom',[1,2]));
% $$$ dc.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% $$$ dc.hbang.data = circ_dist(dc.hbang.data(:,2),dc.hbang.data(:,1));
% $$$ dc.hbang.data = dc.hbang.data(dc.ind);
% $$$ % hbang correction set in metadata
% $$$ dc.hbang.data = dc.hbang.data + Trial.meta.correction.headBody;


% GENERATE state collection matrix    
% $$$ dc.stcm = stc2mat(Trial.stc,xyz,states);
% $$$ dc.stcm = dc.stcm(dc.ind,:);

% COMPUTE head and body relative velocity vectors
hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
% $$$ bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
dc.hvfl = hvfl(dc.ind,:);
% $$$ dc.bvfl = bvfl(dc.ind,:);    

dc.states = states;