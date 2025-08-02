% req20220128
% rhm theta phase difference x gdz rate map
% results meh

MjgER2016_load_data();
Trial = Trials{20};

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');
stc = Trial.load('stc');
ang = create(MTADang,Trial,xyz);

xts = [1:size(xyz,1)]./xyz.sampleRate;

rbm = copy(xyz);
rbm.data = nunity(minus(ang(:,'pelvis_root','spine_upper',3),ang(:,'pelvis_root','spine_middle',3)));
rbm.filter('ButFilter',4,[0.8,16],'bandpass');

rbm = copy(xyz);
rbm.data = nunity(MedianFilter(minus(ang(:,'pelvis_root','spine_upper',3),ang(:,'pelvis_root','spine_middle',3)),300));
rbm.filter('ButFilter',4,[0.8,4],'bandpass');

figure,
subplot(211);
hold('on');
plot(xts,rbm.data*20)
%plot(xts,rpm.data*50)
%plot(xts,nunity(rhm.data)*10)
subplot(212);
plotSTC(stc,1,'staggeredStates',false);
linkaxes(findobj(gcf(),'Type','Axes'),'x');



phz = load_theta_phase(Trial,rhm,69,sessionList(20).subject.correction.thetaPhase);


unitSubset = units{20};
spk = Trial.load('spk',250,'theta-groom-sit',unitSubset);
spk.clu(phz(spk.res)<pi) = [];
spk.res(phz(spk.res)<pi) = [];

pft = pfs_2d_theta(Trial,unitSubset);
ghz = compute_ghz(Trial,unitSubset,pft,[],[],'hcom',[],[],rhm);


for u = 1:numel(unitSubset);
    pfs.purge_savefile();
fet = copy(rhm);
fet.data = cat(2,ghz(:,u),circ_dist(rhmPhz.data,phz.data));

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = unitSubset(u);
pargs.tag          = 'ghz_rhmtheta';
pargs.binDims      = [ 0.1,0.2]; % X Y 
pargs.SmoothingWeights = [2, 2];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-1,1;-pi,pi];
pargs.states       = 'lloc+lpause&theta';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
pargs.xyzp = copy(fet);
pargs.spk = copy(spk);
pargs.compute_pfs = @PlotPFCirc;
pfsArgs = struct2varargin(pargs);
pfs = MTAApfs(Trial,pfsArgs{:});


clf(),
subplot(131);plot(pft,unitSubset(u),1,'colorbar',[],true,'colorMap',@jet);
subplot(132);plot(pfs,unitSubset(u),1,'colorbar',[],false,'colorMap',@jet);
title(unitSubset(u));

    pfs.purge_savefile();
fet = copy(rhm);
fet.data = cat(2,ghz(:,u),circ_dist(rhmPhz.data,circshift(phz.data,50100)));


pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = unitSubset(u);
pargs.tag          = 'ghz_rhmtheta';
pargs.binDims      = [ 0.1,0.2]; % X Y 
pargs.SmoothingWeights = [2, 2];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-1,1;-pi,pi];
pargs.states       = 'lloc+lpause&theta';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
pargs.xyzp = copy(fet);
pargs.spk = copy(spk);
pargs.compute_pfs = @PlotPFCirc;
pfsArgs = struct2varargin(pargs);
pfs = MTAApfs(Trial,pfsArgs{:});

subplot(133);plot(pfs,unitSubset(u),1,'colorbar',[],false,'colorMap',@jet);
waitforbuttonpress();
end