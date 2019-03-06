

MjgER2016_load_data();

t = 20;
Trial = Trials{t};
unitSubset = units{t};
sampleRate = 250;

% PREPROC xyz
xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);

pft = pfs_2d_theta(Trial,unitSubset);

hst = xyz(:,'hcom',:)-xyz(:,'spine_upper',:);

[th,phi] = cart2sph(hst(:,1,1),hst(:,1,2),hst(:,1,3));

thd = circ_dist(circshift(th,-5),circshift(th,5));
hrz = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);

spk = copy(Trial.spk);
spk.create(Trial,sampleRate,'theta-groom-sit',unitSubset,'deburst');

pfs = Trial;
for u = unitSubset,
    xyzp = copy(xyz);
    xyzp.data = [hrz(:,u==unitSubset),thd];
    xyzp.label = 'fetavhrz';
    pfs = MTAApfs(pfs,                ...
                  u,                  ... unit
                  'walk&theta',      ... state
                  true,               ... overwrite
                  'angvelhrzwalk',   ... tag
                  [0.1,0.0025],        ... binDims 
                  [1.8,1.8],          ... SmoothingWeights
                  'xy',               ... type 
                  false,              ... spkShuffle
                  false,              ... posShuffle
                  1,                  ... numIter
                  xyzp,               ... xyzp
                  [-1,1;-0.1,0.1],    ... boundaryLimits
                  'autoSaveFlag',false,...
                  'spk',spk);
    if u==unitSubset(1),pfs.save();end
end
pfs = MTAApfs(Trial,'tag','angvelhrzpause');

for u = unitSubset,
    mr = pft.maxRate(u);
    subplot(121);
    plot(pft,u,'mean','text',[mr],true);
    subplot(122);    
    plot(pfs,u,1,'text',[mr*2],false);
    title(num2str(u));
    waitforbuttonpress();
end
    