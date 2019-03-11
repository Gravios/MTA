

MjgER2016_load_data();

t = 20;
Trial = Trials{t};
unitSubset = units{t};
sampleRate = 250;
stc = Trial.stc.copy();

pft = pfs_2d_theta(Trial,unitSubset);

xyz = preproc_xyz(Trial,'trb');

figure,PlotSessionErrors(Trial,gcf,xyz)
figure,PlotSessionErrors(Trial,gcf)


s = MTASession.validate(Trial.filebase)


% PREPROC xyz

xyz.resample(sampleRate);


try, lfp = Trial.load('lfp',sessionList(t).thetaRef);
catch, lfp = Trial.load('lfp',sessionList(t).thetaRef);
end
phz = lfp.phase([5,13]);    
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);


roll = fet_roll(Trial,sampleRate);

stcm = stc2mat(stc,xyz,states);


fxyz = filter(copy(xyz),'ButFilter',3,2,'low');

hst = xyz(:,'head_front',:)-xyz(:,'head_back',:);
[thHead,phiHead] = cart2sph(hst(:,1,1),hst(:,1,2),hst(:,1,3));

fhst = fxyz(:,'head_front',:)-fxyz(:,'head_back',:);
[fthHead,fphiHead] = cart2sph(fhst(:,1,1),fhst(:,1,2),fhst(:,1,3));


%figure,plot((1:size(xyz,1))./xyz.sampleRate,circ_dist(fthHead,thHead))

hsway = copy(xyz);
hsway.data = circ_dist(fthHead,thHead);
fhsway = hsway.copy();
fhsway.filter('ButFilter',3,4,'low');

dhsway = copy(xyz);
dhsway.data = circshift(fhsway.data,-1)-circshift(fhsway.data,1);


% COMPUTE head frame of reference vectors for all time points
fxyz = filter(copy(xyz),'ButFilter',3,30,'low');
hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);

tvec = circshift(fxyz(:,'hcom',[1,2]),-10)-circshift(fxyz(:,'hcom',[1,2]),10);
tvec = sq(bsxfun(@rdivide,tvec,sqrt(sum(tvec.^2,3))));
tvec = cat(3,tvec,sq(tvec)*[0,-1;1,0]);


pfhr = nan([size(xyz,1),numel(unitSubset),2]);
pfdr = nan([size(xyz,1),numel(unitSubset),2]);
pfhrTh = nan([size(xyz,1),numel(unitSubset),1]);
pfdrTh = nan([size(xyz,1),numel(unitSubset),1]);
pfdha = nan([size(xyz,1),numel(unitSubset),1]);
for u = 1:numel(unitSubset),%&pft.data.spar>0.15&pft.data.spar<0.3),
    [mxr,mxp] = pft.maxRate(unitSubset(u));
    pfhr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),hvec,2,[2,3]);
    pfdr(:,u,:) = multiprod(bsxfun(@minus,mxp,sq(fxyz(:,'hcom',[1,2]))),tvec,2,[2,3]);    

    [pfhrTh(:,u),pfhrRd(:,u)] = cart2pol(pfhr(:,u,1),pfhr(:,u,2));
    [pfdrTh(:,u),pfdrRd(:,u)] = cart2pol(pfdr(:,u,1),pfdr(:,u,2));
    pfdha(:,u) = circ_dist(pfhrTh(:,u)-0.174532925199433,pfdrTh(:,u));    
end

ahd = acos(dot(pfhr(:,u,:),pfdr(:,u,:),3)./(sqrt(sum(pfhr(:,u,:).^2,3)).*sqrt(sum(pfdr(:,u,:).^2,3))));



hst = xyz(:,'hcom',:)-xyz(:,'spine_upper',:);
[th,phi] = cart2sph(hst(:,1,1),hst(:,1,2),hst(:,1,3));
thd = circ_dist(circshift(th,-5),circshift(th,5));
hrz = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);





spk = copy(Trial.spk);
spk.create(Trial,sampleRate,'theta-groom-sit-rear',unitSubset,'deburst');

tper = stcm(:,1)==1 & stcm(:,2)~=2;

figure();
for u = 1:numel(unitSubset)
%u = 20;
unit = unitSubset(u);
res = spk(unit);
phzspk = phz(res,spk.map(unit==spk.map(:,1),2));            
res(phzspk>-pi/4&(stcm(res,1)~=1|stcm(res,2)==2)) = [];
shft = 0;
subplot(521);plot(pft,unit,1,'text',[],true);title(num2str(unit));
subplot(522);plot(hsway(res),thd(res),'.');xlim([-0.1,0.1]);ylim([-0.1,0.1]);
subplot(523);plot(pfdha(tper,u),pfhr(tper,u,2),'.');grid('on');xlim([-pi,pi]);ylim([-300,300]);
subplot(524);hist2([pfdha(tper,u),pfhr(tper,u,2)],linspace(-pi,pi,21),linspace(-200,200,21));grid('on');
subplot(525);plot(pfdha(res,u),pfhr(res+shft,u,2),'.');grid('on');xlim([-pi,pi]);ylim([-300,300]);
subplot(526);hist2([pfdha(res,u),pfhr(res+shft,u,2)],linspace(-pi,pi,21),linspace(-200,200,21));grid('on');
% $$$ subplot(523);plot(hsway(res),pfhr(res+shft,u,2),'.');grid('on');xlim([-0.1,0.1]);ylim([-300,300]);
% $$$ subplot(524);hist2([hsway(res),pfhr(res+shft,u,2)],linspace(-0.1,0.1,21),linspace(-200,200,21));grid('on');
% $$$ subplot(525);plot(thd(res),pfhr(res+shft,u,2),'.');grid('on');xlim([-0.1,0.1]);ylim([-300,300]);
% $$$ subplot(526);hist2([thd(res),pfhr(res+shft,u,2)],linspace(-0.1,0.1,21),linspace(-200,200,21));grid('on');
subplot(527);plot(roll(tper),pfhr(tper,u,2),'.');grid('on');xlim([-pi/2,pi/2]);ylim([-300,300]);
subplot(528);hist2([roll(tper),pfhr(tper,u,2)],linspace(-pi/2,pi/2,41),linspace(-200,200,21));grid('on');
subplot(529);plot(roll(res),pfhr(res+shft,u,2),'.');grid('on');xlim([-pi/2,pi/2]);ylim([-300,300]);
subplot(5,2,10);hist2([roll(res),pfhr(res+shft,u,2)],linspace(-pi/2,pi/2,41),linspace(-200,200,21));grid('on');
waitforbuttonpress();
end


phzBins = linspace(-pi,pi,7);
%pfdhaBins = linspace(-pi,pi,10);
pfdhaBins = linspace(-pi,pi,16);

figure();
sp = tight_subplot(numel(phzBins)-1,numel(pfdhaBins),0.01,0.01);
sp = reshape(sp,[numel(pfdhaBins),numel(phzBins)-1])';

for u = 1:numel(unitSubset)

    %u = 20;
    unit = unitSubset(u);
    res = spk(unit);
    %res = res(stcm(res,1)==1&(stcm(res,3)==3|stcm(res,5)==5|stcm(res,6)==6));
    %res = res(stcm(res,1)==1&(stcm(res,3)==3));    
    %res = res(stcm(res,1)==1&(stcm(res,5)==5));        
    res = res(stcm(res,1)==1&(stcm(res,6)==6));            
    %res = res(stcm(res,1)==1);    
    
    shft = 0;
    phzInds   = discretize(phz(res,spk.map(unit==spk.map(:,1),2)),phzBins);
    pfdhaInds = discretize(pfdha(res,u),pfdhaBins);
    axes(sp(1,1));
    plot(pft,unit,1,'text',[],true);title(num2str(unit));
    for s = 1:numel(pfdhaBins)-1,
        for p = 1:numel(phzBins)-1,
            %cla();
            hold('on');            
            axes(sp(p,s+1));
            ind = pfdhaInds==s&phzInds==p;
            if sum(ind),
                plot(pfhr(res(ind),u,1),pfhr(res(ind),u,2),'.b');
                grid('on');
                xlim([-300,300]);
                ylim([-300,300]);
            end
        end
    end
    %waitforbuttonpress();
end

%subplot(526);hist2([pfdha(res,u),pfhr(res+shft,u,2)],linspace(-pi,pi,21),linspace(-200,200,21));grid('on');





figure,plot(dhsway(res),pfhr(res,u,2),'.')
figure,hist2([hsway(res),pfhr(res,u,2)],linspace(-0.1,0.1,21),linspace(-200,200,21));


figure,plot3(fhsway(res),pfhr(res,u,2),thd(res),'.')
figure,plot3(fhsway(res),pfhr(res,u,2),pfhr(res,u,1),'.')



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
    