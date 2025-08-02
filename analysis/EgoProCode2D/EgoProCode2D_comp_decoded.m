function d = EgoProCode2D_comp_decoded(dca,state,trlgrp,mfun,beta)

switch trlgrp
  case 'CA1'
    tInds = [1:3,5:8,11];
  case 'ER06'
    tInds = [1:3];
  case 'jg05'
    tInds= [5:8];
  case 'FS03'
    tInds = [11];
  otherwise
    tInds = trlgrp;
end

% ACCUMULATE vars for specific subset of decoding
d = struct('fwd',[],...
                 'lat',[],...
                 'tfwd',[],...
                 'tlat',[],...
                 'xyz',[],...
                 'dst',[],...
                 'hvf',[],...
                 'hvl',[],...
                 'bvl',[],...
                 'bvf',[],...
                 'hvlF',[],...
                 'hvlP',[],...
                 'hav',[],...
                 'hbv',[],...           
                 'hba',[],...
                 'phz',[],...
                 'trl',[]);

for t = tInds,
    switch state
      case 'theta'
        stateInd = (dca{t}.stcm(:,3)==3 | dca{t}.stcm(:,4)==4 | dca{t}.stcm(:,5)==5);
      case 'loc'
        stateInd = (dca{t}.stcm(:,6)==6 | dca{t}.stcm(:,8)==8);
      case 'pause'
        stateInd = (dca{t}.stcm(:,7)==7 | dca{t}.stcm(:,9)==9);
    end
    mind =   dca{t}.stcm(:,1)==1        ... Theta state
           & stateInd                   ... states
           & dca{t}.ucnt>=3             ... coactive units
           & dca{t}.ucnt<=10            ...
           & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))>0 ...
           & sqrt(sum(dca{t}.xyz(:,'hcom',[1,2]).^2,3))<300;
    
    d.fwd  = cat(1, d.fwd, (dca{t}.ecom(mind,1)-22.5));
    d.lat  = cat(1, d.lat, (dca{t}.ecom(mind,2)-10*double(t>4)+10*double(t<4)));
    d.tfwd = cat(1, d.tfwd, dca{t}.tcom(mind,1));
    d.tlat = cat(1, d.tlat,dca{t}.tcom(mind,2));%+20*double(t>4));
    d.xyz  = cat(1, d.xyz, dca{t}.xyz(mind,{'hcom','nose'},:));
    d.dst  = cat(1, d.dst, dca{t}.hdist(mind));
    d.hvf  = cat(1, d.hvf, dca{t}.hvfl(mind,1));
    d.hvl  = cat(1, d.hvl, dca{t}.hvfl(mind,2));
    d.bvl  = cat(1, d.bvl, dca{t}.bvfl(mind,2));
    d.bvf  = cat(1, d.bvf, dca{t}.bvfl(mind,1));
    d.hvlF = cat(1, d.hvlF,dca{t}.hvflF(mind,2));
    d.hvlP = cat(1, d.hvlP,dca{t}.hvflP(mind,2));
    d.hav  = cat(1, d.hav, dca{t}.hvang(mind,1));
    d.hbv  = cat(1, d.hbv, dca{t}.hbavl(mind,1));
    d.hba  = cat(1, d.hba, dca{t}.hbang(mind,1));
    d.phz  = cat(1, d.phz, dca{t}.phz(mind,1));
    d.trl  = cat(1, d.trl, t*ones([sum(mind),1]));
end
headAngle = sq(d.xyz(:,2,[1,2])-d.xyz(:,1,[1,2]));
headAngle = atan2(headAngle(:,2), headAngle(:,1));
mazeAngle =  sq(d.xyz(:,1,[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
d.hma = circ_dist(headAngle,mazeAngle);
[Xe,Ye] = pol2cart(d.hma,d.dst);
d.correction.lat = mfun(beta,[Xe(:),Ye(:)]);
d.clat = (d.lat-d.correction.lat) ./ 10;
d.clt  = (d.lat-d.correction.lat) ./ 10;
d.fwd  = d.fwd ./ 10;
[d.ang, d.rad] = cart2pol( d.fwd, d.clat);
d.nsamples = numel(d.fwd);