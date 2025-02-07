
% CA1
tind = [3,4,5,17,18,19,20,21,22,23,29];
% CA3
%tind = [6,7,26,27,30];
sampleRate = 250;
halfSpikeWindow = 0.020;
global AP
% compute_ratemaps ---------------------------------------------------------------------------------
AP.compute_ratemaps =                                                                            ...
    struct('get_featureSet',            @fet_xy,                                                 ...
           'sampleRate',                16,                                                      ...
           'pfsArgs',                   struct('states',           'theta-groom-sit-rear',       ...
                                               'binDims',          [50,50],                      ...
                                               'SmoothingWeights', [2.4,2.4],                    ...
                                               'numIter',          1,                            ...
                                               'boundaryLimits',   [-500,500;-500,500],          ...
                                               'halfsample',       false)                        ...
           );
%---------------------------------------------------------------------------------------------------    % LOAD decoded positions
dca = cf(@(T,U) accumulate_decoding_vars( T, U, sampleRate, ...
                                           halfSpikeWindow), Trials(tind),units(tind));

medD = [];
skwD = [];
stdD = [];
medR = [];
skwR = [];
stdR = [];
medL = [];
skwL = [];
stdL = [];
medC = [];
skwC = [];
stdC = [];

rdists = 5:5:45;
for r = 1:numel(rdists)
    clear('xcomp','ycomp','zcomp','ccomp');
    xcomp.label = 'hba';            xcomp.data = [];    xcomp.edgs = [-1.25,-0.2,0.2,1.25];
    ycomp.label = 'theta phase';    ycomp.data = [];    ycomp.edgs = linspace( 0.5, 2*pi-0.5,6 );
    ccomp.label = 'ego lat';        ccomp.data = [];    ccomp.clim = [-80,80];
    fcomp.data = [];
    mcomp.data = [];
    for t = [1:3,5:8,11],%{ER06,jg05,FS03}
        dc = dca{t};
        headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
        headAngle = atan2(headAngle(:,2),headAngle(:,1));
        mazeAngle =  sq(dca{t}.xyz(:,'hcom',[1,2]));
        mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
        headMazeAngle = circ_dist(headAngle, mazeAngle);
        mind =  dc.stcm(:,1)==1                                             ... State - Theta
                & (dc.stcm(:,3)==3|dc.stcm(:,4)==4|dc.stcm(:,5)==5)         ... State - Walk,Turn,Pause
                & dc.hvfl(:,1)>-2                                           ... FwdHeadSpeed
                & dc.ucnt>=4 & dc.ucnt<=10                                  ... UnitCount
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)<rdists(r)+5  ... DistanceFromMazeCenter
                & (sqrt(sum(dc.xyz(:,'hcom',[1,2]).^2,3))./10)>rdists(r)-5;   % DistanceFromMazeCenter
        %mind(mind==true) = randn([sum(mind),1])>0; % Subsample
        xcomp.data = cat(1, xcomp.data, dc.hbang(mind,1));
        ycomp.data = cat(1, ycomp.data, dc.phz(mind));
        ccomp.data = cat(1, ccomp.data, dc.esax(mind,2)./10);
        mcomp.data = cat(1, mcomp.data, headMazeAngle(mind));
    end
    [xcomp,ycomp,zcomp] = compute_2d_discrete_stats(xcomp,ycomp,ccomp);
    zmean   = zcomp.mean;    zmean(zcomp.count<5)  = nan;
    zmedian = zcomp.median;  zmean(zcomp.count<5)  = nan;
    zstd    = zcomp.std;     zstd(zcomp.count<5)   = nan;
    zcount  = zcomp.count;   zcount(zcomp.count<5) = nan;
    sectors = linspace(-pi,pi,11);
    sectorsL = linspace(-pi,pi,11);
    for s = 1:numel(sectors)-1
        if r>1
            bind = mcomp.data<sectors(s+1) & mcomp.data>sectors(s) ...
                   & ycomp.data>4.5& ycomp.data<5;
        else
            bind = ycomp.data>4.5& ycomp.data<5;
        end
        ind  = bind & xcomp.data<-0.2;        
          medR(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
          skwR(r,s) = skew(ccomp.data(ind,1));
          stdR(r,s) = std(ccomp.data(ind,1));
        ind = bind & xcomp.data<0.2& xcomp.data>-0.2;
          medC(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
          skwC(r,s) = skew(ccomp.data(ind,1));
          stdC(r,s) = std(ccomp.data(ind,1));
        ind = bind & xcomp.data>0.2;
          medL(r,s) = median(ccomp.data(ind,1));%-rmodel(r,s);
          skwL(r,s) = skew(ccomp.data(ind,1));
          stdL(r,s) = std(ccomp.data(ind,1));
        medD(r,s) = medR(r,s)-medL(r,s);
        skwD(r,s) = skwR(r,s)-skwL(r,s);
        stdD(r,s) = stdR(r,s)-stdL(r,s);
    end
end

sectorc = mean([sectors(2:end);sectors(1:end-1)]);
rdistc = mean([rdists(2:end);rdists(1:end-1)]);
rdiste = [0,rdists+2.5];
%[THETA,RR] = meshgrid(sectors,[rdists,rdists(end)]);
[THETA,RR] = meshgrid(sectors,[rdiste]);
[THETAC,RRC] = meshgrid(sectorc,rdists);
%[A,B] = pol2cart(circ_dist(THETA,diff(sectors([1,2]))),RR);
[X,Y] = pol2cart(THETA,RR);
[Xc,Yc] = pol2cart(THETAC,RRC);
mfun = @(beta,x) beta(2).*(sum(x.^2,2)).*sin(atan2(x(:,2),x(:,1)))+beta(1);
[beta,R,J,COVB,MSE,EMI] = nlinfit([Xc(:),Yc(:)].*10,medC(:).*10,mfun,[0,0.0001]);
