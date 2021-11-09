% req20210629
% 
% DECODE: XY 
% WINDOW: 300ms
% BINNING: phz
% 
% SECTIONS: 
%    XY Decoding
%    CONVERT decoding cell to array
%    PLOT Asc and Dsc Theta phase: JPDF of decoding VAR (ESAX) and behvaior VAR ( HVFL, HBANG ) 




global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

sampleRate = 250;

phzBins = linspace(0,2*pi,25);


%%%<<< XY Decoding ---------------------------------------------------------------------------------
tind = [3:5,17:25];
cind = 1;
dc = {};
for t = tind;
    Trial = Trials{t};
    tunits = units{t};
    thetaRefChan = sessionList(t).thetaRefGeneral;
    phzCorr = phzCorrection(t);
    halfSpkWindow = 0.15;
    ufrWindow = 0.02;
    %smoothingWeights = [250.^2,250.^2];
    smoothingWeights = [250.^2,250.^2];

    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    phz = load_theta_phase(Trial,xyz,thetaRefChan,phzCorr);

    pfs = pfs_2d_theta(Trial,tunits);
    mask = create_tensor_mask(pfs.adata.bins(1:2));

    spk = Trial.load('spk', sampleRate, 'theta-sit-groom', tunits, 'deburst');
    ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

    dc{cind} = decode_ufr_binnedPhz(Trial, ...
                                    tunits, ...
                                    sampleRate, ...
                                    ufr, ...
                                    phz, ...
                                    pfs, ...
                                    mask,...
                                    [],[],[],...
                                    smoothingWeights,...
                                    'overwrite',false);

    dc{cind}.pfs = pfs;

    dc{cind}.iphz = discretize(dc{cind}.phz,phzBins);

    dc{cind}.xyz  = copy(xyz);
    dc{cind}.xyz.data  = xyz.data(dc{cind}.ind,:,:);

    dc{cind}.hRot = 0.17;
    dc{cind}.vxy  = filter(copy(xyz),'ButFilter',4,1.5,'low');;
    dc{cind}.vxy  = dc{cind}.vxy.vel('hcom',[1,2]);
    dc{cind}.vxy.data  = dc{cind}.vxy(dc{cind}.ind,1);
    dc{cind}.lvxy = copy(dc{cind}.vxy);
    dc{cind}.lvxy.data(dc{cind}.lvxy.data<=0) = 0.001;
    dc{cind}.lvxy.data = log10(dc{cind}.lvxy.data);

    dc{cind}.hvec = dc{cind}.xyz(:,'nose',[1,2])-dc{cind}.xyz(:,'hcom',[1,2]);
    dc{cind}.hvec = sq(bsxfun(@rdivide,dc{cind}.hvec,sqrt(sum(dc{cind}.hvec.^2,3))));
    dc{cind}.hvec = cat(3,dc{cind}.hvec,sq(dc{cind}.hvec)*[0,-1;1,0]);
    dc{cind}.hvec = multiprod(dc{cind}.hvec,                           ...
                              [cos(dc{cind}.hRot),-sin(dc{cind}.hRot); ...
                        sin(dc{cind}.hRot), cos(dc{cind}.hRot)],...
                              [2,3],                                   ...
                              [1,2]);

    dc{cind}.tvec = circshift(xyz(:,'hcom',[1,2]),-round(250*0.1)) - ...
                    circshift(xyz(:,'hcom',[1,2]),round(250*0.1));
    dc{cind}.tvec = sq(bsxfun(@rdivide,dc{cind}.tvec,sqrt(sum(dc{cind}.tvec.^2,3))));
    dc{cind}.tvec = cat(3,dc{cind}.tvec,sq(dc{cind}.tvec)*[0,-1;1,0]);
    dc{cind}.tvec = dc{cind}.tvec(dc{cind}.ind,:,:);

    dc{cind}.fet  = fet_HB_pitchB(Trial,dc{cind}.sampleRate);
    dc{cind}.fet.data  = dc{cind}.fet.data(dc{cind}.ind,:);


    dc{cind}.hdist = filter(copy(xyz),'ButFilter',3,5,'low');
    xycoor =    xyz(:,'hcom',[1,2]);
    [~,dc{cind}.hdist.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hdist.data = dc{cind}.hdist.data(dc{cind}.ind);


    dc{cind}.ecom = cf(@(e,x,h,f)  ...
                       sq(multiprod(bsxfun(@minus,...
                                           e(:,[1,2]),...
                                           sq(x(:,'hcom',[1,2]))),...
                                    h,...
                                    2,...
                                    [2,3])...
                          ), ...
                        {dc{cind}.com},{dc{cind}.xyz},{dc{cind}.hvec});
    dc{cind}.ecom = dc{cind}.ecom{1};
    % SAX 
    dc{cind}.esax = cf(@(e,x,h,f)  ...
                       sq(multiprod(bsxfun(@minus,...
                                           e(:,[1,2]),...
                                           sq(x(:,'hcom',[1,2]))),...
                                    h,...
                                    2,...
                                    [2,3])...
                          ), ...
                       {dc{cind}.sax},{dc{cind}.xyz},{dc{cind}.hvec});
    dc{cind}.esax = dc{cind}.esax{1};
    % MAX 
    dc{cind}.emax = cf(@(e,x,h,f)  ...
                       sq(multiprod(bsxfun(@minus,...
                                           e(:,[1,2]),...
                                           sq(x(:,'hcom',[1,2]))),...
                                    h,...
                                    2,...
                                    [2,3])...
                          ), ...
                       {dc{cind}.max},{dc{cind}.xyz},{dc{cind}.hvec});
    dc{cind}.emax = dc{cind}.emax{1};

    % HVANG
    dc{cind}.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
    xycoor = cat(2,...
                 dc{cind}.hvang(:,'spine_upper',[1,2])-dc{cind}.hvang(:,'bcom',[1,2]),...
                 dc{cind}.hvang(:,'nose',[1,2])-dc{cind}.hvang(:,'hcom',[1,2]));
    dc{cind}.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    % Positive: CCW (Left)     Negative: CW (Right)
    dc{cind}.hvang.data = circ_dist(circshift(dc{cind}.hvang.data(:,2),-10),...
                                    circshift(dc{cind}.hvang.data(:,2),+10));
    dc{cind}.hvang.data = dc{cind}.hvang.data(dc{cind}.ind);

    

% COMPUTE corrected hbang
    dc{cind}.hbang = filter(copy(xyz),'ButFilter',3,30,'low');    
    xycoor = cat(2,...
                 dc{cind}.hbang(:,'spine_upper',[1,2])-dc{cind}.hbang(:,'bcom',[1,2]),...
                 dc{cind}.hbang(:,'nose',[1,2])-dc{cind}.hbang(:,'hcom',[1,2]));
    dc{cind}.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hbang.data = circ_dist(dc{cind}.hbang.data(:,2),dc{cind}.hbang.data(:,1));
    dc{cind}.hbang.data = dc{cind}.hbang.data(dc{cind}.ind);
    dc{cind}.hbang.data = dc{cind}.hbang.data+0.2+-0.4*double(cind>4);
    
% GENERATE state collection matrix    
    dc{cind}.stcm = stc2mat(Trial.stc,xyz,states);
    dc{cind}.stcm = dc{cind}.stcm(dc{cind}.ind,:);

% COMPUTE head and body relative velocity vectors
    hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    dc{cind}.hvfl = hvfl(dc{cind}.ind,:);
    dc{cind}.bvfl = bvfl(dc{cind}.ind,:);    
    
    cind  = cind +1;
end                     

%%%>>>



%%%<<< CONVERT decoding cell to array --------------------------------------------------------------
cind = 1;
dca = dc{cind};
dca.sessionId = tind(cind).*ones(size(dc{cind}.ind));

%dca.xyz = dca.xyz.data;
dca.pos = sq(dca.xyz(:,'hcom',:));
dca.fet = dca.fet.data;
dca.vxy = dca.vxy.data;
dca.lvxy = dca.lvxy.data;
dca.hvang = dca.hvang.data;
dca.hbang = dca.hbang.data;
dca.hdist = dca.hdist.data;


for cind = 2:numel(dc),
    dca.sessionId = cat(1,dca.sessionId,tind(cind).*ones(size(dc{cind}.ind)));
    dca.ind = cat(1,dca.ind,dc{cind}.ind);
    dca.max = cat(1,dca.max,dc{cind}.max);
    dca.com = cat(1,dca.com,dc{cind}.com);    
    dca.sax = cat(1,dca.sax,dc{cind}.sax);        
    dca.post = cat(1,dca.post,dc{cind}.post);
    dca.phz = cat(1,dca.phz,dc{cind}.phz);    
    dca.iphz = cat(1,dca.iphz,dc{cind}.iphz);        
    dca.ucnt = cat(1,dca.ucnt,dc{cind}.ucnt);    
    %dca.uinc = cat(1,dca.uinc,dc{cind}.uinc);    
    %dca.xyz = cat(1,dca.xyz,dc{cind}.xyz.data);    
    dca.pos = cat(1,dca.pos,sq(dc{cind}.xyz(:,'hcom',:)));
    dca.vxy = cat(1,dca.vxy,dc{cind}.vxy.data);        
    dca.lvxy = cat(1,dca.lvxy,dc{cind}.lvxy.data);            
    dca.hvec = cat(1,dca.hvec,dc{cind}.hvec);    
    dca.tvec = cat(1,dca.tvec,dc{cind}.tvec);        
    dca.fet = cat(1,dca.fet,dc{cind}.fet.data);
    dca.hvang = cat(1,dca.hvang,dc{cind}.hvang.data);    
    dca.hbang = cat(1,dca.hbang,dc{cind}.hbang.data);        
    dca.hdist = cat(1,dca.hdist,dc{cind}.hdist.data);        
    
    dca.ecom  = cat(1,dca.ecom,dc{cind}.ecom);        
    dca.esax  = cat(1,dca.esax,dc{cind}.esax);        
    dca.emax  = cat(1,dca.emax,dc{cind}.emax);
    
    dca.stcm = cat(1,dca.stcm,dc{cind}.stcm);    
    dca.hvfl = cat(1,dca.hvfl,dc{cind}.hvfl);        
    dca.bvfl = cat(1,dca.bvfl,dc{cind}.bvfl);            
    
    dca.hRot  = cat(1,dca.hRot,dc{cind}.hRot);
    dca.pfs   = cat(1,dca.pfs,dc{cind}.pfs);
end

%%%>>>



%%%<<< PLOT Asc and Dsc Theta phase: JPDF of decoding VAR (ESAX) and behvaior VAR ( HVFL, HBANG ) --

figure
norm = 'xprob';
%norm = '';
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3 & dca.hvang<0.001;
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3 & dca.hdist<300 & ~ismember(dca.sessionId,[3,4,5]);;
sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<350 &  ~ismember(dca.sessionId,[3,4,5]) & dca.ucnt>=2;
%sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<350 &  dca.ucnt>=2;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & dca.sessionId==tind(cind);
subplot(221);
    ind = ismember(dca.iphz,[3:10]) &  sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('forward egocentric decoding');%    xlabel('ego forward');    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
    title('Decending Theta Phase');    
subplot(222);
    ind = ismember(dca.iphz,[16:23]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('forward egocentric decoding');%    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
    title('Ascending Theta Phase');            
subplot(223);
    ind = ismember(dca.iphz,[3:10]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('lateral egocentric decoding');%    xlabel('ego lateral');
    ylabel('head body Angle');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(224);
    ind = ismember(dca.iphz,[16:23]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('lateral egocentric decoding');%    xlabel('ego lateral');
    ylabel('head body Angle');
    Lines([],0,'k');
    Lines(0,[],'k');

%%%>>>    
    
    
%%%<<< PLOT Asc and Dsc Theta phase: JPDF of decoding VAR (ESAX) and behvaior VAR ( HVFL, HBANG ) --

figure();
norm = 'xprob';
%norm = '';
sind =    dca.stcm(:,1)==1 ...
       & (dca.stcm(:,5)==5 | dca.stcm(:,4)==4 | dca.stcm(:,3)==3 | dca.stcm(:,6)==6) ...
       &  dca.ucnt>=2 ...
       & dca.hbang<0;

subplot(221);
    ind = ismember(dca.iphz,[1:8])  & sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,50),...
           linspace(-10,80,20),norm);
    xlabel('forward egocentric decoding');
    ylabel('forward head speed');
    Lines([],0,'k');
    Lines(0,[],'k');
    title('Decending Theta Phase');
subplot(222);
    ind = ismember(dca.iphz,[10:18]) & sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,50),...
           linspace(-10,80,20),norm);
    xlabel('forward egocentric decoding');%    xlabel('ego forward');
    ylabel('forward head speed');
    Lines([],0,'k');
    Lines(0,[],'k');
    title('Ascending Theta Phase');        
subplot(223);
    ind = ismember(dca.iphz,[1:8]) & sind;
    hist2([dca.esax(ind,2),...
           dca.hvfl(ind,2)],...
           linspace(-200,200,50),...
           linspace(-50,50,20),norm);          
    %linspace(-1.5,1.5,20),norm);
    xlabel('lateral egocentric decoding');%    xlabel('ego lateral');
    ylabel('lateral Head speed');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(224);
    ind = ismember(dca.iphz,[10:18])  & sind;
    hist2([dca.esax(ind,2),...
           dca.hvfl(ind,2)],...
           linspace(-200,200,50),...
           linspace(-50,50,20),norm);
           %linspace(-1.5,1.5,20),norm);    
    xlabel('lateral egocentric decoding');%    xlabel('ego lateral');    xlabel('ego lateral');
    ylabel('lateral Head speed');
    Lines([],0,'k');
    Lines(0,[],'k');
ForAllSubplots('caxis([0,0.08]);');

%%%>>>
    
    

figure
norm = 'xprob';
%norm = '';
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3 & dca.hvang<0.001;
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3 & dca.hdist<300 & ~ismember(dca.sessionId,[3,4,5]);;
%sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<300 &  ~ismember(dca.sessionId,[3,4,5]) & dca.ucnt>=2;
sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<350 &  ~ismember(dca.sessionId,[3,4,5]) & dca.ucnt>=2& dca.hbang>0;;
%sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<350 &  dca.ucnt>=2;
%sind = (dca.stcm(:,5)==5 | dca.stcm(:,3)==3) & dca.hdist<350 &  dca.ucnt>=2 & dca.hbang>0;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & dca.sessionId==tind(cind);
subplot(221);
    ind = ismember(dca.iphz,[3:10]) &  sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(222);
    ind = ismember(dca.iphz,[16:23]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,1),...
           dca.hvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(223);
    ind = ismember(dca.iphz,[3:10]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,2),...
           dca.hvfl(ind,2)],...
           linspace(-200,200,30),...
           linspace(-60,60,20),norm);
    xlabel('ego lateral');
    ylabel('hvfl lateral');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(224);
    ind = ismember(dca.iphz,[16:23]) & dca.ucnt>=2  & sind;
    hist2([dca.esax(ind,2),...
           dca.hvfl(ind,2)],...
           linspace(-200,200,30),...
           linspace(-60,60,20),norm);
    xlabel('ego lateral');
    ylabel('hvfl lateral');
    Lines([],0,'k');
    Lines(0,[],'k');


phzBinc = (phzBins(2:end)+phzBins(1:end-1))./2;

msaxf = [];
for p = 1:24,
    ind = dca.iphz==p & dca.ucnt>=2  & dca.stcm(:,2)~=2;
    msaxf(p) = median(dca.esax(ind,1));
end

figure,
plot(phzBinc,msaxf);

figure
plot(msaxf);




figure,
hold('on');
subplot(121);
hist2([dc.ego.com(dc.iphz==10,2),...
      dc.ego.com(dc.iphz==10,1)],...
      linspace(-300,300,50),...
      linspace(-300,300,50));
subplot(122);
hist2([dc.ego.com(dc.iphz==20,2),...
      dc.ego.com(dc.iphz==20,1)],...
      linspace(-300,300,50),...
      linspace(-300,300,50));

figure,
norm = 'xprob';
ind = dc.ucnt>1 & dc.post>0.007;
%ind = ~ind;
subplot(121);
hist2([dc.ego.com(ind,1),...
      dc.phz(ind)],...
      linspace(-200,200,50),...
      phzBins,...
      norm);
Lines([],pi,'k');
xlabel('Forward');
ylabel('Theta Phase');
subplot(122);
hist2([dc.ego.com(ind,2),...
      dc.phz(ind)],...
      linspace(-200,200,50),...
      phzBins,...
      norm);
Lines([],pi,'k');
xlabel('Lateral');
ylabel('Theta Phase');

figure,
norm = 'xprob';
norm = '';
subplot(1,2,1);
ind = pi/8<dc.phz& dc.phz<pi & dc.ucnt>2 & dc.post>0.01;
hist2([dc.ego.sax(ind,1),...
      dc.lvel(ind,1)],...
      linspace(-200,200,50),...
      linspace(0,1.8,20),...
      norm);
xlabel('ego forward');
ylabel('log10(speed) cm/s');
subplot(1,2,2);
hist2([dc.ego.sax(ind,2),...
       dc.lvel(ind,1)],...
      linspace(-200,200,50),...
      linspace(0,1.8,20),...
      norm);
xlabel('ego lateral');
ylabel('log10(speed) cm/s');


figure,
ind = pi/8<dc.phz& dc.phz<pi & dc.ucnt>2 & dc.post>0.007;
hist2([dc.ego.sax(ind,1),...
      dc.vel(ind,1)],...
      linspace(-200,200,50),...
      linspace(0,60,20),'xprob');
xlabel('ego forward');
ylabel('speed cm/s');

figure,
ind = pi<dc.phz& dc.phz<1.8*pi & dc.ucnt>2 & dc.post>0.007; %& dc.hdist.data< 350;
hist2([dc.ego.sax(ind,2),...
      dc.hbang(ind,1)-0.2],...
      linspace(-200,200,50),...
      linspace(-1.5,1.5,20),'xprob');
xlabel('ego lateral');
ylabel('hba rad');
Lines([],0,'k');
Lines(0,[],'k');


