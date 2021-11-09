% req20210705
% 
% XYHB 30sps 300ms

global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

sampleRate = 30;

%%%<<< XYHB Decoding -------------------------------------------------------------------------------

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

tind = [3:5,17:25];
cind = 1;
dc = {};
for t = tind;
    Trial = Trials{t};
    tunits = units{t};
    thetaRefChan = sessionList(t).thetaRefGeneral;

    halfSpkWindow = 0.15;
    ufrWindow = 0.15;

    smoothingWeights = [800.^2,800.^2, 1.2.^2, 1.2.^2];
    %smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];

    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    pfs = compute_xyhb_ratemaps( Trial, tunits);

    ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
    mask = ds.mask;

    spk = Trial.load('spk', sampleRate, '', tunits, 'deburst');
    ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

    dc{cind} = decode_ufr_boxcar(Trial, ...
                                 tunits, ...
                                 sampleRate, ...
                                 ufr, ...
                                 pfs, ...
                                 mask,...
                                 halfSpkWindow,...
                                 smoothingWeights,...
                                 'overwrite',true);

    dc{cind}.pfs = pfs;

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

    % COM 
    dc{cind}.ecom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.com},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.ecom = dc{cind}.ecom{1};
    % SAX 
    dc{cind}.esax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.sax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.esax = dc{cind}.esax{1};
    % MAX 
    dc{cind}.emax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.max},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.emax = dc{cind}.emax{1};
    % COM 
    dc{cind}.elom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lom},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elom = dc{cind}.elom{1};
    % SAX 
    dc{cind}.elax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elax = dc{cind}.elax{1};

    % HVANG
    dc{cind}.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
    xycoor = cat(2,...
                 dc{cind}.hvang(:,'spine_upper',[1,2])-dc{cind}.hvang(:,'bcom',[1,2]),...
                 dc{cind}.hvang(:,'nose',[1,2])-dc{cind}.hvang(:,'hcom',[1,2]));
    dc{cind}.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    % Positive: CCW (Left)     Negative: CW (Right)
    dc{cind}.hvang.data = circ_dist(circshift(dc{cind}.hvang.data(:,2),-10),...
                                    circshift(dc{cind}.hvang.data(:,2),+10));
    dc{cind}.hvang = dc{cind}.hvang.data(dc{cind}.ind);

    
    
% COMPUTE corrected hbang
    dc{cind}.hbang = filter(copy(xyz),'ButFilter',3,9,'low');    
    xycoor = cat(2,...
                 dc{cind}.hbang(:,'spine_upper',[1,2])-dc{cind}.hbang(:,'bcom',[1,2]),...
                 dc{cind}.hbang(:,'nose',[1,2])-dc{cind}.hbang(:,'hcom',[1,2]));
    dc{cind}.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hbang.data = circ_dist(dc{cind}.hbang.data(:,2),dc{cind}.hbang.data(:,1));
    dc{cind}.hbang = dc{cind}.hbang.data(dc{cind}.ind)+0.2+-0.4*double(cind>4);

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

ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,2) == 2 |dc{cind}.stcm(:,3) == 3 | ...
                                 dc{cind}.stcm(:,4) == 4 |dc{cind}.stcm(:,5) == 5 | ...
                                 dc{cind}.stcm(:,6) == 6);
figure();
hist2([dc{cind}.elax(ind,3),dc{cind}.esax(ind,3)],linspace(-1.25,1.25,30),linspace(-1.25,1.25,30),[],'mud')

figure();
hist2([dc{cind}.elax(ind,4),dc{cind}.esax(ind,4)],linspace(-1.25,1.25,30),linspace(-1.25,1.25,30),[],'mud')

ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,2) == 2);
ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,3) == 3);
ind = dc{cind}.stcm(:,1) == 1 & (dc{cind}.stcm(:,5) == 5);

% $$$  |dc{cind}.stcm(:,3) == 3 | ...
% $$$                                  dc{cind}.stcm(:,4) == 4 |dc{cind}.stcm(:,5) == 5 | ...
% $$$                                  dc{cind}.stcm(:,6) == 6
figure();
subplot(221)
histogram([dc{cind}.esax(ind,3)],linspace(-1.25,1.25,30));
ylim([0,2000]);
subplot(222)
histogram([dc{cind}.elax(ind,3)],linspace(-1.25,1.25,30));
ylim([0,2000]);
subplot(223)
histogram([dc{cind}.esax(ind,4)],linspace(-1.25,1.25,30));
ylim([0,5000]);
subplot(224)
histogram([dc{cind}.elax(ind,4)],linspace(-1.25,1.25,30));
ylim([0,5000]);

std(dc{cind}.esax(ind,4))
std(dc{cind}.elax(ind,4))
mean(dc{cind}.esax(ind,4))
mean(dc{cind}.elax(ind,4))

hberrSax = sqrt(sum(dc{cind}.esax(:,3:4).^2,2));
hberrLax = sqrt(sum(dc{cind}.elax(:,3:4).^2,2));

[mean(hberrSax(ind)),std(hberrSax(ind))]
[mean(hberrLax(ind)),std(hberrLax(ind))]


%%%>>>

