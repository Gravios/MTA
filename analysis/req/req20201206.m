% Current problem redefinition of place cell selection criteria
%
% Interneuron classification criteria
%  - theta phase-locking strength
%  - theta phase-preference
%
figure();
plot(xyzp(tper,1),xyzp(tper,2),'.');
Lines([],[-500:100:500],'k');
Lines([-500:100:500],[],'k');

%req20201117(Trial)

% Interneuron summary
figure()
for u = unitsInt,
    
end


MjgER2016_load_data();

phzCorrection = [pi*1.25,pi*1.25,                             ... er01
                 pi/2,pi/2,pi/2,                    ... ER06
                 pi/1.25,pi/1.25,                             ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04                 
                 pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4 ... jg05
                 pi/4,pi/4,pi,pi/1.25,pi*1.25]; % new units - jg05, jg05, ER06, Ed10, er01

cf(@(t,p) report_interneuron_summary(t,'phzCorrection',p), Trials(end), num2cell(phzCorrection(end)));




unitsPyr = {...
    [161],...           er01 20110719
    [],...              er01 20110721
    [],...              ER06 20130612
    [],...              ER06 20130613
    [115,118,152],...   ER06 20130614
    [],...              Ed10 20140816
    [],...              Ed10 20140817
    [],...              jg04 20120128
    [1],...             jg04 20120129
    [18],...            jg04 20120130
    [22],...            jg04 20120131
    [],...              jg04 20120201
    [],...              jg04 20120210
    [],...              jg04 20120211    
    [],...              jg04 20120212
    [19],...            jg04 20120213
    [],...              jg05 20120309
    [90],...            jg05 20120310
    [],...              jg05 20120311
    [171],...           jg05 20120312
    [],...              jg05 20120315
    [],...              jg05 20120316
    [65],...            jg05 20120317
    [49,65],...         jg05 20120323    
    [3],...             jg05 20120324
    [],...              ER06 20130624
    [],...              Ed10 20140815
    [115,116,184] ...   er01 20110722
};

unitsInts = {...
    [ 31, 78, 82,125,169,195,203],...                                    er01 20110719
    [ 31, 76,105,147],...                                                er01 20110721
    [ 27, 32, 68, 69,105,124,125,126,128,182,220,221,222,225],...        ER06 20130612
    [ 10, 13, 15, 18, 60, 93,112,113,192,213,215,216,220],...            ER06 20130613
    [  5, 10, 11, 12, 22, 43, 49, 91, 94, 95,110,124,144,145,146,155,... ER06 20130614
       156,168,169,170,171,172,176,177,179,181,182],...                  
    [  9, 11, 12, 29, 30, 43, 45, 46, 51, 52, 54, 55, 56, 87],...        Ed10 20140816
    [ 11, 12, 31, 47, 48, 51, 56, 80, 81, 89],...                        Ed10 20140817
    [  8,  9, 16],...                                                    jg04 20120128
    [ 21],...                                                            jg04 20120129
    [ 24],...                                                            jg04 20120130
    [ 10, 24, 27],...                                                    jg04 20120131
    [ 10],...                                                            jg04 20120201
    [  4,  5],...                                                        jg04 20120210
    [  4],...                                                            jg04 20120211
    [  6],...                                                            jg04 20120212
    [  2],...                                                            jg04 20120213
    [  5, 10, 15, 27, 28, 38, 64, 66,100,114,116,117,121,122],...        jg05 20120309
    [  4,  6,  7,  8, 27, 28, 43, 59, 71, 86, 99,100,101,102],...        jg05 20120310    
    [  2,  7,  9, 25, 41, 44, 46, 70, 71, 72, 98,112,113,118,...         jg05 20120311
     119,135,144,148,188,193,203,204,205],...
    [  3,  7,  8, 15, 16, 43, 45, 50, 76, 77, 92,106,124,184],...        jg05 20120312
    [  1,  5, 34, 60, 69],...                                            jg05 20120315
    [ 17, 28, 49, 55],...                                                jg05 20120316
    [ 15, 16, 17, 52],...                                                jg05 20120317
    [  2,  4, 20, 30, 39, 40, 51],...                                    jg05 20120323
    [  7, 14, 16, 31, 34],...                                            jg05 20120324
    [  5, 13, 27, 66, 71, 72, 75, 91, 94,105,127,157,209,232,233,251,254],...ER06 20130624
    [  4, 17, 20, 22, 23, 25, 29 ,30, 31, 40, 66, 69, 72],...            Ed10 20140815
    [ 93,175] ...                                                        er01 20110722
};


tppRUN = [];
tprRUN = [];
tppREM = [];
tprREM = [];
for t = 1:numel(Trials)
    Trial = Trials{t};    
    lfp = Trials{t}.load('lfp',sessionList(t).thetaRefGeneral);
    phz = lfp.phase([5,13]);    
    phz.data = unwrap(phz.data);
    phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection(t); 
    phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;
    
    spk = Trial.spk.copy();
    spk.create(Trial,              ...  
               phz.sampleRate,     ...
               '',                 ...
               unitsInts{t},       ...
               ''                  ...
    );    
    stc = Trials{t}.stc.copy();    
    runPer = [stc{'theta-groom-sit',phz.sampleRate}];
    remPer = [stc{'theta&sit'}];              
    for u = 1:numel(unitsInts{t}),
        res = spk(unitsInts{t}(u));               
        resRUN = res(WithinRanges(res,runPer.data));        
        resREM = res(WithinRanges(res,remPer.data));
        if numel(resRUN)>3,
            tppRUN(end+1) = circ_mean(phz(res(WithinRanges(res,runPer.data))));
            tprRUN(end+1) = circ_r(phz(res(WithinRanges(res,runPer.data))));
        else,
            tppRUN(end+1) = nan;
            tprRUN(end+1) = nan;            
        if numel(resREM)> 3,
            tppREM(end+1) = circ_mean(phz(res()));
            tprREM(end+1) = circ_r(phz(res(WithinRanges(res,remPer.data))));        
        else,
            tppREM(end+1) = nan;            
            tprREM(end+1) = nan;
        end
    end%for u
end%for t









t = 20;    
Trial = Trials{t};
spk = Trial.spk.copy();
spk.create(Trial,              ...  
           sampleRate,         ...
           '',                 ...
           unitsInts{t},       ...
           ''                  ...
);    
pft  = pfs_2d_theta(Trial);
%pfis = pfs_2d_states(Trial.

stc = Trial.stc.copy();    
rper = [stc{'R&s',sampleRate}];
sper = [stc{'s',sampleRate}];
rper = [stc{'r',sampleRate}];
rper.data(diff(rper.data,1,2)<200,:) = [];
%sper = [Trial.stc{'t',sampleRate}];
rper = stc.get_state_transitions(Trial,{'rear','pause'}, 0.2, xyz);    
rper = rper-80;    
ron  = stc.get_state_transitions(Trial,{'walk','rear' }, 0.2, xyz);  ron  = mean(ron, 2);
roff = stc.get_state_transitions(Trial,{'rear','pause'}, 0.2, xyz);  roff = mean(roff,2);
ronp  = stc.get_state_transitions(Trial,{'pause','rear' }, 0.2, xyz);  ronp  = mean(ronp, 2);
rofft = stc.get_state_transitions(Trial,{'rear','turn'}, 0.2, xyz);  rofft = mean(rofft,2);

figure();
    hold('on');
    plot(xyz(:,'spine_upper',3));
    Lines(rofft,[],'r');
    

TrigRasters(rofft,1000,res,92*ones(size(res)),250);
%pch = fet_HB_pitchB(Trial,sampleRate);

figure();
plot(xyz(:,'spine_upper',3));
Lines(rper(:),[],'k');

frq = copy(phz);
frq.data = circ_dist(frq.data,circshift(frq.data,1))/(2*pi)*sampleRate;    
figure,plot(frq.data)


figure,
    t = 20;
for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,rper));
clf();
scatter(frq(res),phz(res),10,xyz(res,'spine_upper',3),'filled');
colormap('jet');    
xlim([0,15]);
waitforbuttonpress();
end

figure();
for u = unitsInts{20};
res = spk(u);
%res = res(WithinRanges(res,sper));
%[mccg,tbins] = CCG([mean(rper.data,2);res],[ones([size(rper,1),1]);2*ones(size(res))],1,100,sampleRate,[1,2]);
subplot2(2,3,1,1);
    plot(pft,u,1,'text');title(num2str(u));

subplot2(2,3,1,2);
    [mccg,tbins] = CCG([ron(:,1);res],[ones([size(ron,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
    bar(tbins,mccg(:,1,2));
    title('walk to rear');
    
subplot2(2,3,1,3);
    [mccg,tbins] = CCG([ronp(:,1);res],[ones([size(ronp,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
    bar(tbins,mccg(:,1,2));
    title('pause to rear');

subplot2(2,3,2,2);
    [mccg,tbins] = CCG([roff(:,1);res],[ones([size(roff,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
    bar(tbins,mccg(:,1,2));
    title('rear to pause');
    
subplot2(2,3,2,3);
    [mccg,tbins] = CCG([rofft(:,1);res],[ones([size(rofft,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
    bar(tbins,mccg(:,1,2));
    title('rear to turn');
waitforbuttonpress();
end


sper = [Trial.stc{'t-m-s',sampleRate}];
sper = [Trial.stc{'s-t',sampleRate}];
sper = [Trial.stc{'s-t',sampleRate}];

sper = [Trial.stc{'R-t',sampleRate}]+[-0.2,0.2];
u = [ 92,106];
figure();
for x = 1:14,
    for y = 1:14,
        u = [ unitsInts{t}(x),unitsInts{t}(y)];
        res1 = spk(u(1));res1 = res1(WithinRanges(res1,sper));
        res2 = spk(u(2));res2 = res2(WithinRanges(res2,sper));
        [mccg,tbins] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))],1,70,sampleRate,[1,2]);
        subplot2(14,14,y,x);
        bar(tbins,mccg(:,1,2));
        if x == y,
            ylim([0,75]);
        end
    end
end



figure,
gfor u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,sper));
pshift = -8*pi:pi/4:8*pi;
%pshift = -250:25:250;
mccg = zeros([101,numel(pshift)]);
for p = 1:numel(pshift)
rper = LocalMinima(circ_dist(phz.data,pshift(p))+pi,10,0.1);
%rper = LocalMinima(circshift(phz.data,pshift(p))+pi,10,0.1);
[tccg,tbins] = CCG([rper(:,1);res],[ones([size(rper,1),1]);2*ones(size(res))],4,50,sampleRate,[1,2]);
mccg(:,p) = tccg(:,1,2);
end;
clf();
imagesc(mccg');
Lines(51,[],'k');
Lines([],numel(pshift)/2,'k');
colormap('jet');
axis('xy');
waitforbuttonpress();
end


figure,
for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,sper));
tic
rper = LocalMinima(circ_dist(phz.data,-circ_mean(phz(res))),10,0.1);
toc
[mccg,tbins] = CCG([rper(:,1);res],[ones([size(rper,1),1]);2*ones(size(res))],4,50,sampleRate,[1,2]);
clf();
subplot(121);
rose(phz(res),36);
subplot(122);hold('on');
bar(tbins,mccg(:,1,2)-mean(mccg(:,1,2)));
plot(tbins,normalize(RectFilter(conv((mccg(:,1,2)-mean(mccg(:,1,2))).^2,ones([11,1]),'same'),3,5),'range'),'r');
Lines(tbins(51),[],'k');
waitforbuttonpress();
end





% Speed X Interneurons

clear('vel');
vel = vel(filter(copy(xyz),               ... % copy <= position data
                 'ButFilter',             ... % filter <= filter class
                 4,                       ... % filter <= num pole 
                 1.5,                     ... % filter <= freq boundaries
                 'low'),                  ... % filter <= filter mode
          {'spine_lower','hcom'},         ... % vel <= markers
          [1,2]                           ... % vel <= xyz subspace
);


% $$$ velccg.hst.eds = [0:5:80];
% $$$ velccg.hst.ctr = mean([velccg.hst.eds(1:end-1);velccg.hst.eds(2:end)]);
% $$$ velccg.hst.ind = discretize(vel(:,2),velccg.hst.eds);

velccg.hst.eds = [-2:0.2:1.8];
velccg.hst.ctr = mean([velccg.hst.eds(1:end-1);velccg.hst.eds(2:end)]);
velccg.hst.ind = discretize(log10(vel(:,1)),velccg.hst.eds);

velccg.ccg.opts.halfBins = 100;
velccg.ccg.opts.binSize = 1;
velccg.ccg.opts.sampleRate = spk.sampleRate;
velccg.ccg.data = zeros([velccg.ccg.opts.halfBins*2+1,numel(velccg.hst.ctr)]);
velccg.ccg.time = zeros([velccg.ccg.opts.halfBins*2+1,1]);

%u = [92,106];
u = [15];

for v = 1:numel(velccg.hst.ctr),
    res = spk(u);
    res = res(velccg.hst.ind(res)==v);
    [velccg.ccg.data(:,v),velccg.ccg.time] = ... % CCG( T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
        CCG(res,                             ... % CCG <= T
            ones(size(res)),                 ... % CCG <= G
            velccg.ccg.opts.binSize,         ... % CCG <= BinSize
            velccg.ccg.opts.halfBins,        ... % CCG <= HalfBins
            velccg.ccg.opts.sampleRate,      ... % CCG <= SampleRate
            1                                ... % CCG <= GSubset
    ); % CCG
end

figure();
imagesc(velccg.ccg.time,velccg.hst.ctr,velccg.ccg.data');
axis('xy');
colormap('jet');




figure();
imagesc(velccg.ccg.time,velccg.hst.ctr,bsxfun(@rdivide,velccg.ccg.data,sum(velccg.ccg.data))');
axis('xy');
colormap('jet');
caxis([0,0.01]);




% hvf(l) X Interneurons
hfl = fet_href_HXY(Trial,sampleRate,false,'trb',2);

hflccg.hst.eds = [0:5:70];
%hflccg.hst.eds = [-22.5:5:70];
%hflccg.hst.eds = [-62.5:5:62.5];
hflccg.hst.ctr = mean([hflccg.hst.eds(1:end-1);hflccg.hst.eds(2:end)]);
%hflccg.hst.ind = discretize(hfl(:,1),hflccg.hst.eds);

hflccg.hst.ind = discretize(sqrt(sum(hfl(:,:).^2,2)),hflccg.hst.eds);

hflccg.ccg.opts.halfBins = 100;
hflccg.ccg.opts.binSize = 1;
hflccg.ccg.opts.sampleRate = spk.sampleRate;
hflccg.ccg.data = zeros([hflccg.ccg.opts.halfBins*2+1,numel(hflccg.hst.ctr)]);
hflccg.ccg.time = zeros([hflccg.ccg.opts.halfBins*2+1,1]);

%u = [92,106];
u = [92];

for vf = 1:numel(hflccg.hst.ctr),
    res = spk(u);
    res = res(hflccg.hst.ind(res)==v);
    [hflccg.ccg.data(:,v),hflccg.ccg.time] = ... % CCG( T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
        CCG(res,                             ... % CCG <= T
            ones(size(res)),                 ... % CCG <= G
            hflccg.ccg.opts.binSize,         ... % CCG <= BinSize
            hflccg.ccg.opts.halfBins,        ... % CCG <= HalfBins
            hflccg.ccg.opts.sampleRate,      ... % CCG <= SampleRate
            1                                ... % CCG <= GSubset
    ); % CCG
end

figure();
imagesc(hflccg.ccg.time,hflccg.hst.ctr,hflccg.ccg.data');
axis('xy');
colormap('jet');

figure();
imagesc(hflccg.ccg.time,hflccg.hst.ctr,bsxfun(@rdivide,hflccg.ccg.data,sum(hflccg.ccg.data))');
axis('xy');
colormap('jet');
caxis([0,0.01]);




% hvfl X Interneurons
hfl = fet_href_HXY(Trial,sampleRate,false,'trb',2);

rottheta = -0.17;
rotmat = [cos(rottheta),-sin(rottheta);sin(rottheta),cos(rottheta)];
ind = 10000:10200;
figure,
    axes();
    hold('on');
    plot(hfl(ind,1),hfl(ind,2));
    temphfl = multiprod(rotmat,hfl(ind,:),[1,2],[2]);
    plot(temphfl(:,1),temphfl(:,2));
    
figure
ind = stc{'t'};
subplot(121);
    rottheta = 0.0;
    rotmat = [cos(rottheta),-sin(rottheta);sin(rottheta),cos(rottheta)];
    hist2(multiprod(rotmat,hfl(ind,:),[1,2],[2]),linspace(-80,80,40),linspace(-80,80,40));
    caxis([0,600]);
    Lines([],0,'r');
subplot(122);
    rottheta = -0.17;
    rotmat = [cos(rottheta),-sin(rottheta);sin(rottheta),cos(rottheta)];
    hist2(multiprod(rotmat,hfl(ind,:),[1,2],[2]),linspace(-80,80,40),linspace(-80,80,40));
    caxis([0,600]);
    Lines([],0,'r');



hfl.data = multiprod(rotmat,hfl.data,[1,2],[2]);


hvfccg.hst.eds = [-30:20:90];
hvfccg.hst.ctr = mean([hvfccg.hst.eds(1:end-1);hvfccg.hst.eds(2:end)]);
hvfccg.hst.ind = discretize(hfl(:,1),hvfccg.hst.eds);

hvlccg.hst.eds = [-55:10:55];
hvlccg.hst.ctr = mean([hvlccg.hst.eds(1:end-1);hvlccg.hst.eds(2:end)]);
hvlccg.hst.ind = discretize(hfl(:,2),hvlccg.hst.eds);

hvflccg.ccg.opts.halfBins = 100;
hvflccg.ccg.opts.binSize = 1;
hvflccg.ccg.opts.sampleRate = spk.sampleRate;
hvflccg.ccg.data = zeros([hvflccg.ccg.opts.halfBins*2+1,numel(hvfccg.hst.ctr),numel(hvlccg.hst.ctr)]);
hvflccg.ccg.time = zeros([hvflccg.ccg.opts.halfBins*2+1,1]);

%u = [92,106];
u = [92];
u = [76];
u = [50];

for vf = 1:numel(hvfccg.hst.ctr),
    for vl = 1:numel(hvlccg.hst.ctr),
    res = spk(u);
    res = res(hvfccg.hst.ind(res)==vf & hvlccg.hst.ind(res)==vl);
    [hvflccg.ccg.data(:,vf,vl),hvflccg.ccg.time] = ... % CCG( T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
        CCG(res,                               ... % CCG <= T
            ones(size(res)),                   ... % CCG <= G
            hvflccg.ccg.opts.binSize,          ... % CCG <= BinSize
            hvflccg.ccg.opts.halfBins,         ... % CCG <= HalfBins
            hvflccg.ccg.opts.sampleRate,       ... % CCG <= SampleRate
            1                                  ... % CCG <= GSubset
    ); % CCG
    end
end


figure();
sp = flipud(tight_subplot(numel(hvfccg.hst.ctr),1,0.01));
for vf = 1:numel(hvfccg.hst.ctr),
    axes(sp(vf));
    imagesc(hvflccg.ccg.time,hvlccg.hst.ctr,sq(hvflccg.ccg.data(:,vf,:))');
     axis('xy');
    colormap('jet');
    caxis([0,50]);
end

figure();
sp = flipud(tight_subplot(numel(hvfccg.hst.ctr),1,0.01));
for vf = 1:numel(hvfccg.hst.ctr),
    axes(sp(vf));
    imagesc(hvflccg.ccg.time,hvlccg.hst.ctr,bsxfun(@rdivide,sq(hvflccg.ccg.data(:,vf,:)),sum(sq(hvflccg.ccg.data(:,vf,:))))');
     axis('xy');
    colormap('jet');
    caxis([0,0.01]);
end

figure();
imagesc(hflccg.ccg.time,hflccg.hst.ctr,bsxfun(@rdivide,hflccg.ccg.data,sum(hflccg.ccg.data))');
axis('xy');
colormap('jet');


% Speed X Interneurons
rbm = {'hcom','nose'};
hed = xyz.copy();
hed.data = hed(:,rbm,:);
hed.model = xyz.model.rb(rbm);
hed.filter('ButFilter',4,2,'low');
ang = create(MTADang,Trial,hed);

hav = ang.copy();
hav.data = circ_dist(circshift(ang(:,1,2,1),-1),ang(:,1,2,1)).*sampleRate;
hav.data = circ_dist(circshift(ang(:,1,2,2),-1),ang(:,1,2,2)).*sampleRate;

% if theta is driven by vestibular input 
% if theta is driven by visual path
% if theta is driven by trajectory



havccg.hst.eds = [-2:0.3:1];
havccg.hst.ctr = mean([havccg.hst.eds(1:end-1);havccg.hst.eds(2:end)]);
havccg.hst.ind = discretize(log10(abs(hav(:,1))),havccg.hst.eds);

havccg.ccg.opts.halfBins = 100;
havccg.ccg.opts.binSize = 1;
havccg.ccg.opts.sampleRate = spk.sampleRate;
havccg.ccg.data = zeros([havccg.ccg.opts.halfBins*2+1,numel(havccg.hst.ctr)]);
havccg.ccg.time = zeros([havccg.ccg.opts.halfBins*2+1,1]);

%u = [92,106];
u = [45];

for v = 1:numel(havccg.hst.ctr),
    res = spk(u);
    res = res(havccg.hst.ind(res)==v);
    [havccg.ccg.data(:,v),havccg.ccg.time] = ... % CCG( T, G, BinSize, HalfBins, SampleRate, GSubset, Normalization, Epochs)
        CCG(res,                             ... % CCG <= T
            ones(size(res)),                 ... % CCG <= G
            havccg.ccg.opts.binSize,         ... % CCG <= BinSize
            havccg.ccg.opts.halfBins,        ... % CCG <= HalfBins
            havccg.ccg.opts.sampleRate,      ... % CCG <= SampleRate
            1,                               ... % CCG <= GSubset
            'hz'                             ... % CCG <= Normalization
    ); % CCG
end

figure();
imagesc(havccg.ccg.time,havccg.hst.ctr,havccg.ccg.data');
axis('xy');
colormap('jet');
colorbar();

figure();
imagesc(havccg.ccg.time,havccg.hst.ctr,bsxfun(@rdivide,havccg.ccg.data,sum(havccg.ccg.data))');
axis('xy');
colormap('jet');
caxis([0,0.01]);





figure,
%for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,sper));




    
cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials(1),unitsInts(1));
tag = 'interneurons_xyhb_2020';
tag = 'interneurons_xyhb_2020_final';    overwrite = false;
tag = 'interneurons_xyhb_2020_loc';      overwrite = false;
tag = 'interneurons_xyhb_2020_loc_rear'; overwrite = true;
tag = 'interneurons_xyhb_2020_theta';    overwrite = true;
pfi = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);
[rmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(pfi,unitsInts);

pfti = cf(@(t,u)   pfs_2d_theta(t,u),                             Trials, unitsInts);
pfsi = cf(@(t,u)   pfs_2d_states(t,u,'pfsArgsOverride',struct('numIter',1,'halfsample',false)), Trials, unitsInts);
bfsi = cf(@(t,u)   compute_bhv_ratemaps(t,u,[],[],[],[],1,1000),  Trials, unitsInts);

pfsiStates = {'Loc','L Loc','H Loc','Rear','Pause','L Pause','H Pause'};

MjgER2016_load_data();
configure_default_args();

% LOAD place restricted behavior fields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);
%pfsa = cf(@(s,t) cat(2,{t},s), pfss,pfts);
% COMPUTE bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units);
numComp = size(eigVecs,2);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims));
    fpc{i}(validDims) = eigVecs(:,i);
end
fpcLims = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

figure,
imagesc(reshape_eigen_vector(fpc{4},bfs(1)))
axis('xy');



bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';
bhvLims = [-1.6, 0.6; ...
           -0.5, 1.7];


whos rmaps

bins = pfi{1}.adata.bins;
binDims = pfi{1}.adata.binSizes';

% RESHAPE rmaps into [ xPosition, yPosition, headPitch, bodyPitch, unit ]
rmapa = nan([binDims,size(rmaps,2)]);
for u = 1:size(rmaps,2),
    rmapa(:,:,:,:,u) = reshape(rmaps(:,u),binDims);
end

% PLOT example
rmap = rmapa(:,:,:,:,201);
figure();
for x = 1:6,
    for y = 1:6,
        subplot2(6,6,7-y,x);
        imagescnan({bins{3:4},sq(rmap(x,y,:,:))'});
        axis('xy');
    end
end

% bmaps = reshape(rmapa(3,4,:,:,:),[],size(rmapa,length(binDims)+1));

% COLLECT inner spatial bins
bmaps = [];
for k = 2:4,
    for j = 2:4,
        bmaps(:,:,k,j) = reshape(rmapa(k,j,:,:,:),[],size(rmapa,length(binDims)+1));
    end
end
bmaps = reshape(bmaps,[size(bmaps,1),prod(size(bmaps))/size(bmaps,1)]);

% PLOT valid elements of the behavior space
figure();
imagesc(reshape(validDims,binDims(end-1:end))');
axis('xy');


% COMPRESS bmaps to only valid behavior space elements
vmaps = bmaps(validDims,:);
vmaps(:,sum(isnan(vmaps))>0)= [];
vmaps(isnan(vmaps(:)))=0;

% COMPUTE factors
[LU,LR,FSr,VT] = erpPCA(vmaps',5);

% PLOT first 5 eigen vectors
figure()
for v = 1:5,
    eigVec = nan(binDims(end-1:end));
    eigVec(validDims) = LU(:,v);
    subplot(1,6,v);
    imagescnan({bins{end-1:end},eigVec'});
    axis('xy');
end
    subplot(1,6,6);
    plot(VT(1:5,4),'-+')


% END overall erpPCA 
% START individual unit non-negative matrix factorization (nnmf)

[accg,tbins] = cf(@(t,u)  autoccg(t,u), Trials,unitsInts);

tbins = tbins{1};
accg = cf(@(a,u) a(:,u),   accg, unitsInts);
accg = cat(2,accg{:});



imaps = bmaps(validDims,:);
imaps = reshape(imaps, [sum(validDims), 4*4, 225]);
%imaps(isnan(imaps)) = 0;
umaps = imaps(:,:,u);
%umaps(:,sum(isnan(umaps(:,:)))>0)= [];
invalidDims = sum(isnan(umaps),2)>0;
umaps(invalidDims,:) = [];

interMapCorrCoef = [];
for x = 1:16,
    for y = 1:16,

        tcc = corrcoef(umaps(:,x),umaps(:,y));
        interMapCorrCoef(x,y) = tcc(1,2);
    end
end
figure(),
imagesc(reshape(interMapCorrCoef(1,:),4,4)');,axis('xy');
figure,hist(nonzeros(triu(interMapCorrCoef,1)),40)

%% Main Plot -----------------------------------------------------

cf(@(t) t.load('nq'), Trials);
nq = cf(@(t) rmfield(t.nq,'AvSpk'), Trials);
nq = cf(@(n) StructArray(n), nq);
nq = cf(@(n,u) n(u), nq, unitsInts);
nq = cat(1,nq{:});

u = 158;
u = 127;
u = 166;
[sesId,cluId] = deal(cluSessionMap(u,1),cluSessionMap(u,2));


umaps = imaps(:,:,u);
%umaps(:,sum(isnan(umaps(:,:)))>0)= [];
invalidDims = sum(isnan(umaps),2)>0;
umaps(invalidDims,:) = [];

% LOAD Session data
% $$$ stc = Trials{sesId}.stc.copy();    
% $$$ lfp = Trials{sesId}.load('lfp',sessionList(sesId).thetaRefGeneral);
% $$$ phz = lfp.phase([5,13]);    
% $$$ phz.data = unwrap(phz.data);
% $$$ %phz.resample(xyz);    
% $$$ phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; 
% $$$ phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;
% $$$ clear('lfp');
% $$$ 
% $$$ xyz = preproc_xyz(Trial,'trb');
% $$$ xyz.data = xyz(:,{'head_back','head_left','head_front','head_right'},:);
% $$$ xyz.model = xyz.model.rb({'head_back','head_left','head_front','head_right'});
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ 
% $$$ spkL = Trial.load('spk',phz.sampleRate,'theta-sit-groom',unitsInts{sesId});
% $$$ spkX = Trial.load('spk',xyz.sampleRate,'loc&theta-sit-groom',unitsInts{sesId});
% $$$ 
% $$$ tper = resample(cast([stc{'t&x'}],'TimeSeries'),xyz);
% $$$ tper = logical(tper.data);
% $$$ 
% $$$ hdBinEds = linspace(-pi,pi,37);
% $$$ hdBinCtr = circ_mean([hdBinEds(2:end);hdBinEds(1:end-1)]);
% $$$ hdBinInd = discretize(ang(:,'head_back','head_front',1),hdBinEds);
% $$$ 

% $$$ spkR = Trial.load('spk',xyz.sampleRate,'rear&theta-sit-groom',unitsInts{sesId});
% $$$ sper = resample(cast([stc{'t&r'}],'TimeSeries'),xyz);
% $$$ sper = logical(sper.data);

% $$$ spkP = Trial.load('spk',xyz.sampleRate,'pause&theta-sit-groom',unitsInts{sesId});
% $$$ pper = resample(cast([stc{'t&p'}],'TimeSeries'),xyz);
% $$$ pper = logical(pper.data);



%[LUi,LRi,FSri,VTi] = erpPCA(umaps',5);
[FSri,LRi] = nnmf(umaps',5); LRi=LRi';

% START figure
figure()


% PRINT unit info
subplot2(4,8,1,1);
text(0.1,1,{['Unit: ',num2str(cluId)],...
          Trials{sesId}.filebase,...
          ['eDist:   ',num2str(nq(u).eDist)],...
          ['Refrac:  ',num2str(log10(nq(u).Refrac))],...
          ['SNR:     ',num2str(nq(u).SNR)],...
          ['AmpSym:  ',num2str(nq(u).AmpSym)],...
          ['SpkWidthR:  ',num2str(nq(u).SpkWidthR)]...
         },...
     'VerticalAlignment','top');


% PLOT auto CCG
subplot2(4,8,1,2);
bar(tbins,accg(:,u));
axis('tight');

% PLOT theta phase distribution
subplot2(4,8,1,3);
rose(phz(spkL(cluId)),36);

% PLOT head direction distribution
subplot2(4,8,1,4);
    hdBinOcc = accumarray(hdBinInd(tper),1,[36,1],@sum);
    hdBinSpk = accumarray(hdBinInd(spkX(cluId)),1,[36,1],@sum);
    polar(hdBinCtr',hdBinSpk./hdBinOcc.*xyz.sampleRate);
    title('Loc');
    %ylim([0,40]);
subplot2(4,8,1,5);
    hdBinOcc = accumarray(hdBinInd(sper),1,[36,1],@sum);
    hdBinSpk = accumarray(hdBinInd(spkR(cluId)),1,[36,1],@sum);
    polar(hdBinCtr',hdBinSpk./hdBinOcc.*xyz.sampleRate);
    title('Rear');
    %ylim([0,40]);
subplot2(4,8,1,6);
    hdBinOcc = accumarray(hdBinInd(pper),1,[36,1],@sum);
    hdBinSpk = accumarray(hdBinInd(spkP(cluId)),1,[36,1],@sum);
    polar(hdBinCtr',hdBinSpk./hdBinOcc.*xyz.sampleRate);
    title('Pause');
    %ylim([0,40]);


% PLOT eigVec
for v = 1:5,
    eigVec = nan(binDims(end-1:end));
    pev = zeros([sum(validDims),1]);
    %pev(~invalidDims) = LUi(:,v);
    pev(~invalidDims) = LRi(:,v);
    eigVec(validDims) = pev;
    subplot2(4,8,3,v+2);
    %imagescnan({bins{end-1:end},eigVec'});
    imagesc(bins{end-1},bins{end},eigVec');
    colormap('jet');
    colorbar();
    axis('xy');
end
% PLOT eigScr
%subplot2(4,8,3,8);
%plot(VTi(1:5,4),'-+');
for v = 1:5,
    subplot2(4,8,4,v+2);    
    imagesc(reshape(FSri(:,v),4,4)');
    colormap('jet');
    colorbar();
    axis('xy');
end
% PLOT spatially binned behavior space ratemaps
subplot2(4,8,[3,4],[1,2])
spcBhvMap = reshape(permute(rmapa(2:5,2:5,:,:,u),[3,1,4,2,5]),112,112);
spcBhvMap(~repmat(reshape(validDims,28,28)',[4,4])) = nan;
%spcBhvMap(~repmat(validDims,[16,1])) = nan;
imagesc(spcBhvMap');
mrate = prctile(spcBhvMap(~isnan(spcBhvMap(:))),[0.1,99.9]);
caxis([mrate]);
axis('xy');
colormap('jet')
colorbar();

% PLOT theta behavior space ratemap
subplot2(4,8,2,1)
plot(bfsi{sesId},cluId,1,'colorbar',[mrate],true,'colorMap',@jet,'mazeMask',reshape(validDims,bfsi{1}.adata.binSizes'));

% PLOT theta spatial ratemap
subplot2(4,8,2,2);
plot(pfti{sesId},cluId,1,'text',[mrate],true,'colorMap',@jet);
title('theta');

% PLOT state spatial ratemap
stid = [2,3,4,6,7];
for s = 1:5,
    subplot2(4,8,2,s+2);
    plot(pfsi{sesId}{stid(s)},cluId,1,'text',[mrate],true,'colorMap',@jet);
    title(pfsiStates{stid(s)});
end

%% END Main Plot -----------------------------------------------------


% NEED spk/cycle maps 
lowSampleRate = 5;
fxyz = copy(xyz);
fxyz = resample(fxyz,lowSampleRate);

ufr = Trial.load('ufr',fxyz,[],Trial.spk.map(:,1),1/lowSampleRate,'count');



uphz = unwrap(phz.data);
uphz = interp1(1:size(uphz,1),uphz,linspace(round(250/lowSampleRate),size(uphz,1),size(fxyz,1)));

cyc = (circshift(uphz,-1) - uphz)' ./ (2*pi);

spc = bsxfun(@rdivide,ufr.data,cyc);
vxy = fxyz.vel('hcom',[1,2]);

xbinEds = linspace(-500,500,22);
xbinCtr = round(mean([xbinEds(2:end);xbinEds(1:end-1)]));
ybinEds = linspace(-500,500,22);
ybinCtr = round(mean([ybinEds(2:end);ybinEds(1:end-1)]));



figure,
u = 166;
sts = 'w+n&t';
state = cast([stc{sts,lowSampleRate}],'TimeSeries');
state.data(end) = [];
state.data = logical(state.data);
subplot(3,3,1);
plot(log10(vxy(state.data)),ufr(state.data,u)*lowSampleRate+randn(sum(state.data),1),'.');
xlim([-1,2]);
subplot(3,3,2);
hist2([log10(vxy(state.data)),ufr(state.data,u)*lowSampleRate],linspace(-0.5,2,15),linspace(0,70,15),'xprob');
caxis([0,0.2])
subplot(3,3,4);
plot(log10(vxy(state.data)),spc(state.data,u),'.');
xlim([-1,2]);
state.data = state.data&vxy.data>0;
[P,S] = polyfit(log10(vxy(state.data)),spc(state.data,u),1);
subplot(3,3,5);
hist2([log10(vxy(state.data)),spc(state.data,u)],linspace(-0.5,2,15),linspace(0,10,15),'xprob');
caxis([0,0.2])

figure,
u = 166;
typ = 'int';
sts = 'w+n&t';
state = cast([stc{sts,lowSampleRate}],'TimeSeries');
state.data(end) = [];
state.data = logical(state.data);
xind = discretize(fxyz(state.data,'hcom',1),xbinEds);
yind = discretize(fxyz(state.data,'hcom',2),ybinEds);
subplot(221);
    outufrW = accumarray([xind,yind],ufr(state.data,u)*lowSampleRate,[numel(xbinEds)-1,numel(ybinEds)-1],@mean);
    imagesc(xbinCtr,ybinCtr,outufrW');axis('xy');
    colorbar();
    switch typ,case 'int',caxis([0,70]),otherwise,caxis([0,15]);end
    colormap('jet');
subplot(222);
    outspcW = accumarray([xind,yind],spc(state.data,u),[numel(xbinEds)-1,numel(ybinEds)-1],@mean);
    imagesc(xbinCtr,ybinCtr,outspcW');axis('xy');
    colorbar();
    switch typ,case 'int',caxis([0,7]),otherwise,caxis([0,3]);end    
    colormap('jet');
sts = 'p&t';
state = cast([stc{sts,3}],'TimeSeries');
state.data(end) = [];
state.data = logical(state.data);
xind = discretize(fxyz(state.data,'hcom',1),xbinEds);
yind = discretize(fxyz(state.data,'hcom',2),ybinEds);
subplot(223);
    outufrP = accumarray([xind,yind],ufr(state.data,u)*lowSampleRate,[numel(xbinEds)-1,numel(ybinEds)-1],@mean);
    imagesc(xbinCtr,ybinCtr,outufrP');axis('xy');
    colorbar();
    switch typ,case 'int',caxis([0,70]),otherwise,caxis([0,15]);end    
    colormap('jet');
subplot(224);
    outspcP = accumarray([xind,yind],spc(state.data,u),[numel(xbinEds)-1,numel(ybinEds)-1],@mean);
    imagesc(xbinCtr,ybinCtr,outspcP');axis('xy');
    colorbar();
    switch typ,case 'int',caxis([0,7]),otherwise,caxis([0,3]);end        
    colormap('jet');


% $$$ figure,
% $$$ subplot(211);hist(outufrP(:)-outufrW(:),100)
% $$$ subplot(212);hist((outspcP(:)-outspcW(:))*8,100)


figure
subplot(121);
ind = nniz(outufrW(:))&nniz(outspcW(:));
hist(outufrW(ind)-outspcW(ind)*8,30);
subplot(122);
ind = nniz(outufrW(:))&nniz(outspcW(:));
hist(outufrP(ind)-outspcP(ind)*8,30);

sper = [Trial.stc{'t-m-s',sampleRate}];

thpks = LocalMinima(abs(phz.data),10,0.1);
thpks = thpks(WithinRanges(thpks(:,1),sper),:);
thpks = [thpks,circshift(thpks,-1)];
thpks([1,end],:) = [];
thpks(diff(thpks,1,2)>50,:) = [];

figure();
for s = 1:5,
    subplot(1,5,s);
    hist2([1./(diff(thpks,1,2)./sampleRate),                    ... XData
           circshift(1./(diff(thpks,1,2)./sampleRate),s-3)],    ... YData
          linspace(5,13,30),                                    ... XBin edges
          linspace(5,13,30),                                    ... YBin edges
          'xyprob',                                             ... Normalization(none,(x)yprob)
          ''                                                    ... Transform (none,mud)
    );
    caxis([0,0.02]);
    line([5,13],[5,13],'Color','r');
end

%figure(); hist(diff(thpks,1,2),1000)
nui = numel(unitsInts{20});
thspkPPC = nan([size(thpks,1),nui]);
thspkN = nan([size(thpks,1),nui]);
thspkCM = nan([size(thpks,1),nui]);
thspkVel = nan([size(thpks,1),nui]);
thspkHz = nan([size(thpks,1),nui]);
thspkHvf = nan([size(thpks,1),nui]);
thspkHvl = nan([size(thpks,1),nui]);
for u = 1:nui,
    res = spk(unitsInts{20}(u));
    for t = 1:size(thpks,1),
        tres = res(WithinRanges(res,thpks(t,:)));
        thspkN(t,u) = numel(tres);
        if thspkN(t,u)>0,
            thspkHz(t,u)  = mean(xyz(tres,'hcom',3));            
            thspkVel(t,u)  = mean(vel(tres,2));
            thspkHvf(t,u)  = mean(hfl(tres,1));        
            thspkHvl(t,u)  = mean(hfl(tres,2));                
            thspkCM(t,u)  = circ_mean(phz(tres));
            thspkPPC(t,u) = PPC(phz(tres));
        end
    end
end


figure();
for u = 1:numel(unitsInts{20}),
clf();
    ind = ~isnan(thspkPPC(:,u))&~isnan(thspkN(:,u));
nmz = ''; %'xprob'
subplot2(2,3,1,1);
    hist2([log10(sqrt(sum([thspkHvf(ind,u),thspkHvl(ind,u)].^2,2))),thspkN(ind,u)],linspace(-2,2,30),0.5:1:20,nmz);
    %    hist2([log10(thspkVel(ind,u)),thspkN(ind,u)],linspace(-2,2,30),0.5:1:20,nmz);
    %caxis([0,0.2]);
    colormap('jet');
subplot2(2,3,1,2);
    hist2([(thspkHvf(ind,u)),thspkN(ind,u)],linspace(-20,80,30),0.5:1:20,nmz);
    %caxis([0,0.2]);
    colormap('jet');
subplot2(2,3,1,3);
    hist2([(thspkHvl(ind,u)),thspkN(ind,u)],linspace(-60,60,30),0.5:1:20,nmz);
    %caxis([0,0.2]);
    colormap('jet');
nmz = 'xprob';
subplot2(2,3,2,1);
%    hist2([log10(thspkVel(ind,u)),thspkN(ind,u)],linspace(-2,2,30),0.5:1:20,nmz);
    hist2([log10(sqrt(sum([thspkHvf(ind,u),thspkHvl(ind,u)].^2,2))),thspkN(ind,u)],linspace(-2,2,30),0.5:1:20,nmz);
    caxis([0,0.2]);
    colormap('jet');
subplot2(2,3,2,2);
    hist2([(thspkHvf(ind,u)),thspkN(ind,u)],linspace(-20,80,30),0.5:1:20,nmz);
    caxis([0,0.2]);
    colormap('jet');
subplot2(2,3,2,3);
    hist2([(thspkHvl(ind,u)),thspkN(ind,u)],linspace(-60,60,30),0.5:1:20,nmz);
    caxis([0,0.2]);
    colormap('jet');
    title(num2str(unitsInts{20}(u),'%d'))
waitforbuttonpress();
end


figure();
ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>3;
out = hist2([thspkPPC(ind),thspkCM(ind)],30,15);
imagesc([out,out]')

figure();
ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>2;
out = hist2([thspkN(ind),thspkCM(ind)],15,30);
imagesc([out,out]')

dthpks = diff(thpks,1,2);

figure();
for s  = -2:2,
    subplot(1,5,3+s);
    tdthpks = circshift(dthpks,s);
    ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>=2;
    out = hist2([tdthpks(ind),thspkCM(ind)],40,60);
    imagesc([out,out]');
    caxis([0,40]);
    title(num2str(u));
end
colormap('jet');


%hist2([diff(thpks,1,2),thspkPPC'],30,30);

%---------------------------------------------------------------------------------------------------

% fet
t = 20;
clearvars('-GLOBAL','AP');

pfsArgs = struct('states',           'theta-groom-sit',      ...
                 'binDims',          [0.1,0.1,0.2],          ...
                 'SmoothingWeights', [2,2,2],                ...
                 'numIter',          1,                      ...
                 'boundaryLimits',   [-2,0.8;-0.8,2;-2,2],   ...
                 'halfsample',       false);
ifs = compute_bhv_ratemaps(Trials{t},unitsInts{t},'fet_HPBPHS',[],'none',pfsArgs,1,inf,false);

ifs = cf(@(T,U)  compute_bhv_ratemaps(T,U,'fet_HPBPHS',[],'none',pfsArgs,1,inf,false), Trials,unitsInts);

[irmaps,~] = decapsulate_and_concatenate_mtaapfs({ifs},unitsInts(20));

ibins = ifs.adata.bins;
ibinDims = ifs.adata.binSizes';

% RESHAPE rmaps into [ xPosition, yPosition, headPitch, bodyPitch, unit ]
irmapa = nan([ibinDims,size(irmaps,2)]);
for u = 1:size(irmaps,2),
    irmapa(:,:,:,u) = reshape(irmaps(:,u),ibinDims);
end

figure,
u = find(45==unitsInts{20});
for s = 1:20,
    subplot(5,4,s)
    imagesc(irmapa(:,:,s,u)')
    axis('xy');
    colormap('jet')
end
ForAllSubplots('caxis([0,40]);');




%% 
cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials(1),unitsInts(1));
tag = 'interneurons_xyhb_2020';
tag = 'interneurons_xyhb_2020_final';    overwrite = false;
tag = 'interneurons_xyhb_2020_loc';      overwrite = false;
tag = 'interneurons_xyhb_2020_loc_rear'; overwrite = true;
tag = 'interneurons_xyhb_2020_theta';    overwrite = true;
pfi = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);

