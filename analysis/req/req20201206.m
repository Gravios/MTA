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
    [  4,  5],...                                                        jg04 20120211
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
    [  3,  7, 14, 16, 34],...                                            jg05 20120324
    [  5, 13, 27, 66, 71, 72, 75, 91, 94,105,127,157,209,232,233,251,254],...ER06 20130624
    [  4, 17, 20, 22, 23, 25, 29 ,30, 31, 40, 66, 69],...                Ed10 20140815
    [ 93,175] ...                                                         er01 20110722
};
    

t = 20;    
Trial = Trials{t};
spk = Trial.spk.copy();
spk.create(Trial,              ...  
           sampleRate,         ...
           '',                 ...
           unitsInts{t},       ...
           ''                  ...
);    
pft = pfs_2d_theta(Trial);

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
tag = 'interneurons_xyhb_2020_final'; overwrite = true;
pfi = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);



[rmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(pfi,unitsInts);

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

rmapa = nan([binDims,size(rmaps,2)]);
for u = 1:size(rmaps,2),
    rmapa(:,:,:,:,u) = reshape(rmaps(:,u),binDims);
end


rmap = rmapa(:,:,:,:,200);
figure();
for x = 1:6,
    for y = 1:6,
        subplot2(6,6,7-y,x);
        imagescnan({bins{3:4},sq(rmap(x,y,:,:))'});
        axis('xy');
    end
end

% bmaps = reshape(rmapa(3,4,:,:,:),[],size(rmapa,length(binDims)+1));

bmaps = [];
for k = 2:4,
    for j = 2:4,
        bmaps(:,:,k,j) = reshape(rmapa(k,j,:,:,:),[],size(rmapa,length(binDims)+1));
    end
end
bmaps = reshape(bmaps,[size(bmaps,1),prod(size(bmaps))/size(bmaps,1)]);

figure();
imagesc(reshape(validDims,binDims(end-1:end))');
axis('xy');


vmaps = bmaps(validDims,:);
vmaps(:,sum(isnan(vmaps))>0)= [];
vmaps(isnan(vmaps(:)))=0;

[LU,LR,FSr,VT] = erpPCA(vmaps',5);

figure()
for v = 1:5,
    eigVec = nan(binDims(end-1:end));
    eigVec(validDims) = LU(:,v);
    subplot(1,5,v);
    imagescnan({bins{end-1:end},eigVec'});
    axis('xy');
end



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



% fet

function [fet,featureTitles,featureDesc] = fet_HB_pitchB(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     referenceTrial'  , 'Ed05-20140529.ont.all'
%     referenceFeature', ''
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'   , Trial.xyz.sampleRate,                                       ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , {{}},                                                       ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', '',                                                         ...
                 'overwriteFlag'   , false                                                       ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,overwriteFlag] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  

if overwriteFlag,
% XYZ preprocessed 
    xyz = Trial.load('xyz');
% LOAD pitches 
    pch = fet_HB_pitch(Trial);
    if ~isempty(referenceTrial),
% MAP to reference trial
        pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
    end
    pch.resample(newSampleRate);
    xyz.resample(newSampleRate);
% CONCATENATE features
    fet.data = [circ_dist(pch(:,3),pch(:,1)),pch(:,1)];
    fet.data(~nniz(xyz),:)=0;
    fet.save();
else
    load(fet,Trial);
end


featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['upper body pitch relative to xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------