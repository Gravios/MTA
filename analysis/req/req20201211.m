%req20201211
%    Tags: interneuron pfs ego ratemaps
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: BhvPlaceCode
%    Description: Interneuron interaction with placefield and egocentric positions


Trial = MTATrial.validate('jg05-20120312.cof.all');
sampleRate = 250;

spk = copy(Trial.spk);
spk.create(Trial,sampleRate,'',[],'');

ints = select_units(Trial,'int');
apyrs = select_units(Trial,'pyr');
pyrs = select_placefields(Trial);
[accg,tbins] = autoccg(Trial);

halfBins = 50;
binSize = 1;
normalization = 'count';


accg = []
for int = ints,
i = find(ints==int);
for pyr = pyrs,
    p = find(pyr==pyrs);
    pRes = spk(pyr);
    iRes = spk(int);
    [tccg,tbin] = CCG([iRes;pRes],[i.*ones(size(iRes));p.*ones(size(pRes))],binSize,halfBins,sampleRate,[i;p],normalization);
    accg(:,p,i) = tccg(:,1,2);
end
end

figure();
for i = 1:19,
    subplot(4,5,i);
    bar(tbin,accg(:,pyrs==104,i));
    title(num2str(i));
end


[106,103]


xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
pft = pfs_2d_theta(Trial);

rot = 0.17;
hbaCorrection = -0.25;
thetaPhzChan = 70;
phzCorrection = pi/4;
headCenterCorrection = [-25,-8];
overwrite = true;

units = ints;

binPhzs = linspace(0,2*pi,6);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
pfs = cell([1,numel(binPhzc)]);


% TRANSFORM Local Field Potential -> theta phase
Trial.lfp.filename = [Trial.name,'.lfp'];
phz = load(Trial,'lfp',thetaPhzChan).phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; % mv phzCorrection -> Trial prop
phz.data(phz.data<0) = phz.data(phz.data<0) + 2*pi;
phz.data(phz.data>2*pi) = phz.data(phz.data>2*pi) - 2*pi;


hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(rot),-sin(rot);sin(rot),cos(rot)],...
                 [2,3],...
                 [1,2]);

% GET theta state behaviors, minus rear
thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20];                           % X Y
pargs.SmoothingWeights = [3, 3];                         % X Y
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-410,410;-410,410];
pargs.states       = '';
pargs.overwrite    = false;
pargs.autoSaveFlag = false;    
electrode = 0;

% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = ['_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = '';
end

for phase = 1:numel(binPhzc)
% CHECK existence of pfs object
    for p = 1:numel(pyrs)
    pargs.tag = ['egofield_int_pyr_',num2str(pyrs(p)),'_theta_phase_',num2str(phase),stag];


    filepath = fullfile(Trial.spath,...
                        [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
    
    if exist(filepath,'file'),
        pfs{p,phase} = load(filepath).Pfs;
        if overwrite,
            pfs{p,phase}.purge_savefile();
        else,
            continue;
        end;% if
    end;% if exist
    
        
    for unit = 1:numel(units),
        if unit==1 | electrode~=spk.map(spk.map(:,1)==units(unit),2), % update phase state            
            pargs.spk = copy(spk);
            electrode = 1;
            %electrode = spk.map(spk.map(:,1)==units(unit),2);
            pargs.states = copy(thetaState);
            pargs.states.label = ['thetaPhz_',num2str(phase)];
            pargs.states.data((phz(:,electrode) < binPhzs(phase) )    ...
                              | (phz(:,electrode) >= binPhzs(phase+1)) ) = 0;
            cast(pargs.states,'TimePeriods');
            resInd = WithinRanges(pargs.spk.res,pargs.states.data);
            pargs.spk.res = pargs.spk.res(resInd);
            pargs.spk.clu = pargs.spk.clu(resInd);
        end;% if
        
        [mxr,mxp] = pft.maxRate(pyrs(p));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                          [bsxfun(                                         ...
                                              @plus,                                       ...
                                              multiprod(bsxfun(@minus,                     ...
                                                               mxp,                        ...
                                                               sq(xyz(:,'hcom',[1,2]))),   ...
                                                         hvec,2,[2,3]),                    ...
                                              headCenterCorrection)],                      ...
                                          sampleRate,                                      ...
                                          'egocentric_placefield',                         ...
                                          'egopfs',                                        ...
                                          'p'                                              ...
                                          );
        pargs.xyzp = pfsCenterHR;
        pargs.units  = units(unit);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});
        if unit==1,
            try
                pfTemp.purge_savefile();
            end
            pfTemp.save();        
        end% if 
    end;% for unit
    pfTemp.save();
    pfs{p,phase} = pfTemp;    
    pfTemp = Trial;
    end;%for pyrs    
end;% for phase


figure,
for i = ints(19:end),
    for p = pyrs(20:end),
    clf();
    subplot2(2,5,1,1);
        plot(pft,p,1,'colorbar','colorMap',@jet);
        title(num2str(p));
    subplot2(2,5,1,2);
        plot(pft,i,1,'colorbar','colorMap',@jet);
        title(num2str(i));        
    subplot2(2,5,1,3);
        bar(tbin,accg(:,p==apyrs,i==ints));
        axis('tight');
    mrate = max(cell2mat(cf(@(x) x.maxRate(i,false),pfs(find(pyrs==p),1:4))));
    for phase = 1:5,
    subplot2(2,5,2,phase);
        plot(pfs{find(pyrs==p),phase},i,1,'colorbar',[mrate],false,'colorMap',@jet);
    end
    waitforbuttonpress();
    end;% for pyrs
end;%for ints



plot(pft,p,1,'colorbar','colorMap',@jet);

stc = Trial.stc;
rper = stc{'r',xyz.sampleRate};
rperOn = rper(:,1);
rperOff = rper(:,2);
hed = copy(xyz);
hed.data = hed(:,'hcom',3);
ron = hed.segs(rper(:,1)-200,400,0);
figure();
plot(ron)


set(0, 'defaultFigurePosition',  [0  1200   600    400]);

i =4;
figure,plot(pft,pyrs(i),1,'text',[],true);


res = spk(pyrs(i));
[rccg,tbins] = CCG([res;rperOn], ...
                   [ones(size(res));2.*ones(size(rperOn))],...
                   4,...
                   50, ...
                   sampleRate,...
                   [1,2],'hz');
[wccg,tbins] = CCG([res;rperOff], ...
                   [ones(size(res));2.*ones(size(rperOff))],...
                   4,...
                   50, ...
                   sampleRate,...
                   [1,2],'hz');


figure();
subplot(121);
bar(tbins,rccg(:,1,2));
subplot(122);
bar(tbins,wccg(:,1,2));


ufr = Trial.load('ufr',xyz,spk,[],1,true);

[mxr,mxy] = pft.maxRate(25);

ghz = compute_ghz(Trial,pyrs,pft,[],[],'hcom','sampleRate',sampleRate,'sigma',150);


diff(rper.data,1,2);



u = 79;


rper = stc{'r',xyz.sampleRate};
rper.data(diff(rper.data,1,2)<200,:) = [];

figure,
for o = 1:2
gsegs = sq(mean(abs(GetSegs(ghz,rper(:,o)-100,200))));
rsegs = ufr.segs(rper(:,o)-1000,2000);
[rdst,sind] = sort(gsegs(:,pyrs==u),'ascend');
srsegs = rsegs(:,sind,u);

subplot2(2,2,o,1);
imagesc(srsegs');
Lines(1000,[],'r');
subplot2(2,2,o,2);
hold('on');
plot(mean(srsegs(:,rdst<0.5)'));
plot(mean(srsegs(:,rdst<0.5)')+std(srsegs(:,rdst<0.5)'),'r');
plot(mean(srsegs(:,rdst<0.5)')-std(srsegs(:,rdst<0.5)'),'r');
end

markers = {'pelvis_root','spine_upper','hcom','nose'};
hed = copy(xyz);
hed.data = xyz(:,markers,:);
hed.model = hed.model.rb(markers);
had = create(MTADang,Trial,hed);
hba = circ_dist(had(:,1,2,1),had(:,3,4,1)+0.17)+0.42;

lfp = load(Trial,'lfp',72);
% PHZ - LFP phase within theta band 
phz = lfp.phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;

u = 119;
rper = [stc{'r',xyz.sampleRate}];
rper.data = [rper.data(:,1),rper.data(:,1)+100];
rper.data = [rper.data(:,2)-100,rper.data(:,2)];
u = 119;
res = spk(u);
res = res(WithinRanges(res,rper.data));

figure();
subplot(131);
    plot(pft,u,'mean','text');
subplot(132);
    %scatter(ghz(res,pyrs==u),phz(res),10,had(res,1,2,1),'filled');
    %caxis([-pi,pi]);
    %colormap hsv
    scatter(repmat(ghz(res,pyrs==u),[2,1]),[phz(res);phz(res)+2*pi],15,abs(repmat(hba(res),[2,1])),'filled');
    caxis([0,0.9]);
    colormap('jet');
subplot(133);
    rose(had(res,3,4,1));

    
    
res = spk(u);
[rccg,tbins,pron] = CCG([res;rperOn], ...
                   [ones(size(res));2.*ones(size(rperOn))],...
                   4,...
                   50, ...
                   sampleRate,...
                   [1,2],'hz');
[wccg,tbins] = CCG([res;rperOff], ...
                   [ones(size(res));2.*ones(size(rperOff))],...
                   4,...
                   50, ...
                   sampleRate,...
                   [1,2],'hz');

u = 119;
res = spk(u);

rper = stc{'r',xyz.sampleRate};
rper.data(diff(rper.data,1,2)<200,:) = [];
o = 1;
gsegs = sq(mean(abs(GetSegs(ghz,rper(:,o)-100,200))));
[rdst,sind] = sort(gsegs(:,pyrs==u),'ascend');
rper.data = rper(sind,:);
rper.data(rdst>0.5,:) = [];

rres = [];
tres = [];
for r = 1:size(rper,1),
    rres = [rres;res(abs(res-rper(r,1))<100)-rper(r,1)];
    tres = [tres;res(abs(res-rper(r,1))<100)];
end

figure();
hold('on');
plot(xyz(:,'spine_upper',3));
plot(res,xyz(res,'spine_upper',3),'r*');
Lines(rper(:,1),[],'m');
Lines(rper(:,2),[],'k');


exper = [73000,75000];

[mxr,mxy] = pft.maxRate(u);
figure,
for r = rper.data',
exper = r'+[-1000,1000];
subplot(121);
    cla();
    hold('on');
    plot(xyz(exper(1):exper(2),'hcom',1),xyz(exper(1):exper(2),'hcom',2));
    plot(xyz(exper(1),{'hcom','nose'},1),xyz(exper(1),{'hcom','nose'},2),'k');
    plot(xyz(exper(1)+1000,{'hcom','nose'},1),xyz(exper(1)+1000,{'hcom','nose'},2),'k');
% $$$     plot(xyz(647500,{'hcom','nose'},1),xyz(647500,{'hcom','nose'},2),'k');
% $$$     plot(xyz(647750,{'hcom','nose'},1),xyz(647750,{'hcom','nose'},2),'k');
    plot(mxy(1),mxy(2),'*r');
    xlim([-500,500]);
    ylim([-500,500]);
    circle(mxy(1),mxy(2),100);
    scatter(xyz(res(WithinRanges(res,exper)),'hcom',1),...
            xyz(res(WithinRanges(res,exper)),'hcom',2),...
            15,...
            phz(res(WithinRanges(res,exper))),...
            'filled');
    colormap('hsv');
    caxis([-pi,pi]);
subplot(122);
    cla();
    hold('on');
    plot(xyz(exper,'hcom',3));
% $$$     plot(2500,xyz(647500,'hcom',3),'*r');
% $$$     plot(2750,xyz(647750,'hcom',3),'*r');
    scatter(res(WithinRanges(res,exper))-exper(1)+1,...
            xyz(res(WithinRanges(res,exper)),'hcom',3),...
            15,...
            phz(res(WithinRanges(res,exper))),...
            'filled');
    colormap('hsv');
    caxis([-pi,pi]);
    colorbar
waitforbuttonpress();
end
    
vz = filter(copy(xyz),'ButFilter',4,1,'low')
vz.data = [0;diff(vz(:,'hcom',3))];


rres = res(WithinRanges(res,bsxfun(@plus,rper(:,2),[-200,200])));
% $$$ figure();
% $$$ hold('on');
% $$$ plot(xyz(rres,'hcom',3),...
% $$$      vz(rres));
% $$$ scatter(xyz(rres,'hcom',3),...
% $$$         vz(rres),...
% $$$         15,...
% $$$         phz(rres),...
% $$$         'filled');
% $$$ caxis([-pi,pi]);
% $$$ colorbar();
% $$$ colormap('hsv');


%rres = res(WithinRanges(res,bsxfun(@plus,rper(:,:),[-100,100])));

u = 34;
res = spk(u);

rper = stc{'r',xyz.sampleRate};
rper.data(diff(rper.data,1,2)<200,:) = [];
o = 1;
gsegs = sq(mean(abs(GetSegs(ghz,rper(:,o)-100,200))));
[rdst,sind] = sort(gsegs(:,pyrs==u),'ascend');
rper.data = rper(sind,:);
rper.data(rdst>0.5,:) = [];


figure();
% $$$ plot(xyz(rres,'hcom',3),...
% $$$         phz(rres));
%     vz(rres));
rres = res(WithinRanges(res,bsxfun(@plus,rper(:,:),[-150,150])));
subplot(121);
hold('on');
rrres = rres(vz(rres)>0.1);
scatter(xyz(rrres,'hcom',3),...
        phz(rrres),...        
        15,...
        vz(rrres),...
        'filled');
scatter(xyz(rrres,'hcom',3),...
        phz(rrres)+2*pi,...        
        15,...
        vz(rrres),...
        'filled');
caxis([-1,1]);
colorbar();
colormap('jet');
subplot(122);
hold('on');
rrres = rres(vz(rres)<-0.1);
scatter(xyz(rrres,'hcom',3),...
        phz(rrres),...        
        15,...
        vz(rrres),...
        'filled');
scatter(xyz(rrres,'hcom',3),...
        phz(rrres)+2*pi,...        
        15,...
        vz(rrres),...
        'filled');
caxis([-1,1]);
colorbar();
colormap('jet');

figure();
%rrres = rres(vz(rres)>0.5);
rrres = rres(vz(rres)<-0.5&xyz(rres,'hcom',3)>100);
%rrres = rres(xyz(rres,'hcom',3)>150);
hist(phz(rrres),linspace(-pi,pi,16));

rres = res(WithinRanges(res,bsxfun(@plus,rper.data,[-200,200])));
figure();
scatter(...
    xyz(rres,'hcom',3),...
    vz(rres),...
    15,...
    [circ_dist(phz(rres(1:end-1)),phz(rres(2:end)));0],...
    'filled');
caxis([-pi/4,pi/4]);
colorbar();
colormap('cool');


figure,
scatter(rres,phz(tres),15,abs(ghz(tres,pyrs==u)),'filled')
colormap('jet')


figure();
subplot(121);
bar(tbins,rccg(:,1,2));
subplot(122);
bar(tbins,wccg(:,1,2));

fxyz = copy(xyz);
fxyz.filter('ButFilter',4,[1],'low');
fvxy = fxyz.vel('spine_lower',[1,2]);
%fxyz.data = fxyz(:,{'pelvis_root','spine_upper'},3);
rphz = fxyz.phase([0.5,3]);




figure();
subplot(211);
    hold('on');
    plot(xyz(:,'spine_upper',3));
    plot(xyz(:,'spine_middle',3));    
    plot(xyz(:,'pelvis_root',3));        
    scatter(res,xyz(res,'pelvis_root',3),15,phz(res),'filled');    
    scatter(res,xyz(res,'spine_upper',3),15,phz(res),'filled');
    scatter(res,xyz(res,'spine_middle',3),15,phz(res),'filled');
    colormap('hsv');
    caxis([-pi,pi]);
subplot(212);
    hold('on');
    plot(abs(ghz(:,u==pyrs)));
    scatter(res,ghz(res,u==pyrs),15,phz(res),'filled');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

figure();
for u = ints;
    clf();
res = spk(u);
bper = [stc{'s&a'}];
subplot(121);
rose(rphz(res(WithinRanges(res,bper.data) ...
              & rphz(res)~=0 ...
              & fvxy(res)<0.5),1),18);
subplot(122);
rose(rphz(res(WithinRanges(res,bper.data) ...
              & rphz(res)~=0 ...
              & fvxy(res)<0.5),2),18);
title(num2str(u));
waitforbuttonpress();
end


figure();
for u = pyrs;
    clf();
    res = spk(u);
    bper = [stc{'s&a',sampleRate}];
    subplot(121);
        sphz = rphz(res(WithinRanges(res,bper.data) ...
                       & rphz(res)~=0 ...
                       & fvxy(res)<0.05),1);
        rose(sphz,18);
    subplot(122);
        sphz = rphz(res(WithinRanges(res,bper.data) ...
                       & rphz(res)~=0 ...
                       & fvxy(res)<0.05),2);
        rose(sphz,18);
% $$$     bar(linspace(-pi,3*pi,100),...
% $$$         histc([sphz;sphz+2*pi],linspace(-pi,3*pi,100)),'histc');
    title(num2str(u));
    waitforbuttonpress();
end




ang = create(MTADang,Trial,xyz);


figure,
subplot(311);hold('on')
plot(ang(stc{'p'},'bcom','hcom',2),ang(stc{'p'},'bcom','hcom',3),'.');
plot(ang(stc{'x'},'bcom','hcom',2),ang(stc{'x'},'bcom','hcom',3),'.');
subplot(312);hold('on');
plot(ang(stc{'p'},'bcom','hcom',2),ang(stc{'p'},'bcom','hcom',3),'.');
plot(ang(stc{'h'},'bcom','hcom',2),ang(stc{'h'},'bcom','hcom',3),'.');
plot(ang(stc{'w'},'bcom','hcom',2),ang(stc{'w'},'bcom','hcom',3),'.c');
subplot(313);hold('on');
plot(ang(stc{'p'},'bcom','hcom',2),ang(stc{'p'},'bcom','hcom',3),'.');
plot(ang(stc{'l'},'bcom','hcom',2),ang(stc{'l'},'bcom','hcom',3),'.');
plot(ang(stc{'w'},'bcom','hcom',2),ang(stc{'w'},'bcom','hcom',3),'.c');
linkaxes(findobj(gcf,'Type','Axes'),'xy');


figure();
hold('on');
plot(ang(stc{'w'},'bcom','hcom',2),ang(stc{'w'},'hcom','nose',2),'.c');
plot(ang([stc{'w&h'}],'bcom','hcom',2),ang([stc{'w&h'}],'hcom','nose',2),'.g');
plot(ang([stc{'w&l'}],'bcom','hcom',2),ang([stc{'w&l'}],'hcom','nose',2),'.b');




% how does the rate on the descending phase of theta depend upon the rate on the ascending phase of theta?
%ccg of theta troughs vs spikes.


vxy = vel(filter(copy(xyz),'ButFilter',4,1.5,'low'),'hcom',[1,2]);


thetaTro = LocalMinima(abs(phz.data+pi/2),0.1,10);
%bper = [stc{'w&a',sampleRate}];
bper = [stc{'x+p&a',sampleRate}];
thetaTro = thetaTro(WithinRanges(thetaTro,bper.data));
ifreq = abs(circ_dist(phz(thetaTro-1),phz(thetaTro))).*sampleRate./(2*pi);
nbins = 4;

vBinEdg = linspace(0,40,nbins);
vBinCtr = (vBinEdg(2:end)+vBinEdg(1:end-1))/2;
ivelBinInd = discretize(vxy(thetaTro),vBinEdg);

fBinEdg = linspace(5,12,nbins);
fBinCtr = (fBinEdg(2:end)+fBinEdg(1:end-1))/2;
ifreqBinInd = discretize(ifreq,fBinEdg);

hbins = 100;
u = [20,21,25,61,79,104,110,119];
u = [79];
u = pyrs;
res = spk(u);
fccg = zeros([hbins*2+1,nbins-1,nbins-1]);
for f = 1:nbins-1,
for v = 1:nbins-1,    
    ind = ifreqBinInd==f&ivelBinInd==v;
    [rccg,tbins] = CCG([res;thetaTro(ind)], ...
                       [ones(size(res));2.*ones([sum(ind),1])],...
                       1,...
                       hbins, ...
                       sampleRate,...
                       [1,2],'hz');
    fccg(:,f,v) = rccg(:,2,1);
end
end

fccg = RectFilter(fccg,3,3);

figure();
for v = 1:nbins-1,
    subplot(nbins-1,1,v);
imagesc(tbins,fBinCtr,bsxfun(@rdivide,fccg(:,:,v),sum(fccg(:,:,v)))');
axis('xy');
colormap('jet');
caxis([0,0.01]);
    hold('on');
for f = 1:nbins-1,
    plot(tbins,-cos((2*pi*fBinCtr(f))*(tbins./1000))+fBinCtr(f));
end
end


Lines(tbins(LocalMinima(-fccg(:,3),20)),[],'r');

figure,plot(fccg(LocalMinima(-fccg(:,3),20),3))
