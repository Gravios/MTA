
MjgER2016_load_data();

for tind = 17:23; % jg05-20120312.cof.all

    Trial = Trials{tind};
    unitSubset = units{tind};

    sampleRate = 250;
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    lfp = load(Trial,'lfp',sessionList(tind).thetaRef);
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    


    pfTemp = Trial;
    xyzp = copy(xyz);
    xyzp.data = sq(xyz(:,'hcom',[1,2]));
    thetaState = resample(cast([Trial.stc{'theta-groom-sit'}],'TimeSeries'),xyz);
    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.units        = unitSubset;
    pargs.tag          = 'thetaEXPec3';
    pargs.halfsample   = true;
    pargs.numIter      = 1001;
    pargs.states       = state;
    pargs.xyzp         = xyzp;
    pargs.overwrite    = true;
    pargs.autoSaveFlag = false;
    electrode = 0;
    for u = 1:numel(unitSubset);
        if u==1 | electrode~=spk.map(spk.map(:,1)==unitSubset(u),2), % update phase state
            electrode = spk.map(spk.map(:,1)==unitSubset(u),2);
            state = copy(thetaState);
            state.label = 'thetaEXPec3';
            state.data( phz(:,electrode)<-2.5 | phz(:,electrode)>-1) = 0;
            state = cast(state,'TimePeriods');
        end
        pargs.units  = unitSubset(u);
        pargs.states = state;
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        if u==1,
            pfTemp.purge_savefile();
            pfTemp.save();        
        end    
    end
    pfTemp.save();        
end



for tind = 19:20; % jg05-20120312.cof.all

    Trial = Trials{tind};
    unitSubset = units{tind};

    sampleRate = 250;
    xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
    lfp = load(Trial,'lfp',sessionList(tind).thetaRef);
    phz = lfp.phase([6,12]);
    phz.data = unwrap(phz.data);
    phz.resample(xyz);    
    phz.data = mod(phz.data+pi,2*pi)-pi;
    lfp.resample(xyz);    

    pfTemp = Trial;
    xyzp = copy(xyz);
    xyzp.data = sq(xyz(:,'hcom',[1,2]));
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
    pargs = get_default_args('MjgER2016','MTAApfs','struct');        
    pargs.units        = unitSubset;
    pargs.tag          = 'CA1thetaCA3inputPhase';
    pargs.halfsample   = true;
    pargs.numIter      = 1001;
    pargs.states       = '';
    pargs.xyzp         = xyzp;
    pargs.overwrite    = true;
    pargs.autoSaveFlag = false;
    electrode = 0;
    
    spk = Trial.spk.copy;
    spk.create(Trial,sampleRate,[Trial.stc{'theta-groom-sit-rear'}],unitSubset,'deburst');
    
    for u = 1:numel(unitSubset);
        if u==1 | electrode~=spk.map(spk.map(:,1)==unitSubset(u),2), % update phase state
            electrode = spk.map(spk.map(:,1)==unitSubset(u),2);
            pargs.spk = copy(spk);            
            pargs.states = copy(thetaState);
            pargs.states.label =  'CA1thetaCA3inputPhase';
            pargs.states.data( phz(:,electrode)<0.5 | phz(:,electrode)>2.8) = 0;
            cast(pargs.states,'TimePeriods');
            resInd = WithinRanges(pargs.spk.res,pargs.states.data);
            pargs.spk.res = pargs.spk.res(resInd);
            pargs.spk.clu = pargs.spk.clu(resInd);
        end
        pargs.units  = unitSubset(u);
        pfsArgs = struct2varargin(pargs);
        pfTemp = MTAApfs(pfTemp,pfsArgs{:});            
        if u==1,
            pfTemp.purge_savefile();
            pfTemp.save();        
        end    
    end
    pfTemp.save();        
end


pfTemp = MTAApfs(Trial,[],[],[],'thetaEXPec3');

xyzp = copy(xyz);
xyzp.data = sq(xyz(:,'nose',[1,2]));

state = resample(cast([Trial.stc{'theta-groom-sit'}],'TimeSeries'),xyz);
state.label = 'thetaEXPec3';
state.data(-2.5>phz.data|phz.data>-1) = 0;
state = cast(state,'TimePeriods');

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units  = unitSubset;
pargs.tag    = 'thetaEXPec3';
pargs.halfsample = false;
pargs.numIter = 1;
pargs.states = state;
pargs.xyzp   = xyzp;
pargs.overwrite = false;

pfsArgs = struct2varargin(pargs);
pfsEC3 = MTAApfs(Trial,pfsArgs{:});

xyzp = copy(xyz);
xyzp.data = sq(xyz(:,'hcom',[1,2]));
pargs.xyzp   = xyzp;
pargs.tag    = 'thetaEXPec3HCOM';
pargs.trackingMarker = 'hcom';
pfsArgs = struct2varargin(pargs);
pfsEC3hcom = MTAApfs(Trial,pfsArgs{:});


state = resample(cast([Trial.stc{'theta-groom-sit'}],'TimeSeries'),xyz);
state.label = 'thetaEXPca3';
state.data(phz.data<0|phz.data>2.5) = 0;
state = cast(state,'TimePeriods');

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units  = unitSubset;
pargs.tag    = 'thetaEXPca3';
pargs.states = state;
pargs.halfsample = false;
pargs.numIter = 1;
pargs.xyzp   = xyzp;
pargs.overwrite = false;

pfsArgs = struct2varargin(pargs);
pfsCA3 = MTAApfs(Trial,pfsArgs{:});

xyzp = copy(xyz);
xyzp.data = sq(xyz(:,'hcom',[1,2]));
pargs.xyzp   = xyzp;
pargs.tag    = 'thetaEXPca3HCOM';
pargs.trackingMarker = 'hcom';
pfsArgs = struct2varargin(pargs);
pfsCA3hcom = MTAApfs(Trial,pfsArgs{:});

end


tind = 20;
Trial = Trials{tind};
unitSubset = units{tind};

pfsEC3 = MTAApfs(Trial,[],[],[],'thetaEXPec3');
pfsCA3 = MTAApfs(Trial,[],[],[],'thetaEXPca3');
pfsEC3hcom = MTAApfs(Trial,[],[],[],'thetaEXPec3HCOM');
pfsCA3hcom = MTAApfs(Trial,[],[],[],'thetaEXPca3HCOM');


pft = pfs_2d_theta(Trial,unitSubset);

figure,
for u = unitSubset    

    subplot(332);
    plot(pfsEC3hcom,u,'mean',true,[],true);
    title(num2str(u));
    title('EC3 Input');    
    subplot(333);
    plot(pfsCA3hcom,u,'mean',true,[],true);
    title('CA3 Input');

    subplot(334);
    plot(pft,u,'mean',true,[],true);
    title('theta state')

    subplot(335);
    plot(pfsEC3,u,'mean',true,[],true);
    title(num2str(u));
    subplot(336);
    plot(pfsCA3,u,'mean',true,[],true);
    title(num2str(u));

    subplot(338);
    imagesc(pft.adata.bins{:},(plot(pfsEC3,u,'mean',true,[],true)-plot(pft,u,'mean',true,[],true))');
    caxis([repmat(max(abs(caxis)),[1,2]).*[-1,1]]);
    colorbar
    axis('xy');
    subplot(339);
    imagesc(pft.adata.bins{:},(plot(pfsCA3,u,'mean',true,[],true)-plot(pft,u,'mean',true,[],true))');
    caxis([repmat(max(abs(caxis)),[1,2]).*[-1,1]]);
    colorbar
    axis('xy');

    waitforbuttonpress();
end


[hrzt,~,hrangt] = compute_hrz(Trial,unitSubset,pft,'sampleRate',sampleRate);
[drzt,~,drangt] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
[ddzt]          = compute_ddz(Trial,unitSubset,pft,'sampleRate',sampleRate);

[hrz,~,hrang] = compute_hrz(Trial,unitSubset,pfsCA3,'sampleRate',sampleRate);
[drz,~,drang] = compute_drz(Trial,unitSubset,pfsCA3,'sampleRate',sampleRate);
[ddz]         = compute_ddz(Trial,unitSubset,pfsCA3,'sampleRate',sampleRate);

[hrzEC3,~,hrangEC3] = compute_hrz(Trial,unitSubset,pfsEC3,'sampleRate',sampleRate);
[drzEC3,~,drangEC3] = compute_drz(Trial,unitSubset,pfsEC3,'sampleRate',sampleRate);
[ddzEC3]            = compute_ddz(Trial,unitSubset,pfsEC3,'sampleRate',sampleRate);

spk = Trial.load('spk',sampleRate,'',unitSubset,'deburst');

% COMPUTE polar coordinates of horizontal position
mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                atan2(diff(xyz(:,{'hcom','nose'},2),1,2),diff(xyz(:,{'hcom','nose'},1),1,2)));
thresholds.mazeCenterDist = 380;
thresholds.mazeCenterAng = pi/2;
cind =    mazeCenterDist < thresholds.mazeCenterDist  ...
          |  abs(mazeCenterAng) < thresholds.mazeCenterAng;

stc = Trial.load('stc','msnn_ppsvd_raux');

figure();

for u = unitSubset    
    clf();
    
    res = spk(u);
    res = res(WithinRanges(res,get([stc{'x+p&t',sampleRate}],'data'))&cind(res));
    res(abs(ddz(res,u==unitSubset))>300) = [];

    
    subplot(331);
    plot(pft,u,'mean',true,[],true);

    subplot(332);
    plot(pfsEC3,u,'mean',true,[],true);
    title(num2str(u));
    subplot(333);
    plot(pfsCA3,u,'mean',true,[],true);
    title(num2str(u));

    subplot(334);
    hold('on');
    plot(drzt(res,u==unitSubset),phz(res),'.b');
    plot(drzt(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    subplot(335);
    hold('on');
    plot(drzEC3(res,u==unitSubset),phz(res),'.b');
    plot(drzEC3(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    subplot(336);
    hold('on');
    plot(drz(res,u==unitSubset),phz(res),'.b');
    plot(drz(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    subplot(337);
    hold('on');
    plot(hrzt(res,u==unitSubset),phz(res),'.b');
    plot(hrzt(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    subplot(338);
    hold('on');
    plot(hrzEC3(res,u==unitSubset),phz(res),'.b');
    plot(hrzEC3(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    subplot(339);
    hold('on');
    plot(hrz(res,u==unitSubset),phz(res),'.b');
    plot(hrz(res,u==unitSubset),phz(res)+2*pi,'.b');    
    xlim([-1,1]);
    ylim([-pi,3*pi]);

    waitforbuttonpress();
end

figure,
hold('on');
scatter(abs(ddz(res,u==unitSubset)).*sign(hrz(res,u==unitSubset)),...
        phz(res),10,hrang(res,u==unitSubset),'filled');
scatter(abs(ddz(res,u==unitSubset)).*sign(hrz(res,u==unitSubset)),...
        phz(res)+2*pi,10,hrang(res,u==unitSubset),'filled');
colormap(gca,'hsv');




a = 2.5;
figure,plot(abs(ddzt(res,u==unitSubset)).*sum(bsxfun(@times,[1,-1],[-pi/a>hrangt(res)|hrangt(res)>pi/a,-pi/a<hrangt(res)&hrangt(res)<pi/a]),2),phz(res),'.b');

xyzp = copy(xyz);
xyzp.data = sq(xyz(:,'nose',[1,2,3]));

state = resample(cast([Trial.stc{'theta-groom-sit'}],'TimeSeries'),xyz);
state.label = 'thetaEXPec3';
state.data(-2.5<phz.data|phz.data>-1) = 0;
state = cast(state,'TimePeriods');

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units  = unitSubset;
pargs.tag    = 'thetaEXPec3';
pargs.halfsample = false;
pargs.numIter = 1;
pargs.states = state;
pargs.xyzp   = xyzp;
pargs.overwrite = true;

pfsArgs = struct2varargin(pargs);
pfsEC3 = MTAApfs(Trial,pfsArgs{:});
