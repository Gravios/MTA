
% drift tracking using the phase difference between channels in
% the theta band

%%%<<< Load Data ---------------------------------------------------------------
ThetaRC_load_data();
% - sampleRate
% - xyz
% - units
% ... etc

convolution = @conv;
segment = @MTAData.segs;

trialIndex = 19;

Trial   = Trials     { trialIndex };
units   = Units      { trialIndex };
meta    = sessionList( trialIndex );
pft = pfs_2d_theta(Trial,units);
stc = Trial.stc.copy();

txyz = preproc_xyz(Trial,'trb',sampleRate);

hfvxy = copy(txyz);
hfvxy.filter('ButFilter',4,2.4,'low');
hfvxy = hfvxy.vel('hcom',[1,2]);

lfp_hpc = Trial.load('lfp', [65:96]);

figure();
plot(bsxfun(@minus, lfp_hpc.data(1:1e5,:)./10, linspace(0,16000,size(lfp_hpc,2))))


lfp_shank = Trial.load('lfp', 55:64);

phz_shank = lfp_shank.phase();

figure();
hold('on');
plot(circ_dist(phz_shank(:,1),phz_shank(:,8)))

pdiff = [];
for c1 = 1:8
    for c2 = (c1+1):8
        pdiff(:,c1,c2) = circ_dist(phz_shank(:,c1),phz_shank(:,c2));
    end
end

figure,imagesc(pdiff')

figure,
hold on
imagesc([mean(pdiff(:,1,2),2),...
 mean(pdiff(:,2,3),2),...
mean(pdiff(:,3,4),2),...
mean(pdiff(:,4,5),2),...
mean(pdiff(:,5,6),2),...
mean(pdiff(:,6,7),2),...
mean(pdiff(:,7,8),2)]')


figure,
hold on
imagesc([mean(pdiff(:,1,4),2),...
 mean(pdiff(:,2,5),2),...
mean(pdiff(:,4,7),2),...
mean(pdiff(:,5,8),2)]')

figure,
plot([mean(pdiff(:,1,4),2),...
 mean(pdiff(:,2,5),2)])

figure
plot([mean(pdiff(:,4,7),2),...
mean(pdiff(:,5,8),2)]);

figure,
plot(mean(pdiff(:,2,6),2)-mean(pdiff(:,5,8),2));
ylim([0,2.2])
Lines([],0,'r');


lfp_trc = Trial.load('lfp', meta.subject.channelGroup.thetarc);
lfp_trc.data = diff(lfp_trc.data, 1, 2);

lfp_pyr = Trial.load('lfp', 71);
lfp_prc = Trial.load('lfp', [70,72]);
lfp_prc.data = diff(lfp_prc.data, 1, 2);

lfp_rad = Trial.load('lfp', [76]);
lfp_rrc = Trial.load('lfp', [75,77]);
lfp_rrc.data = diff(lfp_rrc.data, 1, 2);

lfp_lcm = Trial.load('lfp', [84]);
lfp_lrc = Trial.load('lfp', [83,86]);
lfp_lrc.data = diff(lfp_lrc.data, 1, 2);

lfp_dgr = Trial.load('lfp', [87]);
lfp_drc = Trial.load('lfp', [86,88]);
lfp_drc.data = diff(lfp_drc.data, 1, 2);

[70,71,72,75,77,82,85,86,88]

figure,hold('on');
plot(lfp_pyr.data/3+4000,'c');
plot(lfp_rad.data/6+2000,'r');
plot(lfp_lcm.data/10,'g');
plot(lfp_dgr.data/10-2000,'m');
plot(lfp_prc.data/2-4000,'c');
plot(lfp_rrc.data/3-6000,'r');
plot(lfp_lrc.data/3-9000,'g');
plot(lfp_drc.data/3-12000,'m');
xlim([3586946,3591439]);
ylim([-23790,9562]);

flfp_prc = copy(lfp_prc);
flfp_prc.data = flfp_prc.data./3000;
flfp_prc.filter('ButFilter',4,4,'high');


template_index = 3661465;
template_halfwindow = 128;

ti  = template_index;
thw = template_halfwindow;
template_indices = [ ti-(thw-1) : ti+thw ];
template_swr = flfp_prc.data( template_indices );
template_length = thw * 2 + 1;

cprc =                          ...
    convolution(                ...
        flfp_prc.data,          ...
        template_swr,           ...
        'same'                  ...
);
gwin = gausswin(128);
gwin = gwin./sum(gwin);
pprc =                          ...
    convolution(                ...
        abs(cprc),              ...
        gwin,                   ...
        'same'                  ...
);


figure();
hold('on');
plot( flfp_prc.data );
plot( pprc./10 );


pinds = ThreshCross(pprc,40,50);
pinds = round(mean(pinds,2));

pinds = LocalMinima( -pprc, 100, -40);

pin = 4;
ind = [pinds(pin)-105:pinds(pin)+105]';
figure();
hold('on');
plot(lfp_pyr(ind),'k');
plot(lfp_prc(ind),'c');
plot(cprc(ind)./10e4,'r'); 
plot(pprc(ind)./10e4,'r'); 

thw = template_halfwindow;
tlen = template_length;
template_dst =                 ... ASSINGMENT
    mean(                      ... mean
        abs(                   ... absolute value
            segment(           ... frist dimension
                flfp_prc,      ... tensor
                pinds - thw,   ... segment start
                tlen,          ... segment length
                0              ... overrun fill
            )                  ...
        )                      ...
    );
template_var =                 ... ASSINGMENT
    var(                       ... variance
        abs(                   ... absolute value
            segment(           ... frist dimension
                flfp_prc,      ... tensor
                pinds - thw,   ... segment start
                tlen,          ... segment length
                0              ... overrun fill
            )                  ...
        )                      ...
    );

figure();
    plot(                      ...
        log10(template_dst),   ...
        log10(template_var),   ...
        '.'                    ...
    );


mind = find( log10(template_dst) < -0.35 , 10, 'first');


% sthw =  template_halfwindow;
mm = 6;
pind = pinds(mind(mm));
ind = [pind-sthw+1:pind+sthw]';
figure();
hold('on');
plot(flfp_prc(ind))
plot(template_swr);

ind = [pind-2*sthw+1:pind+2*sthw]';
tshifts = [2*sthw-sthw+1]:[2*sthw+sthw];
figure();
    hold('on');
    plot(flfp_prc(ind));
    plot(                      ...
        tshifts,               ...
        template_swr           ...
    )

template_ccg =                 ...
    convolution(               ...
        flfp_prc(ind),         ...
        template_swr,          ...
        'same'                 ...
);

figure();
    plot(                      ...
         template_ccg          ...
    );

template_ccg_threshold = -50     
template_ccg_proximity = 100;    % samples

template_ccg_peaks =           ...
    LocalMinima(               ...
       -template_ccg,          ...
        template_ccg_proximity,...
        template_ccg_threshold ...
);
        
template_shifts =              ...
    minus(                     ...
        template_ccg_peaks,    ...
        template_length        ...
    );


figure();
    hold('on');
    plot(                      ...
        flfp_prc(ind)          ...
    );
    ind =                      ...
        plus(                  ...
            [pind-2:pind+210]', ...
            template_shifts(1) ...
        );

    plot(                      ...
        plus(                  ... x values
            tshifts,           ...
            template_shifts(1) ...
        ),                     ...
        template_swr           ... y values
        );
    
    plot(                      ...
        plus(                  ... x values
            tshifts,           ...
            template_shifts(2) ...
        ),                     ...
        template_swr           ... y values
    );

template_reconstrution =       ...
    nan(                       ...
        size(template_ccg)     ...
    );

template_reconstrution_count = ...
    zeros(                     ...
        [size(template_ccg,1), ...
         numel(template_shifts)] ...
);

template_rec_ind = [];
for tsh = 1 : numel(template_shifts),
    template_rec_ind(:,tsh) = tshifts + template_shifts(tsh);
    template_reconstrution(template_rec_ind(:,tsh)) = template_swr;
    template_reconstrution_count(template_reconstrution_index(:,tsh)) = 1;
end
template_reconstrution = template_reconstrution ./ sum(template_reconstrution_count,2);

ind = [pind-sthw*2+1:pind+sthw*2]';
figure();
    hold('on');
    plot(                      ...
        flfp_prc(ind)          ...
    );
    plot(                      ...
        template_reconstrution ...
    );


ind = [pinds(mind(mm))-105:pinds(mind(mm))+105]';
sum(abs(template_swr-flfp_prc(ind)))
ind = [pinds(mind(mm))-105:pinds(mind(mm))+105]' + template_shifts(1);
sum(abs(template_swr-flfp_prc(ind)))
ind = [pinds(mind(mm))-105:pinds(mind(mm))+105]' + template_shifts(2);
sum(abs(template_swr-flfp_prc(ind)))

figure();




ind = [pinds(mind)-105:pinds(mind)+105]';
figure,plot(abs(template_swr-flfp_prc(ind)))

figure,plot(lfp_prc.data(3661360:3661570))

lfp_prc = Trial.load('lfp', [67:96]);
lfp_prc = Trial.load('lfp', [57,64]);
figure,plot(bsxfun(@minus, lfp_drc.data(1:1e5,:)./10, linspace(0,16000,size(lfp_drc,2))))


elfp = copy(lfp_pyr);
elfp.data = [lfp_pyr.data,... 1
             lfp_rad.data,... 2
             lfp_lcm.data,... 3
             lfp_dgr.data,... 4
             lfp_prc.data,... 5
             lfp_rrc.data,... 6
             lfp_lrc.data,... 7
             lfp_drc.data];
             
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  elfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[0,25]);
[lys,lfs,lts] = fet_spec(Trial,elfp,'mtcsdglong',false,[],specArgsTheta,[],true);

fhfvxy = copy(hfvxy);
fhfvxy.resample(lys);

hvel = log10(clip(fhfvxy.data,0.001,150));


figure,
for c = 1:9
    subplot(9,1,c);
    imagesc(lts,lfs,lys(:,:,c,c)');
end
ForAllSubplots([...
    'set(gca,''ColorScale'',''log'');'...
    'ylim([0,25]);'...
    'colormap(''jet'');'...
    'axis(''xy'');'...
]);
linkxy();


c = 2;
tdr = log10(mean(lys(:,6<lfs&lfs<12,c,c),2)) ...
      - log10(mean(lys(:,(lfs<4),c,c),2));

figure,
binscatter(log10(clip(fhfvxy.data,0.001,150)), tdr)
ForAllSubplots([...
    'set(gca,''ColorScale'',''log'');'...
    'colormap(''jet'');'...
]);


c = 2;
tdr = log10(mean(lys(:,6<lfs&lfs<12,c,c),2)) ...
      - log10(mean(lys(:,(lfs<4)|(12<lfs&lfs<14),c,c),2));


figure
subplot(121);
ind = hvel < 0 & nniz(lys) & hvel > -2.6;
binscatter(angle(lys(ind,12,3,5)), ...
           tdr(ind),...
           30);
subplot(122);
ind = hvel > 0 & nniz(lys);
binscatter(angle(lys(ind,12,3,5)), ...
           tdr(ind),...
           30);
ForAllSubplots([...
    'set(gca,''ColorScale'',''log'');'...
    'xlim([-pi,pi]);'...
    'ylim([-1.5,2]);'...
    'colormap(''jet'');'...
]);


figure
subplot(131);
ind = [stc{'a-t',lys.sampleRate}];
ind.cast('TimeSeries');
ind = logical(ind.data);
binscatter(angle(lys(ind,12,3,5)), ...
           tdr(ind),...
           30);
subplot(132);
ind = [stc{'s&t',lys.sampleRate}];
ind.cast('TimeSeries');
ind = logical(ind.data);
binscatter(angle(lys(ind,12,3,5)), ...
           tdr(ind),...
           30);
subplot(133);
ind = [stc{'a-s-r&t',lys.sampleRate}];
ind.cast('TimeSeries');
ind = logical(ind.data);
binscatter(angle(lys(ind,12,3,5)), ...
           tdr(ind),...
           30);
ForAllSubplots([...
    'set(gca,''ColorScale'',''log'');'...
    'xlim([-pi,pi]);'...
    'ylim([-1.5,2]);'...
    'grid(''on'');'...
    'colormap(''jet'');'...
]);



