% req20170921 ---------------------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: lateral motions of head and body synchrony during locomotion
%  Bugs: NA


Trial = MTATrial.validate('Ed05-20140529.ont.all');
stc = Trial.load('stc','hand_labeled_rev1_Ed');
xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQI');
features = fet_bref(Trial);
features.data = features.data(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
fetHB = fet_HB_angvel(Trial);



fxyz = xyz.copy();
fxyz.filter('ButFilter',3,2.4,'low');



hbfet = features.copy();
hbfet.data = unity([hbfet.data(:,[21]),fetHB(:,2)]);
hbfet.data = hbfet.data(:,[21,25]);



specParm =    struct('nFFT',2^7,'Fs',features.sampleRate,...
                           'WinLength',2^6,'nOverlap',2^6*.875,...
                           'FreqRange',[1,10]);


[ys,fs,ts] = fet_spec(Trial,hbfet,'mtcsdglong',false,'defspec',specParm);

wper = stc{'w',ys.sampleRate};
vxy = fxyz.vel([1,2],[1,2]);
vxy.resample(ys);
vxy.data(vxy.data<1e-3) = nan;
%vxy.data = log10(vxy.data);


chr = [];
mvx = [];
mys = [];
man = [];
for w = wper.data',
    chr(end+1,:) = nanmean(abs(ys(w,:,1,2)))./nanmean(sqrt(ys(w,:,1,1).*ys(w,:,2,2)));
    mvx(end+1,:) = mean(vxy(w,:));
    mys(end+1,:) = nanmean(ys(w,:,2,2));
    %man(end+1,:) = circ_mean(fetHBp(w,:));
end


% FILTER by walk duration
wdur = diff(wper.data,1,2);
ind = wdur>=14;

out = [];
for f = 1:numel(fs),
    out(:,f) = histc(chr(ind,f),linspace(0,1,30));
end
figure,imagesc(1:30,fs,out');


figure,hist(chr(ind,5),linspace(0,1,40));
figure,plot(chr(ind,4),mvx(ind,2),'.');
figure,plot(chr(ind,5),log10(mys(ind,15))./log10(mys(ind,5)),'.');
figure,plot(chr(ind,5),man(ind,2),'.');


%% Multi Session 

Trials   = af(@(t)  MTATrial.validate(t),  get_session_list('hand_labeled'));
stc      = cf(@(t)  t.load('stc'),         Trials);
xyz      = cf(@(t)  preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI'),  Trials);
features = cf(@(t)  fet_bref(t),         Trials);
           cf(@(f)  set(f,'data',f(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30])),features);

fxyz     = cf(@(x)  x.copy(),  xyz);
           cf(@(x)  x.filter('ButFilter',3,2.4,'low'),  fxyz);

hbfet    = cf(@(f)  f.copy(),                           features);
           cf(@(f)  set(f,'data',f.data(:,[21,25])),    hbfet);
fetp     = cf(@(t)  fet_HB_pitch(t),         Trials);
           cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'), fetp, Trials);
           
specParm = repmat({struct('nFFT',2^7,'Fs',features{1}.sampleRate,...
                         'WinLength',2^6,'nOverlap',2^6*.875,...
                         'FreqRange',[1,10])},...
                  [1,numel(Trials)]);


[ys,fs,ts] = cf(@(t,f,p)  fet_spec(t,f,'mtcsdglong',false,'defspec',p),  Trials,hbfet,specParm);

vxy      = cf(@(x)    x.vel([1,2],[1,2]),   fxyz);
wper     = cf(@(s,y)  s{'w',y.sampleRate},  stc,  ys);
           cf(@(v,y)  v.resample(y),        vxy,  ys);
           cf(@(f,y)  f.resample(y),        fetp, ys);           
for s = 1:numel(vxy), vxy{s}.data(vxy{s}.data<1e-3) = 1e-3; end
ang      = cf(@(t,x)  create(MTADang,t,x),  Trials, fxyz);

chr = [];
mvx = [];
wdur = [];
man = [];
for s = 1:numel(ys),
    for w = wper{s}.data',
        wdur(end+1) = diff(w);        
        w = round(mean(w));
        chr(end+1,:) = (abs(ys{s}(w,:,1,2)))./(sqrt(ys{s}(w,:,1,1).*ys{s}(w,:,2,2)));
        mvx(end+1,:) = (vxy{s}(w,:));
        %man(end+1,:) = (ang{s}(w,5,7,2));
        man(end+1,:) = (fetp{s}(w,2));
        %chr(end+1,:) = nanmean(abs(ys{s}(w,:,1,2)))./nanmean(sqrt(ys{s}(w,:,1,1).*ys{s}(w,:,2,2)));
        %mvx(end+1,:) = nanmean(vxy{s}(w,:));
    end
end

ind = wdur>14;
out = [];
for f = 1:numel(fs),
    out(:,f) = histc(chr(ind,f),linspace(0,1,30));
end


figure,imagesc(1:30,fs{1},out')
figure,plot(chr(ind,3),mvx(ind,2),'.');
figure,plot(chr(ind,5),man(ind,1),'.');


