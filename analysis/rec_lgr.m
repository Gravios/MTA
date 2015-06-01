


%Sessions = SessionList('test_grp',...
%                '/storage/gravio/data/processed/xyz/',...
%                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true;
states = {'walk','rear','sit','turn','shake','groom'};
init_ns = numel(states);

Trial = MTATrial('jg05-20120317','all','cof');    
fet = fet_lgr(Trial);
fet.resample(30);
model_names = {};
State_cat_order = {};%repmat({''},[init_ns,1]);
lgrm = {};
pB = repmat({[]},[init_ns,1]);

 
for s = 1:init_ns,
    disp(['inter: ' num2str(s) ', finding best state'])
    temp_states =  states(cellfun(@isempty,regexp(states,['(',strjoin(State_cat_order,'|'),')'])));
    
    if ~isempty(State_cat_order),
        sws = ['-' strjoin(State_cat_order,'-')];
    else
        sws = '';
    end
    
    for i = 1:numel(temp_states),
        model_names(s,i) = {[Trial.filebase,'-','pop_lgr-' temp_states{i} sws]};%mfilename]};
        bhv_lgr(Trial,train,[temp_states(i),Trial.stc{['a-' temp_states{i} sws]}.label],...
                fet,model_names{s,i},false,false);
    end

    
    for i = 1:numel(temp_states),
        lgrm{s,i} = load([model_names{s,i} '-fet_lgr-model.mat']);
        pB{s} = cat(2,pB{s},lgrm{s,i}.B);
    end
    
    [~,BestStateInd] = max(max(abs(pB{s}(2:end,:))));
    State_cat_order{s} = temp_states{BestStateInd};
    
end

figure,
imagesc(abs(pB{1}'));
colormap jet;

tstates = states;
tstates(2) =[];
d_prime = zeros([fet.size(2),numel(tstates)]);
for i = 1:numel(tstates)
    d_prime(:,i) = ((mean(fet(Trial.stc{tstates{i}},:))-mean(fet(Trial.stc{['gper-rear-' tstates{i}]},:)))...
        ./sqrt(.5*(var(fet(Trial.stc{tstates{i}},:))+var(fet(Trial.stc{['gper-rear-' tstates{i}]},:)))))';
end

wfet = MTADxyz('data',fet(Trial.stc{'w'},d_prime(:,1)>1),'sampleRate',fet.sampleRate);
awfet = MTADxyz('data',fet(Trial.stc{'a-w'},d_prime(:,1)>1),'sampleRate',fet.sampleRate);
figure,hist(wfet(randi(wfet.size(1),20000),1),10000);

ind = randi(wfet.size(1),[40000,1]);
figure,
subplot(211),bar(linspace(-.75,2,1000),histc(mean(wfet(ind,[1]),2),linspace(-.75,2,1000)),'histc');
subplot(212),bar(linspace(-.75,2,1000),histc(mean(awfet(ind,[1]),2),linspace(-.75,2,1000)),'histc');



%lgrm = cat(1,cell2mat(lgrm));
%figure
%figure,imagesc(lgrm(1).stats.coeffcorr'),colorbar,colormap jet

Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120317');

dsx = Trial.load('xyz');

vxy = dsx.vel([1,7],[1,2]);
vxy.filter('ButFilter',3,10,'low');


vl = dsx.vel([1,2,7],[1,2]);
%vl.filter('ButFilter',3,.75,'low');
vl.data = cat(1,diff(vxy.data),zeros([1,vxy.size(2)]));
vl.filter('ButFilter',3,.6,'low');

%vz = dsx.vel(7,3);
vz = dsx.copy;
vz.data = [0;diff(dsx(:,7,3))];
vz.filter('ButFilter',3,.6,'low');

dsx.filter('ButFilter',3,.5,'low');
ang = create(MTADang,Trial,dsx);
tfet = abs(diff(circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1))));

% figure,hold on,plot(tfet,'m')
% Lines(Trial.stc{'n'}(:),[],'g')
% Lines(Trial.stc{'r'}(:),[],'r')
% Lines(Trial.stc{'w'}(:),[],'b')

%% Detect Change Points
tpnts = LocalMinima(-tfet,round(.20*dsx.sampleRate),-7e-6);
tprs = [tpnts,circshift(tpnts,-1)];
tprs(end,:) = [];

vpnts = LocalMinima(-abs(vl(:,1)),round(.20*dsx.sampleRate),-.08);
vprs = [vpnts,circshift(vpnts,-1)];
vprs(end,:) = [];

rpnts = LocalMinima(-abs(vz(:,1)),round(.20*dsx.sampleRate),-2);
rprs = [rpnts,circshift(rpnts,-1)];
rprs(end,:) = [];


tsc = [];
for i = tprs'
    tsc(end+1) = abs(circ_dist(ang(i(1),1,4,1),ang(i(2),1,4,1)))/diff(i);
end
ctsc = [];
for i = Trial.stc{'n'}.data'
    ctsc(end+1) = abs(circ_dist(ang(i(1),1,4,1),ang(i(2),1,4,1)))/diff(i);
end
edg = linspace(-8,-1,100);
figure,hold on
bar(edg,histc(log10(tsc),edg),'histc');
h = bar(edg,histc(log10(ctsc),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


tsc = [];
vsc = [];
ssc = [];
hsc = [];
asc = [];
for i = vprs'
    tsc(end+1) = mean(circ_dist(ang(i(1):i(2),1,4,1),ang(i(1)-1:i(2)-1,1,4,1)));
    vsc(end+1) = median(vxy(i(1):i(2),1));
    ssc(end+1) = std(vxy(i(1):i(2),1));
    hsc(end+1) = mean(dsx(i(1):i(2),1,3));
    asc(end+1) = mean(ang(i(1):i(2),1,2,2));
end
ctsc = [];
cvsc = [];
cssc = [];
chsc = [];
casc = [];
for i = Trial.stc{'s'}.data'
    ctsc(end+1) = mean(circ_dist(ang(i(1):i(2),1,4,1),ang(i(1)-1:i(2)-1,1,4,1)));
    cvsc(end+1) = median(vxy(i(1):i(2),1));
    cssc(end+1) = std(vxy(i(1):i(2),1));
    chsc(end+1) = mean(dsx(i(1):i(2),1,3));
    casc(end+1) = mean(ang(i(1):i(2),1,2,2));
end



figure,
subplot(311),hold on
edg = linspace(-10,0,100);
bar(edg,histc(log10(abs(tsc)),edg),'histc');
h = bar(edg,histc(log10(abs(ctsc)),edg),'histc');h.FaceColor = 'r';h.FaceAlpha= .5;
subplot(312),hold on
edg = linspace(-.5,2,100);
bar(edg,histc(log10(abs(vsc)),edg),'histc');
h = bar(edg,histc(log10(abs(cvsc)),edg),'histc');h.FaceColor = 'r';h.FaceAlpha= .5;
subplot(313),hold on
edg = linspace(-2,2,100);
bar(edg,histc(log10(abs(ssc)),edg),'histc');
h = bar(edg,histc(log10(abs(cssc)),edg),'histc');h.FaceColor = 'r';h.FaceAlpha= .5;

figure,hold on
edg = linspace(0,60,100);
bar(edg,histc(hsc,edg),'histc');
h = bar(edg,histc(chsc,edg),'histc');h.FaceColor = 'r';h.FaceAlpha= .5;


figure,hold on
edg = linspace(0,1.6,100);
bar(edg,histc(asc,edg),'histc');
h = bar(edg,histc(casc,edg),'histc');h.FaceColor = 'r';h.FaceAlpha= .5;



figure,
plot(vxy(:,1)),
Lines(reshape(vprs(log10(abs(vsc))>.75,:),[],1),[],'b');
Lines(Trial.stc{'w'}(:),[],'m');

figure,
hist2([log10(vsc)',tsc'],100,100)




figure,
plot(dsx(:,7,3)),hold on
Lines(reshape(tprs(tsc>.8,1),[],1),[],'g');
Lines(reshape(tprs(tsc>.8,2)-10,[],1),[],'m');
Lines(Trial.stc{'n'}(:),[],'b')


vsc = [];
for i = vprs'
    vsc(end+1,1) = median(dsx(i(1):i(2),7,3));
end


figure,hist2(log10(abs([dsx(nniz(dsx),7,3),circshift(dsx(nniz(dsx),7,3),-40)])),linspace(1.5,2.6,100),linspace(1.5,2.6,100))

figure,hist2(log10(abs([vsc(nniz(vsc)),circshift(vsc(nniz(vsc)),-1)])),linspace(1.5,2.6,100),linspace(1.5,2.6,100))



cvsc = [];
for i =  Trial.stc{'w'}.data'
    cvsc(end+1) = median(vxy(i(1):i(2),2));
end

edg = linspace(-3,2,100);
figure,hold on
bar(edg,histc(log10(vsc),edg),'histc');
h = bar(edg,histc(log10(cvsc),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

edg = linspace(-.5,2,100);
figure,hold on
ind = Trial.stc{'a'};
bar(edg,histc(log10(vxy(ind,1)),edg),'histc');
ind = Trial.stc{'w'};
h = bar(edg,histc(log10(vxy(ind,1)),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

edg = linspace(-3,2,100);
figure,hold on
bar(edg,histc(log10(clip(vxy(Trial.stc{'a-r'},1),.001,200)),edg),'histc');
h = bar(edg,histc(log10(clip(vxy(Trial.stc{'w'}+[-.5,5],1),.001,200)),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

figure,
hist2([log10(vsc)',log10(vsch)'],linspace(-3,2,60),linspace(-3,2,60))

figure, bar(linspace(0,1e-3,100),hist(tsc,linspace(0,1e-3,100)),'histc')
figure, bar(linspace(-9,-3,100),hist(log10(tsc),linspace(-9,-3,100)),'histc')
figure, bar(linspace(-9,-3,100),hist(log10(tsc),linspace(-9,-3,100)),'histc')


dsx = Trial.load('xyz');
rsc = [];
for i = rprs'
    rsc(end+1) = median(dsx(i(1):i(2),7,3));
end

crsc = [];
for i =  Trial.stc{'r'}.data'
    crsc(end+1) = median(dsx(i(1):i(2),7,3));
end

edg = linspace(1.2,2.5,100);
figure,hold on
bar(edg,histc(log10(rsc),edg),'histc');
h = bar(edg,histc(log10(crsc),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;



[ys,fs,ts] = fet_spec(Trial,vl,'mtchglong',false);
ys.data = cat(3,ys(:,:,1,1),ys(:,:,2,2));
ys.resample(30);



ind = Trial.stc{'a'};
mfet = [log10(clip(sum(ys(ind,1:30,1),2),.00001,300)),...
        log10(clip(sum(ys(ind,1:30,2),2),.00001,300))];

figure,
hist2(mfet,...
       linspace(-3,2.5,100),...
       linspace(-3,2.5,100))


   
idx = clusterdata(mfet,'linkage','ward','savememory','on','maxclust',4);

figure
hist2(mfet(idx==2|idx==1,:),...
       linspace(-3,2.5,100),...
       linspace(-3,2.5,100))

vl.resample(ys);

aper = Trial.stc{'a',vl.sampleRate}.cast('TimeSeries');
aper.data = aper.data(1:end-1);
ida = nan(vl.size);
ida(logical(aper.data)) = idx;
figure,plot(log10(clip(vl.data,0.01,200)))
hold on,plot(ida+2)
Lines(Trial.stc{'w',vl.sampleRate}(:),[],'m');               
               

%% Some stuff
ang = create(MTADang,Trial,Trial.load('xyz').filter('ButFilter',3,10,'low'));

wfet = Trial.xyz.copy;
wfet.data= [0;diff(abs(circ_dist(ang(:,'pelvis_root','head_back',1),ang(:,'spine_lower','spine_middle',1))))];

defspec = struct('nFFT',2^8,'Fs',wfet.sampleRate,...
                           'WinLength',2^7,'nOverlap',2^7*.875,...
                           'FreqRange',[1,15]);
[wys,fs,ts] = fet_spec(Trial,wfet,'mtchglong',false);
wys.resample(30);
wys.data(wys.data<=0) =1e-10;
figure,imagesc(ts,fs,log10(wys.data(nniz(wys(:,10)),2:end)')),axis xy, 
caxis([-10,-4])
colormap jet

ind = Trial.stc{'w'};
mfet = [log10(clip(sum(ys(ind,1:20,2),2),.00001,300)),...
        log10(clip(sum(wys(ind,1:20),2),1e-10,300))];

figure,
hist2(mfet,...
       linspace(-3,2.5,100),...
       linspace(-10,-3.5,100))

%% Some more stuff
vl.resample(ys);
wper = Trial.stc{'s',ys.sampleRate};
wfet = [];
for i = wper.data',
    wfet(end+1) = mean(log10(vl(i(1):i(2),1)));
end



%% oh god no

fet = fet_lgr(Trial);
fet.resample(30);
fwin = 120;
hwin = fwin/2;
tn = 10000;

ufet = fet.copy;
ufet.data = nunity(fet.data(:,[1:19]));
ufet.filter('ButFilter',3,1,'low');
figure,plot(mean(abs(diff(ufet.data)),2));
Lines(Trial.stc{'r',fet.sampleRate}(:),[],'r');
Lines(Trial.stc{'w',fet.sampleRate}(:),[],'b');


ufet.data = fet(1:tn,10:19);

msfet = fet.copy;
msfet.data = fet(1:tn,10:19);
msfet.data = reshape(msfet.segs(1:msfet.size(1),fwin,nan),[hwin,2,msfet.size]);
msfet.data = msfet(:,:,1:tn,:);
mfet = mean(msfet.data);

% COV MAT for all points
cfet = reshape(repmat(permute(bsxfun(@minus,msfet.data,mfet),[4,1,2,3]),[ufet.size(2),1,1,1]),[ufet.size([2,2]),hwin,2,tn]);
cmat = sq(sum(cfet.*permute(cfet,[2,1,3,4,5]),3));

% MEAN VEC for all points
mfet = permute(mfet,[4,1,2,3]);


%trind = (1:19).*uufet.size(2)-fliplr(0:uufet.size(2)-1);

kld = zeros([tn,1]);
for i = 1:tn,
kld(i) = .5*(trace(cmat(:,:,1,i)\cmat(:,:,2,i))...
          +diff(mfet(:,:,:,i),1,3)'*(1\cmat(:,:,1,i))*diff(mfet(:,:,:,i),1,3)...
          -ufet.size(2)...
          +log2(det(cmat(:,:,1,i))/det(cmat(:,:,2,i))));
end

figure
plot((1:tn)-15,log10(abs(kld)))
Lines(Trial.stc{'r',fet.sampleRate}(:),[],'r');
Lines(Trial.stc{'w',fet.sampleRate}(:),[],'b');


stc='wrnsmk';
wm = {};
for s = 1:6
wper = Trial.stc{stc(s),ufet.sampleRate};
wm{s} = [];
for i = wper.data',
    wm{s}(end+1,:) = mean(ufet(i(1):i(2),:));
end
end

figure,plot(plot(diff(fet(:,1))))



% Effects of threshold on period count

Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120317');
Trial = MTATrial('Ed05-20140529','all','ont');
Trial = MTATrial('Ed01-20140707');


Trials = SessionList('mypc_test_grp');
nstuff = {};
tc = 1;
for t = Trials,
    
    Trial = MTATrial(t.name,t.trialName,t.mazeName);
    
    
    dsx = Trial.load('xyz');
    
    vxy = dsx.vel([1,7],[1,2]);
    vxy.filter('ButFilter',3,10,'low');
    
    % Walking Feature
    vl = dsx.vel([1,2,7],[1,2]);
    %vl.filter('ButFilter',3,.75,'low');
    vl.data = cat(1,diff(vxy.data),zeros([1,vxy.size(2)]));
    vl.filter('ButFilter',3,.47,'low');
    
    % Rearing Feature
    %vz = dsx.vel(7,3);
    vz = dsx.copy;
    vz.data = [0;diff(dsx(:,7,3))];
    vz.filter('ButFilter',3,.6,'low');
    
    % Turning Feature
    dsx.filter('ButFilter',3,.5,'low');
    ang = create(MTADang,Trial,dsx);
    tfet = abs(diff(circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1))));
    
    fet = {};
    fet{1} = vz(:,1);
    fet{2} = vl(:,1);
    fet{3} = tfet;
    
    thr_vls = -2:.1:10;
    
    % nstuff = {Trial,fet}
    for f = 1:numel(fet),
        nsmp_vl = [];
        n = [];
        ff = [];
        for thr_vl = thr_vls;
            vpnts = LocalMinima(-log10(abs(fet{f})),round(.20*dsx.sampleRate),thr_vl);            
            dvpnts = diff(vpnts);
            
            
            ff = cat(1,ff,mean(dvpnts)/std(dvpnts));
            if numel(vpnts)>10,
                n = cat(2,n,histc(log10(dvpnts),linspace(1,6,100)));
            else
                n = cat(2,n,zeros([100,1]));
            end
            nsmp_vl(end+1) = numel(vpnts)-1;
        end
        nstuff{tc,f} = skewness(n);
    end
    tc = tc+1;
end

figure,plot(thr_vls,nsmp_vl)
figure,plot(thr_vls(1:end-1),diff(nsmp_vl))

polyfit(thr_vls(10:30),nsmp_vl(10:30),1)

%jg05-20120310 pfit: 2485.74025974026         -1691.14718614719


nn = cell2mat(nstuff);
nn = reshape(nn,9,121,3);
figure,plot(thr_vls,mean(nn(:,:,2))')

thr_vl = 2.3;
f = 1;

vpnts = LocalMinima(-log10(abs(fet{f})),round(.20*dsx.sampleRate),thr_vl);
vprs = [vpnts,circshift(vpnts,-1)];
vprs(end,:) = [];
figure, plot(log10(abs(fet))),Lines(vprs(:),[],'c');Lines(Trial.stc{'r'}(:),[],'r');

ofet = [xyz(:,7,3),...
       vxy(:,1),...
       circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-30))];

%% Thresholds found by peak in the skewness of iei of distributions over increasing feature thresholds
thr_fet = [2.3,1.8,7.1];

xyz = Trial.load('xyz');
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,.6,'low');
ang = create(MTADang,Trial,fxyz);


%wf = sqrt(sum(sum(circshift(xyz(:,1:4,[1,2]),-6)-xyz(:,1:4,[1,2]),2).^2,3));
%wf = MTADxyz('data',wf,'sampleRate',xyz.sampleRate);
%wf.filter('ButFilter',3,1,'low');

vsc = zeros([xyz.size(1),5]);
for f = 1:numel(fet),
    vpnts = LocalMinima(-log10(abs(fet{f})),round(.20*dsx.sampleRate),thr_fet(f));
    vprs = [vpnts,circshift(vpnts,-1)];
    vprs(end,:) = [];
    
    for i = vprs'
        %    vsc(i(1):i(2),f) = median(ofet(i(1):i(2),f));
        if f == 1,
            vsc(i(1):i(2),f) = median(ofet(i(1):i(2),f));
        elseif f==2,
            %vsc(i(1):i(2),f) = median(ofet(i(1):i(2),f));
            %vsc(i(1):i(2),f) = median(ofet(i(1):i(2),f));
            vsc(i(1):i(2),f) = abs(sqrt(sum(diff(dsx(i,1,[1,2])).^2,3)));
            vang = sq(diff(dsx(i,1,[1,2])));
            vang = atan2(vang(2),vang(1));
            vsc(i(1):i(2),4) = cos(circ_dist(circ_mean(ang(i(1):i(2),1,4,1)),vang));    
            
        elseif f==3,
            vsc(i(1):i(2),f) = abs(circ_dist(ang(i(1),1,4,1),ang(i(2),1,4,1)))/diff(i);
            vang = sq(diff(dsx(i,4,[1,2])));
            vang = atan2(vang(2),vang(1));
            vsc(i(1):i(2),5) = sin(circ_dist(circ_mean(ang(i(1):i(2),1,4,1)),vang));            

        end
        
        
    end
end

mvsc = MTADxyz('data',vsc,'sampleRate',xyz.sampleRate);


%% Confirm Quality of segmentation visually
edgs = linspace(-5,3,100);
figure,hold on
    bar(edgs,histc(unique(log10(abs(mvsc(Trial.stc{'a-w-r-n'},4)))),edgs),'histc');
h = bar(edgs,histc(unique(log10(abs(mvsc(Trial.stc{'w'},4)))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

edgs = linspace(-5,3,100);
figure,hold on
bar(edgs,histc((log10(abs(mvsc(Trial.stc{'a'},2)))),edgs),'histc');
h = bar(edgs,histc((log10(abs(mvsc(Trial.stc{'w'},2)))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


edgs = linspace(-1,2,100);
figure,hold on
bar(edgs,histc(log10(abs(mvsc(Trial.stc{'a-w-n-r-m'},2).*mvsc(Trial.stc{'a-w-n-r-m'},4))),edgs),'histc');
h = bar(edgs,histc(log10(abs(mvsc(Trial.stc{'w'},2).*mvsc(Trial.stc{'w'},4))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;
h = bar(edgs,histc(log10(abs(mvsc(Trial.stc{'n'},2).*mvsc(Trial.stc{'n'},4))),edgs),'histc');
h.FaceColor = 'g';
h.FaceAlpha = .5;

figure,plot(abs(vsc(:,4))),
hold on,plot(log10(vsc(:,2)),'r')
Lines(Trial.stc{'w'}(:),[],'m');
Lines(Trial.stc{'n'}(:),[],'g');

figure,hold on
edgs = linspace(-8,-1,100);
edgs = linspace(-3,2,100);
aind = Trial.stc{'a-w-r'};
sind = Trial.stc{'w'};
afet = log10(abs(mvsc(aind,4)).*abs(mvsc(aind,2)));
sfet = log10(abs(mvsc(sind,4)).*abs(mvsc(sind,2)));
bar(edgs,histc(afet,edgs),'histc');
h = bar(edgs,histc(sfet,edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


%% Recreate for continous situation for speed X coherent movement of LS marker with direction of body movement

% .5 sec lag
lag = 20;
vang = {circshift(dsx(:,1,1),lag)-circshift(dsx(:,1,1),-lag),circshift(dsx(:,1,2),lag)-circshift(dsx(:,1,2),-lag)};
vang = atan2(vang{:});
vang = abs(cos(circ_dist(circ_mean(GetSegs(ang(:,1,3,1),1:ang.size(1),30,nan))',vang)));    

bang = {circshift(dsx(:,3,1),lag)-circshift(dsx(:,3,1),-lag),circshift(dsx(:,3,2),lag)-circshift(dsx(:,3,2),-lag)};
bang = atan2(bang{:});
bang = abs(cos(circ_dist(circ_mean(GetSegs(ang(:,1,3,1),1:ang.size(1),30,nan))',bang)));    




% msx = Trial.load('xyz');
% mxy = msx.vel([],[1,2]);
% mxy.filter('ButFilter',3,2,'low');

figure,hold on
plot(mxy(:,1))
plot(vang*20,'r')
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');            

figure,plot(median(mxy(:,1:3),2).*vang.*bang)
Lines(Trial.stc{'w'}(:),[],'b');
Lines(Trial.stc{'n'}(:),[],'g');            

vang = MTADxyz('data',vang,'sampleRate',dsx.sampleRate);
bang = MTADxyz('data',bang,'sampleRate',dsx.sampleRate);

figure,hold on
edgs = linspace(-3,2.5,100);
aind = Trial.stc{'a-w-r'};
sind = Trial.stc{'w'};
afet = log10(abs(median(mxy(aind,1:3),2).*vang(aind).*bang(aind)));
sfet = log10(abs(median(mxy(sind,1:3),2).*vang(sind).*bang(sind)));
bar(edgs,histc(afet,edgs),'histc');
h = bar(edgs,histc(sfet,edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;



wn = 61;
wh = 30;
t = [-wh:wh]';

x = 60;

sg = dsx(x+1:x+wn,7,3);
[pn,st] = polyfit(t,sg,1);
[pn,st2] = polyfit(t,sg,2);

[pn,st3] = polyfit(t,sg,3);


% IDEA if 

z = dsx.copy;
z.data = dsx(:,7,3);
zs = z.segs(1:z.size(1),350,nan);

zv = circshift(var(zs)',175);
figure,hist(log10(zv),100);

figure,
plot(log10(zv))
Lines(Trial.stc{'r'}(:),[],'r');
Lines([],2.5,'k');



z = vxy.copy;
z.data = vxy(1:4:end,1);
zs = z.segs(1:z.size(1),120,nan);

zv = MTADxyz('data',circshift(var(zs)',60),'sampleRate',vxy.sampleRate/4);
figure,hist(log10(zv.data),100);

zv.data(~nniz(zv)) = 0;
zv.filter('ButFilter',3,.1,'low');

figure,
plot(log10(zv.data))
Lines(Trial.stc{'w',zv.sampleRate}(:),[],'r');
Lines([],1,'k');



edgs = linspace(-1.5,3,100);
figure,hold on,
ind = Trial.stc{'a-w-n-r'};
bar(edgs,histc(log10(zv(ind)),edgs),'histc')
ind = JoinRanges(Trial.stc{'w',zv.sampleRate}.data,Trial.stc{'n',zv.sampleRate}.data);
ind = JoinRanges(ind,Trial.stc{'r',zv.sampleRate}.data);
h = bar(edgs,histc(log10(zv(ind)),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


figure,
ind = Trial.stc{'a',zv.sampleRate};
edgs = linspace(-1.5,3,10);
hist2(log10([zv(ind),circshift(zv(ind),25)]),edgs,edgs)




%lr fet god noooooooooooo

Trial = MTATrial('jg05-20120317');
sampleRate = 30;
fet = fet_lgr(Trial);
fet.resample(sampleRate);
aper = Trial.stc{'a'};


states = {'sit','rear','groom','turn','walk'};
B = {};
new_smat = [];
for s = states
% compute lgr coefficients
[smat] = max(stc2mat(Trial.stc,fet,{s}),[],2);
ind = resample(cast(aper.copy,'TimeSeries'),fet);
ind = logical(ind.data);
[B{end+1}] = mnrfit(fet(ind,:),smat(ind)+1,'model','nominal');

% compute lgr score of original feature
d_state = mnrval(B{end},fet.data);
[~,dind] = max(d_state,[],2);
new_smat(:,end+1) = dind-1;

% remove periods from aper which were scored as target behavior
nper = ThreshCross(dind-1,.5,10);
aper = aper-nper;
end


% hierarchy lies in the order of states
tnew_smat = new_smat;
for i=fliplr(1:size(new_smat,2)-1),
    tnew_smat(:,i+1) = new_smat(:,i+1).*double(~sum(new_smat(:,1:i),2));
end
Stc = Trial.stc.copy;
Stc.updateMode('lgr_labeled');
Stc.states = {};
nStc = mat2stc(tnew_smat,Stc,sampleRate,states,{'s','r','m','n','w'});

smat = ~~stc2mat(Trial.stc,fet,states);

cmat = zeros([numel(states),numel(states)]);
for i = 1:numel(states),
    for j = 1:numel(states),
        cmat(i,j) = sum(double(smat(:,i)==1&tnew_smat(:,j)==1));
    end
end

sacc = bsxfun(@rdivide,diag(cmat),sum(cmat,2));

Data.fillgaps(.2);


Trial = MTATrial('jg05-20120310');
sampleRate = 30;
fet = fet_lgr(Trial);
fet.resample(sampleRate);
aper = Trial.stc{'a'};


test_smat = [];
for i = 1:numel(states),
    d_state = mnrval(B{i},fet.data);
    [~,dind] = max(d_state,[],2);
    test_smat(:,end+1) = dind-1;
end
for i=fliplr(1:size(test_smat,2)-1),
    test_smat(:,i+1) = test_smat(:,i+1).*double(~sum(test_smat(:,1:i),2));
end

figure,imagesc(test_smat')










