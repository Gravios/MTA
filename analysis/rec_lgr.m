


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
vl = dsx.acc([1,2,7],[1,2]);
vl.filter('ButFilter',3,.5,'low');
vz = dsx.vel(4,3);
vz.filter('ButFilter',3,.5,'low');
dsx.filter('ButFilter',3,.5,'low');
ang = create(MTADang,Trial,dsx);
tfet = abs(diff(circ_dist(ang(:,1,4,1),circshift(ang(:,1,4,1),-1))));
figure,hold on,plot(tfet,'m')
Lines(Trial.stc{'n'}(:),[],'g')
Lines(Trial.stc{'r'}(:),[],'r')
Lines(Trial.stc{'w'}(:),[],'b')

tpnts = LocalMinima(-tfet,round(.25*dsx.sampleRate),0);
tprs = [tpnts,circshift(tpnts,-1)];
tprs(end,:) = [];
tsc = [];
for i = tprs'
    tsc(end+1) = abs(circ_dist(ang(i(1),1,4,1),ang(i(2),1,4,1)))/diff(i);
end
ctsc = [];
for i = Trial.stc{'n'}.data'
    ctsc(end+1) = abs(circ_dist(ang(i(1),1,4,1),ang(i(2),1,4,1)))/diff(i);
end
edg = linspace(-5,-1,100);
figure,hold on
bar(edg,histc(log10(tsc),edg),'histc');
h = bar(edg,histc(log10(ctsc),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


tsc = [];
for i = tprs'
    tsc(end+1) = mean(circ_dist(ang(i(1):i(2),1,4,1),ang(i(1)-1:i(2)-1,1,4,1)))/diff(i);
end
ctsc = [];
for i = Trial.stc{'n'}.data'
    ctsc(end+1) = mean(circ_dist(ang(i(1):i(2),1,4,1),ang(i(1)-1:i(2)-1,1,4,1)))/diff(i);
end

edg = linspace(-10,-2.5,100);
figure,hold on
bar(edg,histc(log10(abs(tsc)),edg),'histc');
h = bar(edg,histc(log10(abs(ctsc)),edg),'histc');
h.FaceColor = 'r';
h.FaceAlpha= .5;


dsx = Trial.load('xyz');
vxy = dsx.vel([1,7],[1,2]);
vxy.filter('ButFilter',3,2,'low');

figure,
hist2([log10(vsc)',tsc'],100,100)




figure,
plot(dsx(:,7,3)),hold on
Lines(reshape(tprs(tsc>.8,1),[],1),[],'g');
Lines(reshape(tprs(tsc>.8,2)-10,[],1),[],'m');
Lines(Trial.stc{'n'}(:),[],'b')

vpnts = LocalMinima(-abs(vl(:,1)),round(.25*dsx.sampleRate),0);
vprs = [vpnts,circshift(vpnts,-1)];
vprs(end,:) = [];
vsc = [];
for i = vprs'
    vsc(end+1) = median(vxy(i(1):i(2),1));
end

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
rpnts = LocalMinima(-abs(vz(:,1)),round(.25*dsx.sampleRate),0);
rprs = [rpnts,circshift(rpnts,-1)];
rprs(end,:) = [];
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

