% req20231006
%     Tags: ccg theta power estimation lfp state
%     Status: Active
%     Type: Utility
%     Author: Justin Graboski
%     Final_Forms: N/A
%     Project: General
%     Description: lfp state segmentation

% it works ... for a

xyz = preproc_xyz(Trial,'trb',120);
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
lvxy = copy(vxy);
lvxy.data(lvxy.data<=0.0001) = 0.0001;
lvxy.data = log10(lvxy.data);

stcm = stc2mat(Trial.stc,afet,{'theta','rear','loc','pause','sit','groom'});

lfpPyrC = Trial.load('lfp',[52]);
lfpPyr = Trial.load('lfp',[49,55]);
lfpPyrD = copy(lfpPyrC);
lfpPyrD.data = diff(lfpPyr.data,1,2);
lfpPyr = Trial.load('lfp',[70,73]);
lfpPyrW = Trial.load('lfp',[68,75]);
lfpRad = Trial.load('lfp',[77,78]);
lfpLMo = Trial.load('lfp',[82,83]);

figure();
subplot(211);
hold('on');
plot(diff(lfpPyr.data,1,2))
plot(diff(lfpRad.data,1,2))
plot(diff(lfpLMo.data,1,2))
subplot(212);
hold('on');
plot(lfpPyr.data(:,1))
plot(lfpRad.data(:,1))
plot(lfpLMo.data(:,1))
linkax('x')




figure,
subplot(212)
plot(

flfpRad = copy(lfpRad)
flfpRad.data = flfpRad.data(:,2);
flfpRad.filter('ButFilter',4,[35,55],'bandpass');

flfpLMo = copy(lfpLMo)
flfpLMo.data = flfpLMo.data(:,2);
flfpLMo.filter('ButFilter',4,[75,120],'bandpass');

hamm = hamming(64);

gpowRad = sum(bsxfun(@times, flfpRad.segs(1:size(flfpRad,1),64),hamm).^2)';
gpowLMo = sum(bsxfun(@times, flfpLMo.segs(1:size(flfpLMo,1),64),hamm).^2)';



phz = load_theta_phase(Trial,lfp.sampleRate);
thetaTroughs = abs(phz.data-pi)<0.1&[0;diff(phz.data-pi)]>0;
thetaTroughs = LocalMinima(-convn(thetaTroughs,ones([1,21]),'same'),0,21);



figure();
subplot(5,1,[1]);
    imagesc(ts,fs,log10(ggys.data)');    axis(gca(),'xy');    colormap(gca,'jet');
subplot(5,1,[2]);
    imagesc(ts,fs,log10(gpys.data)');    axis(gca(),'xy');    colormap(gca,'jet');
subplot(5,1,[3,4]);
    hold('on');
    plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpPyr(:,1)+10000,'k')
    plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, diff(lfpPyr.data,1,2),'k')
    plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpLMo.data(:,1)-15000)
subplot(515);
    plotSTC(Trial.stc,1)
linkax('x')

figure();
subplot(411);
hold('on');
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpPyr(:,1),'k')
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, diff(lfpPyr.data,1,2),'k')
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpLMo.data(:,1))
subplot(412);
hold('on');
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpPyr.data(:,1)*3,'k')
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpRad.data(:,1))
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpLMo.data(:,1))
subplot(413);
hold('on');
%plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, nunity(diff(lfpPyr.data,1,2))-nunity(lfpLMo.data(:,1)),'k')
plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, nunity(lfpLMo.data(:,2))-dlfpPyrFilt.data,'m')
subplot(414);
plotSTC(Trial.stc,1)
linkax('x')



lfpRC = lfpPyr.copy();  lfpRC.data = diff(lfpRC.data,1,2);
lfpRW = lfpPyrW.copy(); lfpRW.data = diff(lfpRW.data,1,2);


doCrossSpec = true;

defspec = struct('nFFT',2^12,'Fs',lfpRC.sampleRate,...
                 'WinLength',2^11,'nOverlap',2^11*0.875,...
                 'FreqRange',[1,20]);
[ysCS,fs,ts] = fet_spec(Trial,lfpPyr, 'mtcsdglong',false,[],defspec,[],doCrossSpec);
[ysRC,fs,ts] = fet_spec(Trial,lfpPyrD,          [],false,[],defspec);
[ysPC,fs,ts] = fet_spec(Trial,lfpPyrC,          [],true,[],defspec);

fpyrPhz = copy(ysCS);
%fpyrPhz.data = imgaussfilt(angle(ysCS(:,:,1,2)),0.8);
fpyrPhz.data = angle(ysCS(:,:,1,2));

fysRC = copy(ysRC);
fysRC.data = imgaussfilt(log10(ysRC.data),3);
fysRC.filter('ButFilter',4,0.05,'low');
fysPC = copy(ysPC);
fysPC.data = imgaussfilt(log10(ysPC.data),3);
fysPC.filter('ButFilter',4,0.05,'low');
fysO = copy(ysCS);
fysO.data = imgaussfilt(log10(ysCS.data(:,:,1,1)),3);
fysO.filter('ButFilter',4,0.05,'low');
figure,
%imagesc(ts,fs,log10(ysO.data)');
subplot(211);
%imagesc(ts,fs,angle(ysCS(:,:,1,2))');
imagesc(ts,fs,fysPC(:,:,1,1)');
axis(gca(),'xy');
colormap(gca(),'jet');
subplot(212);
imagesc(ts,fs,log10(ysCS(:,:,1,1))');
axis(gca(),'xy');
colormap(gca(),'jet');
linkax('xy')



pyrPhz=copy(ysCS); pyrPhz.data = mean(angle(ysCS(:,6<fs&fs<10,1,2)),2);
lpow = copy(ysPC); lpow.data = mean(fysPC(:,fs<18),2);
dpow = copy(ysPC); dpow.data = mean(fysPC(:,fs<4),2);
rcr =  copy(ysRC); rcr.data  = mean(fysO(:,fs<9&fs>6),2)./mean(fysRC(:,fs<9&fs>6),2);
rpow = copy(ysRC); rpow.data = mean(fysRC(:,fs<9&fs>6),2);
tpow = copy(ysPC); tpow.data = mean(fysPC(:,fs<9&fs>6),2)-mean(fysPC(:,fs<5|(fs<14&fs>11)),2);
opow = copy(ysCS); opow.data = mean(fysO(:,fs<9&fs>6),2);
spow = copy(ysPC); spow.data = sum(abs(bsxfun(@minus,fysO.data,mean(fysO.data,2))),2);



ind = [Trial.stc{'walk+turn+pause+rear'}];
ind.fillgaps()
ind.resample(spow);
figure,plot(spow(:),tpow(:),'.')
hold('on');,plot(spow(ind),tpow(ind),'.r')



figure,
subplot(511);
    imagesc(ggts,ggfs,log10(ggys.data)');
    axis(gca(),'xy');caxis(gca(),[-1,3]);colormap(gca,'jet')
subplot(512);
    imagesc(wts,wfs,log10(wys.data)');
    axis('xy');colormap(gca,'jet')
subplot(513);
    hold(gca(),'on');
    plot(ggts,nunity(sum(log10(ggys.data)+6,2)),'c');
    plot(ggts,nunity(sum(log10(gpys(:,ggfs<110))-4,2)),'b');
    plot(ggts,nunity(sum(log10(ggys(:,ggfs<10&ggfs>6))+5.5,2)./sum(log10(wys(:,fs<10&fs>6))-1,2)),'r')
    plot(ggts,log10(mean(ggys(:,ggfs<10&ggfs>6),2))-log10(mean(ggys(:,ggfs<5|(ggfs<20&ggfs>10)),2)));
    Lines([],0,'k');
    Lines([],7,'k');
subplot(514);
    plot([1:size(phzRW,1)]./phzRW.sampleRate,circ_dist(phzRW(:),phzRC(:)))
    Lines([],0,'k');    
subplot(515);
    plotSTC(Trial.stc,1)
linkax('x')


bhvStates = {'sit-theta','sit&theta','walk+turn+pause+rear&theta',...
             'groom&theta'};
figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    histogram(circ_dist(phzRW(ind),phzRC(ind)),linspace(-1,1,100))    
    Lines(mean(circ_dist(phzRW(ind),phzRC(ind))),[],'r');
end
%ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
%ForAllSubplots('colormap(gca(),''jet'')');



figure,
hold('on')
ind = [Trial.stc{'theta'}];
plot(wys(ind,12),ggys(ind,12),'.b')
ind = [Trial.stc{'walk&theta'}];
plot(wys(ind,12),ggys(ind,12),'.r')
ind = [Trial.stc{'sit&theta'}];
plot(wys(ind,12),ggys(ind,12),'.c')

figure
subplot(511);
ind = [Trial.stc{'walk&theta'}];
histogram(circ_dist(phzP(ind,1),phzP(ind,2)),linspace(-pi,pi,100))
subplot(512);
ind = [Trial.stc{'pause&theta'}];
histogram(circ_dist(phzP(ind,1),phzP(ind,2)),linspace(-pi,pi,100))
subplot(513);
ind = [Trial.stc{'rear&theta'}];
histogram(circ_dist(phzP(ind,1),phzP(ind,2)),linspace(-pi,pi,100))
subplot(514);
ind = [Trial.stc{'groom&theta'}];
histogram(circ_dist(phzP(ind,1),phzP(ind,2)),linspace(-pi,pi,100))
subplot(515);
ind = [Trial.stc{'sit&theta'}];
histogram(circ_dist(phzP(ind,1),phzP(ind,2)),linspace(-pi,pi,100))

phzP = lfpPyr.copy()
phzP.resample(120);
phzP = phzP.phase([5,12],4);
phzPD = phzP.copy();
phzPD.data = circ_dist(phzPD(:,1),phzPD(:,2));

tpow = ggys.copy();
tpow.data = log10(mean(gpys(:,ggfs<10&ggfs>6),2));%-log10(mean(gpys(:,ggfs<5|(ggfs<20&ggfs>10)),2));
tpow.data(~nniz(tpow.data)) = 1;
tpow.resample(phzPD);

rpow = ggys.copy();
rpow.data = log10(mean(ggys(:,ggfs<10&ggfs>6),2));%-log10(mean(ggys(:,ggfs<5|(ggfs<20&ggfs>10)),2));
rpow.data(~nniz(rpow.data)) = 1;
rpow.resample(phzRC);

rcr = ysCS.copy();
rcr.data = nunity(log10(mean(gpys(:,ggfs<10&ggfs>5),2))-log10(mean(ggys(:,ggfs<10&ggfs>5),2)));
rcr.data(~nniz(rcr.data)) = 1;
rcr.resample(phzPD);


figure,
ind = [Trial.stc{'sit&theta'}];
plot(rcr(ind),phzPD(ind),'.')
hold(gca(),'on');
ind = [Trial.stc{'walk&theta'}];
plot(rcr(ind),phzPD(ind),'.r')


figure,
subplot(411);
ind = [Trial.stc{'sit&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
Lines(-0.5,[],'w');Lines([],-0.5,'w');
subplot(412);
ind = [Trial.stc{'walk&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
Lines(-0.5,[],'w');Lines([],-0.5,'w');
subplot(413);
ind = [Trial.stc{'pause&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
Lines(-0.5,[],'w');Lines([],-0.5,'w');
subplot(414);
ind = [Trial.stc{'groom&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
Lines(-0.5,[],'w');Lines([],-0.5,'w');
colormap('jet');


bhvStates = {'sit-theta','sit&theta','walk+turn+pause+rear&theta',...
             'groom&theta'};
             
figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rcr(ind),phzPD(ind)],linspace(-3,3.5,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([0,1000])');

figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rcr(ind),pyrPhz(ind)],linspace(-3,3.5,50),linspace(-2,1,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([0,1000])');





figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([circ_dist(phzRW(ind),phzRC(ind)),phzPD(ind)],linspace(-0.5,0.5,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.05,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([100,1000])');

figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');



figure,
p = 1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);    
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rpow(ind),pyrPhz(ind)],linspace(-0.25,0.8,50),linspace(-2,1,50))
    Lines(0.25,[],'w');Lines([],-0.75,'w')
    if s==1,title('rpow vs pyrPhz');end
end
p = p+1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rpow(ind),lpow(ind)],linspace(-0.25,0.8,50),linspace(4.75,6,50))
    Lines(0.25,[],'w');
    Lines([],5.6,'w');
    if s==1,title('rpow vs lpow');end            
end
p = p+1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rpow(ind),tpow(ind)],linspace(-0.25,0.8,50),linspace(-2,3,50))
    Lines(0.25,[],'w');
    Lines([],0.2,'w');
    if s==1,title('rpow vs tpow');end        
end
p = p+1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),pyrPhz(ind)],linspace(-2,3,50),linspace(-2,1,50))
    Lines(0.25,[],'w');
    Lines([],-.6,'w');
    if s==1,title('tpow vs pyrPhz');end    
end
p = p+1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),rcr(ind)],linspace(-2,3,50),linspace(-pi,pi,50))
    Lines(-0.25,[],'w');
    Lines([],-.6,'w');
    if s==1,title('tpow vs rcr');end    
end
p = p+1;
for s = 1:numel(bhvStates)
    subplot2(numel(bhvStates),6,s,p);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rpow(ind),rcr(ind)],linspace(-0.25,0.8,50),linspace(-pi,pi,50))
    Lines(-0.25,[],'w');
    Lines([],-.6,'w');
    if s==1,title('rpow vs rcr');end    
end
ForAllSubplots('colormap(gca(),''jet'')');
%ForAllSubplots('caxis([10,800])');



afet = copy(rcr);    
afet.data = nunity([tpow(:,1), rpow(:,1), lpow(:,1), pyrPhz(:,1), rcr(:,1), dpow(:,1),opow(:,1)]); %rcr(:,1),dpow(:,1)
afet.data(~nniz(afet.data),:) = 0;
afetMeanW = median(afet(Trial.stc{'w'},:),'omitnan');
%afet.data = bsxfun(@minus,afet.data,afetMeanW);
afetLims = prctile(afet.data,[0.1,99.9])'


clear('hmm');
updateOM = 1;
hmm.K = 7;
hmm = hmminit(afet.data(:,2:end),hmm,'full');
hmm.train.cyc = 200;
hmm.obsmodel='Gauss';
hmm.train.obsupdate=ones([1,hmm.K])*updateOM;
hmm.train.init = 1;
hmm = hmmtrain(afet.data(:,2:end),size(afet,1),hmm);

hmmStateMatrix = zeros([size(afet,1),7]);

[decode] = hmmdecode(afet.data(:,2:end),size(afet,1),hmm);
hmmStateMatrix(:,1) = decode.q_star';



model_path = create_directory(fullfile(Trial.path.project,'analysis', 'hmm_models'));
for modelInd = 2:7
    clear('hm');
    hm.K = 7;
    hm = hmminit(afet.data(:,[1:7]~=modelInd),hm,'full');
    hm.train.cyc = 200;
    hm.obsmodel='Gauss';
    hm.train.obsupdate=ones([1,hm.K])*updateOM;
    hm.train.init = 1;
    hmm(modelInd) = hmmtrain(afet.data(:,[1:7]~=modelInd),size(afet,1),hm);
    
    % trcpow hrcpow tpow


    %load( fullfile( model_path ,['test_hmm_model.mat']),'hmm');

    diag(hmm(modelInd).P)

    % COMPUTE hmm states
    [decode] = hmmdecode(afet.data(:,[1:7]~=modelInd),size(afet,1),hmm(modelInd));
    hmmStateMatrix(:,modelInd) = decode.q_star';
end
save( fullfile( model_path ,['test_hmm_models.mat']),'hmm');

modelInd = 1;
[decode] = hmmdecode( afet.data(:,[1:7]~=modelInd), size(afet,1), hmm(modelInd));
hmmStateMatrix(:,modelInd) = decode.q_star';

for modelInd = 1:7
    [decode] = hmmdecode(afet.data(:,[1:7]~=modelInd),size(afet,1),hmm(modelInd));
    hmmStateMatrix(:,modelInd) = decode.q_star';
end

cm =cov(afet.data);
mu = zeros([7,7,7]);
for a = 1:7; for p = 1:7
mu(p,:,a) = mean(afet(hmmStateMatrix(nniz(hmmStateMatrix(:,a)),a)==p,:),'omitnan');
end;end

for modelInd = 2:7
    mmat = mfun(hmmStateMatrix(:,1),hmmStateMatrix(:,modelInd));
    [~,mInd] = max(mmat,[],2);
    degMat = (bsxfun(@minus,mInd',mInd)+eye(7))==0;
    assert(numel(nonzeros(degMat))/2==1,'prepare for headache');
    [i,j] = ind2sub([7,7],find(degMat,1,'first'))
    smat = bsxfun(@eq,hmmStateMatrix(:,modelInd),[1:7]);
    smat = smat(:,mInd);
    [~,hmmStateMatrix(:,modelInd)] = max(smat,[],2);
end


hmmStateMatrix(:,8:end) = [];
sortedHmmStateMatrix = hmmStateMatrix;

for modelInd = 2:7
    for p = 1:7,
        for d = 1:7,
            vd(d,p)  = sqrt((mu(p,:,1)-mu(d,:,modelInd))*cm*(mu(p,:,1)-mu(d,:,modelInd))');
        end
    end
    [~,sind] = sort(vd(:));
    [i,j] = ind2sub([7,7],sind);
    dmapInd = [i,j];
    source = 1:7;
    sink =   1:7;
    map = [];
    k =1;
    while ~isempty(source)&&~isempty(sink)
        maprow = dmapInd(1,:);
        if ismember(maprow(1),sink) && ismember(maprow(2),source)
            map(k,:) = maprow;
            k = k+1;
            sink(sink==maprow(1)) = [];
            source(source==maprow(2)) = [];
        end
        dmapInd(1,:) = [];    
    end

    % generate new hmmStateMatrix
    for s = 1:7
        sortedHmmStateMatrix(hmmStateMatrix(:,modelInd)==map(s,1),modelInd) = map(s,2);
    end
end

newHmmStateMatrix = sortedHmmStateMatrix;

mixedStates =  unique(newHmmStateMatrix,'rows','sorted');
msCound = zeros([size(mixedStates,1),1]);
for i = 1:size(mixedStates,1)
    ind = ismember(newHmmStateMatrix,mixedStates(i,:),'rows');
    msCount(i) = sum(ind);
end
hcMixedStates = mixedStates(msCount>100,:);
dd = zeros(size(hcMixedStates,1)*[1,1]);
for i = 1:size(hcMixedStates,1)
    ind = ismember(newHmmStateMatrix,hcMixedStates(i,:),'rows');
    if sum(ind)>1    
        muA = mean(afet(ind,:));
    else
        muA = afet(ind,:);
    end
    
    coA = inv(cov(afet(ind,:)));
    for j = i+1:size(hcMixedStates,1)
        indB = ismember(newHmmStateMatrix,hcMixedStates(j,:),'rows');
        if sum(indB)>1    
            muB = mean(afet(indB,:));
        else
            muB = afet(indB,:);
        end
        muB = mean(afet(indB,:));
        dd(i,j) = sqrt((muB-muA)*coA*(muB-muA)');
    end
end
% dd fill in dd triu
ddt = squareform(dd + rot90(fliplr(dd),1));
Z = linkage(ddt,'complete');
T = cluster(Z,'maxclust',7);
[~,cind] = sort(T);    
currentClu = 1;
stateID = hcMixedStates(cind(1),:);
for t = 1:numel(cind)
    if currentClu ~= T(cind(t))
        currentClu = T(cind(t));
        stateID = hcMixedStates(cind(t),:);
    end
    ind = ismember(newHmmStateMatrix,hcMixedStates(cind(t),:),'rows');
    newHmmStateMatrix(ind,:) = repmat(stateID,[sum(ind),1]);
end



mixedStates =  unique(newHmmStateMatrix,'rows','sorted');
msCount = zeros([size(mixedStates,1),1]);
for i = 1:size(mixedStates,1)
    ind = ismember(newHmmStateMatrix,mixedStates(i,:),'rows');
    msCount(i) = sum(ind);
end
hcMixedStates = mixedStates(msCount>50,:);
coA =inv(cov(afet.data));
dd = zeros(size(hcMixedStates,1)*[1,1]);
for i = 1:size(hcMixedStates,1)
    ind = ismember(newHmmStateMatrix,hcMixedStates(i,:),'rows');
    if sum(ind)>1    
        muA = mean(afet(ind,:));
    else
        muA = afet(ind,:);
    end
    
    coA = inv(cov(afet(ind,:)));
    for j = i+1:size(hcMixedStates,1)
        indB = ismember(newHmmStateMatrix,hcMixedStates(j,:),'rows');
        if sum(indB)>1    
            muB = mean(afet(indB,:));
        else
            muB = afet(indB,:);
        end
        muB = mean(afet(indB,:));
        dd(i,j) = sqrt((muB-muA)*coA*(muB-muA)');
    end
end
% dd fill in dd triu
ddt = squareform(dd + rot90(fliplr(dd),1));
Z = linkage(ddt,'complete')
T = cluster(Z,'maxclust',7);
[~,cind] = sort(T);    
currentClu = 1;
stateID = hcMixedStates(cind(currentClu),:);
for t = 1:numel(cind)
    if currentClu ~= T(cind(t));
        currentClu = T(cind(t));        
        stateID = hcMixedStates(cind(t),:);
    end
    ind = ismember(newHmmStateMatrix,hcMixedStates(cind(t),:),'rows');
    newHmmStateMatrix(ind,:) = repmat(stateID,[sum(ind),1]);
end



mixedStates =  unique(newHmmStateMatrix,'rows','sorted');
msCount = zeros([size(mixedStates,1),1]);
for i = 1:size(mixedStates,1)
    ind = ismember(newHmmStateMatrix,mixedStates(i,:),'rows');
    msCount(i) = sum(ind);
end
hcMixedStates = mixedStates(msCount>20,:);
coA =inv(cov(afet.data));
dd = zeros(size(hcMixedStates,1)*[1,1]);
for i = 1:size(hcMixedStates,1)
    ind = ismember(newHmmStateMatrix,hcMixedStates(i,:),'rows');
    if sum(ind)>1    
        muA = mean(afet(ind,:));
    else
        muA = afet(ind,:);
    end
    
    coA = inv(cov(afet(ind,:)));
    for j = i+1:size(hcMixedStates,1)
        indB = ismember(newHmmStateMatrix,hcMixedStates(j,:),'rows');
        if sum(indB)>1    
            muB = mean(afet(indB,:));
        else
            muB = afet(indB,:);
        end
        muB = mean(afet(indB,:));
        dd(i,j) = sqrt((muB-muA)*coA*(muB-muA)');
    end
end
% dd fill in dd triu
ddt = squareform(dd + rot90(fliplr(dd),1));
Z = linkage(ddt,'complete')
T = cluster(Z,'maxclust',7);
[~,cind] = sort(T);    
currentClu = 1;
stateID = hcMixedStates(cind(currentClu),:);
for t = 1:numel(cind)
    if currentClu ~= T(cind(t));
        stateID = hcMixedStates(cind(t),:);
    end
    ind = ismember(newHmmStateMatrix,hcMixedStates(cind(t),:),'rows');
    newHmmStateMatrix(ind,:) = repmat(stateID,[sum(ind),1]);
end


mixedStates =  unique(newHmmStateMatrix,'rows','sorted');
msCount = zeros([size(mixedStates,1),1]);
for i = 1:size(mixedStates,1)
    ind = ismember(newHmmStateMatrix,mixedStates(i,:),'rows');
    msCount(i) = sum(ind);
end
hcMixedStates = mixedStates(msCount>10,:);
%coA =inv(cov(afet.data));
dd = zeros(size(hcMixedStates,1)*[1,1]);
for i = 1:size(hcMixedStates,1)
    ind = ismember(newHmmStateMatrix,hcMixedStates(i,:),'rows');
    if sum(ind)>1    
        muA = mean(afet(ind,:));
    else
        muA = afet(ind,:);
    end
    
    coA = inv(cov(afet(ind,:)));
    for j = i+1:size(hcMixedStates,1)
        indB = ismember(newHmmStateMatrix,hcMixedStates(j,:),'rows');
        if sum(indB)>1    
            muB = mean(afet(indB,:));
        else
            muB = afet(indB,:);
        end
        muB = mean(afet(indB,:));
        dd(i,j) = sqrt((muB-muA)*coA*(muB-muA)');
    end
end
% dd fill in dd triu
ddt = squareform(dd + rot90(fliplr(dd),1));
Z = linkage(ddt,'complete')
T = cluster(Z,'maxclust',7);
[~,cind] = sort(T);    
currentClu = 1;
stateID = hcMixedStates(cind(currentClu),:);
for t = 1:numel(cind)
    if currentClu ~= T(cind(t));
        currentClu = T(cind(t));
        stateID = hcMixedStates(cind(t),:);
    end
    ind = ismember(newHmmStateMatrix,hcMixedStates(cind(t),:),'rows');
    newHmmStateMatrix(ind,:) = repmat(stateID,[sum(ind),1]);
end


mixedStates =  unique(newHmmStateMatrix,'rows','sorted');
msCount = zeros([size(mixedStates,1),1]);
for i = 1:size(mixedStates,1)
    ind = ismember(newHmmStateMatrix,mixedStates(i,:),'rows');
    msCount(i) = sum(ind);
end
[~,sind] = sort(msCount,'descend');
hcMixedStates = mixedStates(sind,:);
%coA =inv(cov(afet.data));
dd = zeros([7,size(hcMixedStates,1)-7]);
for i = 1:7
    ind = ismember(newHmmStateMatrix,hcMixedStates(i,:),'rows');
    if sum(ind)>1    
        muA = mean(afet(ind,:));
    else
        muA = afet(ind,:);
    end
    
    coA = inv(cov(afet(ind,:)));
    for j = 1:(size(hcMixedStates,1)-7)
        indB = ismember(newHmmStateMatrix,hcMixedStates(j+7,:),'rows');
        if sum(indB)>1    
            muB = mean(afet(indB,:));
        else
            muB = afet(indB,:);
        end
        muB = mean(afet(indB,:));
        dd(i,j) = sqrt((muB-muA)*coA*(muB-muA)');
    end
end




finalStateIDs = hcMixedStates(1:7,:);

[~,ddStateInds] = min(dd);

fhmms = newHmmStateMatrix;
for t = 1:numel(ddStateInds)
    ind = ismember(newHmmStateMatrix,hcMixedStates(t+7,:),'rows');
    fhmms(ind,:) = repmat(finalStateIDs(ddStateInds(t),:),[sum(ind),1]);
end


finalStateIDs = hcMixedStates(1:7,:);

finalShit = zeros([size(fhmms,1),1]);
for t = 1:7
    ind = ismember(fhmms,finalStateIDs(t,:),'rows');
    finalShit(ind,:) = t;
end


figure,
subplot(211);
plot(sortedHmmStateMatrix);ylim([0,8]);
subplot(212);plot(finalShit);ylim([0,8]);
linkax('xy');


finalShit(finalShit==7) = 1;



afet = copy(rcr);    
afet.data = nunity([opow(:,1),tpow(:,1), rpow(:,1), pyrPhz(:,1),dpow(:,1),rcr(:,1)]); %rcr(:,1),dpow(:,1),opow(:,1),hpow(:,1)
afet.data(~nniz(afet.data),:) = 0;
afetMeanW = median(afet(Trial.stc{'w'},:),'omitnan');
%afet.data = bsxfun(@minus,afet.data,afetMeanW);
afetLims = prctile(afet.data,[0.1,99.9])'

afet = copy(rcr);    
afet.data = nunity([ opow(:,1), tpow(:,1), rpow(:,1), pyrPhz(:,1), rcr(:,1), dpow(:,1)]); %rcr(:,1),dpow(:,1)
afet.data(~nniz(afet.data),:) = 0;
afetMeanW = median(afet(Trial.stc{'w'},:),'omitnan');
%afet.data = bsxfun(@minus,afet.data,afetMeanW);
afetLims = prctile(afet.data,[0.1,99.9])'


clear('hmmg');
updateOM = 1;
hmmg.K = 6;
hmmg = hmminit(afet.data,hmmg,'full');
hmmg.train.cyc = 200;
hmmg.obsmodel='Gauss';
hmmg.train.obsupdate=ones([1,hmmg.K])*updateOM;
hmmg.train.init = 1;
hmmg = hmmtrain(afet.data,size(afet,1),hmmg);
diag(hmmg.P)

[decode] = hmmdecode(afet.data,size(afet,1),hmmg);
finalShit = decode.q_star';


finalShit = swap_state_vector_ids(finalShit,1,4);
finalShit = swap_state_vector_ids(finalShit,5,2);
finalShit = swap_state_vector_ids(finalShit,2,3);
finalShit = swap_state_vector_ids(finalShit,6,3);
finalShit = swap_state_vector_ids(finalShit,4,5);
finalShit(finalShit==7) = 6;
finalShit(finalShit==6) = 1;
finalShit(finalShit==3) = 1;
finalShit(finalShit>3) = finalShit(finalShit>3) -1;

figure();
subplot(5,1,[1]);
imagesc(ts,fs,fysRC.data');    axis(gca(),'xy');    colormap(gca,'jet');
subplot(5,1,[2]);
imagesc(ts,fs,fysO.data');    axis(gca(),'xy');     colormap(gca,'jet');   caxis([3.5,6.5]);
subplot(5,1,[3,4]);
hold('on');
plot(ts,finalShit);
ylim(gca(),[0,9])
subplot(515);
plotSTC(Trial.stc,1,[],{'walk','rear','turn','pause','groom','sit','theta'},'brgcmyk')
linkax('x')







%tempHmmStateMatrix = hmmStateMatrix;
hmmStateMatrix = tempHmmStateMatrix;


for modelInd = 1:7
    tempa = hmmStateMatrix(:,modelInd)==2;
    tempb = hmmStateMatrix(:,modelInd)==4;
    hmmStateMatrix(tempa,modelInd) = 4;
    hmmStateMatrix(tempb,modelInd) = 2;
end
for modelInd = 1:7
    tempa = hmmStateMatrix(:,modelInd)==6;
    tempb = hmmStateMatrix(:,modelInd)==4;
    hmmStateMatrix(tempa,modelInd) = 4;
    hmmStateMatrix(tempb,modelInd) = 6;
end
for modelInd = 1:7
    tempa = hmmStateMatrix(:,modelInd)==5;
    tempb = hmmStateMatrix(:,modelInd)==3;
    hmmStateMatrix(tempa,modelInd) = 3;
    hmmStateMatrix(tempb,modelInd) = 5;
end


smat = smat(:,[1,7,3,4,5,6,2]);

smat = smat(:,[1,2,6,4,5,3,7]);

sclr = jet(7);
figure();
subplot(5,1,[1]);
    imagesc(ts,fs,fysRC.data');    axis(gca(),'xy');    colormap(gca,'jet');
subplot(5,1,[2]);
    imagesc(ts,fs,fysO.data');    axis(gca(),'xy');     colormap(gca,'jet');   caxis([3.5,6.5]);
subplot(5,1,[3,4]);
    hold('on');
    imagesc(ts,[1:7],smat')
% $$$     for s = 1:7
% $$$         %plot(ts,newhmmStateMatrix(:,s),'-','Color',sclr(s,:));
% $$$         plot(ts,newHmmStateMatrix(:,s),'-','Color',sclr(s,:));
% $$$     end
% $$$     ylim(gca(),[0,8])
% $$$     plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpPyr(:,1)+10000,'k')
% $$$     plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, diff(lfpPyr.data,1,2),'k')
% $$$     plot([1:size(lfpPyr,1)]./lfpPyr.sampleRate, lfpLMo.data(:,1)-15000)
subplot(515);
    plotSTC(Trial.stc,1,[],{'walk','rear','turn','pause','groom','sit','theta'},'brgcmyk')
linkax('x')

mu = zeros([7,7,7]);
for a = 1:7
for p = 1:7
mu(p,:,a) = mean(afet(hmmStateMatrix(nniz(hmmStateMatrix(:,a)),a)==p,:),'omitnan');
end
end

figure,
for a = 1:7
    subplot(1,7,a);
imagesc(mu(:,:,a)');
end


ttt = decode3.q_star;
ttt(ttt==5) = 1;
ttt(ttt==7) = 6;


tempa = ttt==1;
tempb = ttt==4;
ttt(tempa) = 4;
ttt(tempb) = 1;



figure,
subplot(311);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
    plot([1:size(hmmState,1)]./sampleRate,thmmState)
    %plot([1:size(hmmState,1)]./sampleRate,sq(mode(thmmStateSegs)))
subplot(312);
    hold('on');
    plot([1:size(hmmState,1)]./sampleRate,hmmState)
subplot(313);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

figure,
subplot(211)
imagesc(lts,lfs,log10(lys.data'));
caxis([4.5,6])
axis('xy');
colormap(gca,'jet');
subplot(212)
plotSTC(Trial.stc,1);
linkax('x');

figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),rcr(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');


figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([lvxy(ind,2),phzPD(ind)],linspace(-3,2,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');


figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),circ_dist(phzRW(ind),phzRC(ind))],linspace(-3,3,50),linspace(-0.5,0.5,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],0,''w'');');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([100,1000])');


figure();
hold('on');
ind = [Trial.stc{'sit-theta'}];
plot3(rcr(ind),tpow(ind),circ_dist(phzRW(ind),phzRC(ind)),'.');
ind =  [Trial.stc{'walk+turn+pause+rear&theta'}];
plot3(rcr(ind),tpow(ind),circ_dist(phzRW(ind),phzRC(ind)),'.r');



clear('xcomp','ycomp','zcomp','wcomp');
nbins = 15;
xcomp.data = tipmmv(:,1);
xcomp.edgs = linspace(-1,1,nbins);
ycomp.data = tipmmv(:,2);
ycomp.edgs = linspace(-.4,.4,nbins);
% $$$ zcomp.data = tdpow(:,1);
% $$$ zcomp.edgs = linspace(-5,1,nbins);
%zcomp.data = lvxy(:,2);
%zcomp.edgs = linspace(-3,2,nbins);
zcomp.data = irat(:,1);
zcomp.edgs = linspace(-1,2,nbins);


figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],0,''w'');');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([100,1000])');


figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([rcr(ind),circ_dist(phzRC(ind),phzRW(ind))],linspace(-3,3,50),linspace(-0.5,0.5,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],0,''w'');');
ForAllSubplots('colormap(gca(),''jet'')');
ForAllSubplots('caxis([100,1000])');


figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');

figure,
for s = 1:numel(bhvStates)
    subplot(numel(bhvStates),1,s);
    ind = [Trial.stc{bhvStates{s}}];
    hist2([tpow(ind),rcr(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
end
ForAllSubplots('Lines(0.25,[],''w'');Lines([],-0.75,''w'')');
ForAllSubplots('colormap(gca(),''jet'')');


figure,
subplot(421);
ind = [Trial.stc{'sit&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
subplot(422);
ind = [Trial.stc{'sit-theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
subplot(423);
ind = [Trial.stc{'walk+turn+pause+rear&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
subplot(425);
ind = [Trial.stc{'groom&theta'}];
hist2([rcr(ind),phzPD(ind)],linspace(-3,3,50),linspace(-pi,pi,50))
colormap('jet');
ForAllSubplots('Lines(-0.4,[],''w'');Lines([],-0.75,''w'')');





figure,
subplot(421);
ind = [Trial.stc{'sit&theta'}];
hist2([phzP(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
subplot(422);
ind = [Trial.stc{'sit-theta'}];
hist2([phzP(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
subplot(423);
ind = [Trial.stc{'walk&theta'}];
hist2([phzP(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
subplot(425);
ind = [Trial.stc{'pause&theta'}];
hist2([phzP(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
subplot(427);
ind = [Trial.stc{'groom&theta'}];
hist2([phzP(ind),phzPD(ind)],linspace(-pi,pi,50),linspace(-pi,pi,50))
colormap('jet');
ForAllSubplots('Lines(-0.4,[],''w'');Lines([],-0.75,''w'')');