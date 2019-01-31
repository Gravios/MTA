%% phase precession during grooming and sitting
% spk
% phz
% drz
MjgER2016_load_data();
t = 20; % 'jg05-20120312'
Trial      = Trials{t}
unitSubset = units{t};
sampleRate = 250;
stc = Trial.load('stc','msnn_ppsvd_raux');
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit','ripple'};

spk = Trial.spk.copy();
spk.create(Trial,sampleRate,[],[],'deburst');

xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
fxyz = filter(copy(xyz),'ButFilter',3,20,'low');
vxy = vel(filter(copy(xyz),'ButFilter',3,0.5,'low'),{'spine_lower','head_right'},[1,2]);
stcm = stc2mat(stc,xyz,states);

hvec = fxyz(:,'head_front',[1,2])-fxyz(:,'head_back',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);


lfp = Trial.load('lfp',[68,72,76,82]);
%lfp = load(Trial,'lfp',[69,72,78,84]);%sessionList(tind).thetaRef);
phz = lfp.phase([4,10]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);    

specArgs = struct('nFFT',2^9,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^8,...
                  'nOverlap',2^8*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,40]);
   
[ys,fs,ts] = fet_spec(Trial,lfp,[],[],[],specArgs);

figure,
imagesc(ts,fs,log10(ys(:,:,1))');axis xy

lys = ys.copy;
lys.data = mean(log10(ys(:,fs<10,2)),2);
lys.resample(xyz);
flys = filter(copy(lys),'ButFilter',5,0.25,'low');

lts = ys.copy;
%tspecm = bsxfun(@rdivide,log10(ys(:,fs<15,4)),mean(log10(ys(:,fs<15,4)),2));
tspecm = bsxfun(@rdivide,ys(:,fs<12,4),mean(ys(:,fs<12,4),2));
lts.data = sum(1/numel(fs<12).*tspecm.*log2(tspecm),2);
lts.resample(xyz);
figure,plot((1:numel(lys.data))./sampleRate,log2(lts.data))
flts = filter(copy(lts),'ButFilter',5,0.25,'low');

lds = ys.copy;
lds.data = mean(log10(ys(:,5<fs & fs<12,2)),2)./mean(log10(ys(:,(12<fs&fs<18)|(1<fs&fs<4),2)),2);
lds.resample(xyz);

lrs = ys.copy;
lrs.data = mean(log10(ys(:,5<fs&fs<10,2)),2)./mean(log10(ys(:,5<fs&fs<10,4)),2);
lrs.resample(xyz);

% SPEC LFP 
sitIndPer = ThreshCross(double(flys<3&logical(stcm(:,8))),0.5,0);
%sitIndPer = ThreshCross(double(flys<2.9&logical(stcm(:,8))&flts.data<0.1),0.5,0);
%sitIndPer = ThreshCross(double(flys<2.9&logical(stcm(:,8))&vxy.data<3),0.5,0);
figure,
hax = tight_subplot(6,1,0,0.1);
for s = 1:4,
    axes(hax(s));
    imagesc(ts,fs,log10(ys(:,:,s))'),axis('xy');
    caxis([min(lys(stc{'s'}))-0.5,5]);
    colormap('jet');
    Lines(sitIndPer(:,1)./sampleRate,[],'m'); 
    Lines(sitIndPer(:,2)./sampleRate,[],'k'); 
end
axes(hax(5));
hold('on');
plot((1:size(vxy,1))./sampleRate,vxy.data(:,1));
plot((1:size(vxy,1))./sampleRate,vxy.data(:,2));
plot((1:size(vxy,1))./sampleRate,(lrs.data(:,1).*2).^4,'r');
Lines([],5,'k');
drawnow();
axes(hax(6));
plotSTC(stc,1,'text',{'rear','walk','turn','pause','groom','sit'},'rbgcym');
linkaxes(hax,'x');
ylim(hax(5),[0,10]);
hax(5).YTickLabelMode = 'auto';
hax(6).XTickLabelMode = 'auto';
% $$$ set(hax(5),'YTickLabels',hax(5).YTick)
% $$$ hax(5).YTickMode = 'auto';




% find main region of sitting
ny = 2;
nx = 2;
figure();
i = 1;
subplot(ny,nx,i);i=i+1;
ind = stc{'m'}; % grooming periods
plot(xyz(ind,'nose',1),xyz(ind,'nose',2),'.');
xlim([-500,500]);
ylim([-500,500]);
subplot(ny,nx,i);i=i+1;
ind = stc{'s'}; % grooming periods
plot(xyz(ind,'nose',1),xyz(ind,'nose',2),'.');
xlim([-500,500]);
ylim([-500,500]);

pft = pfs_2d_theta(Trial,unitSubset);
[drz,~,drang] = compute_drz(Trial,unitSubset,pft,'sampleRate',sampleRate);
drz = MTADfet.encapsulate(Trial,drz,sampleRate,'drz','drz','d');

% get units which have peak drz greater than 0.3 where y < -100
unitsGroom = unitSubset(mean(1-abs(drz(stc{'m'},:)))>0.2);
unitsSit   = unitSubset(mean(1-abs(drz(stc{'s'},:)))>0.2);


[posEstCom,posEstMax,posEstSax,posteriorMax] = bhv_decode(Trial,250,unitSubset,'xy',[],[],0.025,false);

posError = MTADfet.encapsulate(Trial,sqrt(sum((posEstCom-sq(xyz(:,'hcom',[1,2]))).^2,2)),sampleRate,'poserr','perr','e');

posError = MTADfet.encapsulate(Trial,sqrt(sum((posEstCom-sq(xyz(:,'hcom',[1,2]))).^2,2)),sampleRate,'poserr','perr','e');

decError = [multiprod(posEstCom(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3])];
decError = [multiprod(posEstSax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3])];
decError = [multiprod(posEstMax(:,[1,2])-sq(xyz(:,'hcom',[1,2])),hvec(:,:,:),2,[2,3])];

figure,
ind = logical(stcm(:,8));
hist2([lys(ind),lrs(ind)],100,100)

% JPDF lys lts lrs
%p_hist = @(ind) hist2([lys(ind),log10(lts(ind))],linspace(2.5,4,100),linspace(-3,0,100));
p_hist = @(ind) hist2([lys(ind),lds(ind)],linspace(2.5,4,100),linspace(0.8,1.3,100));
%p_hist = @(ind) hist2([lys(ind),lrs(ind)],linspace(2.5,4,100),linspace(0.7,1.2,100));
figure,
for s = 1:size(stcm,2),
    ind = logical(stcm(:,s));
    subplot(3,4,s);    p_hist(ind);    grid('on');    title(states{s});
end
ind = nniz(lts);
subplot(3,4,s+1);    p_hist(ind);    grid('on');    title('all');





% DECERROR vs spec feature
ind = logical(stcm(:,8));
figure();hist2([lys(ind),posError(ind)],linspace(2,5,100),100)

% LM Theta phase preference 
figure();
for u = 1:numel(unitsSit)
    res = spk(unitsSit(u));
    res = res(flys(res)<3.2&logical(stcm(res,8)));
    subplot(2,numel(unitsSit),u);
    plot(pft,unitsSit(u),'mean',[],[],true);
    subplot(2,numel(unitsSit),u+numel(unitsSit));
    rose(phz(res,4));
end



ind = logical(stcm(:,8));
figure();hist2([flys(ind),posError(ind)],100,100)

phzCorrLMtoPYR = circ_mean(diff(phz(stc{'w'},[1,4]),1,2));

figure();
sitInd = flys<3&logical(stcm(:,8))&posteriorMax>0.005;%&flts.data<0.1;
tphz = circ_dist(phz(sitInd,4),phzCorrLMtoPYR);
phzBins = linspace(-pi,pi,20);
errBins = linspace(-300,300,20);
subplot(221);hist2([tphz,decError(sitInd,1)],phzBins,errBins);
subplot(222);hist2([tphz,decError(sitInd,2)],phzBins,errBins);
errBins = linspace(0,300,10);
subplot(223);hist2([tphz,abs(decError(sitInd,1))],phzBins,errBins);
subplot(224);hist2([tphz,abs(decError(sitInd,2))],phzBins,errBins);



figure,
subplot(121);plot(phz(sitInd,4),abs(decError(sitInd,1)),'.')
subplot(122);plot(phz(sitInd,4),abs(decError(sitInd,2)),'.')

munits = unitSubset;
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,40,30];
ind = stc{'hpause+rear&t',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    subplotfit(u,numel(munits));
    %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
% $$$     hist2([repmat(ang(res,'spine_middle','spine_upper',2),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi]],...
% $$$           linspace(-pi/2,pi/2,30),...
% $$$           linspace(-pi,pi*3,30));
% $$$     caxis([0,numel(res)/100]);
 
% $$$     plot(repmat(ang(res,'pelvis_root','nose',2),[2,1]),...
% $$$       ([phz(res,1);phz(res,1)+2*pi]),'.')
% $$$     %           flipud([phz(res,1);phz(res,1)+2*pi]),'.')
% $$$     xlim([-0.7,1.5]);

% $$$     %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
    hist2([repmat(ang(res,'spine_middle','hcom',2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace(-0.6,1.3,30),...
          linspace(-pi,pi*3,30));
    caxis([0,numel(res)/150]);
    

% $$$     plot(repmat(circ_dist(ang(res,'head_back','head_front',2),ang(res,'pelvis_root','spine_upper',2)),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')         
    %flipud([phz(res,2);phz(res,2)+2*pi]),'.')
% $$$     xlim([-1.5,1]);
    %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
% $$$     hist2([repmat(circ_dist(ang(res,'head_back','head_front',2),ang(res,'spine_middle','spine_upper',2)),[2,1]),...
% $$$            [phz(res,1);phz(res,1)+2*pi]],...
% $$$           linspace(-1.25,1.25,40),...
% $$$           linspace(-pi,pi*3,40));
% $$$     caxis([0,numel(res)/150]);
    
end

x = 60;
munits = unitSubset(x:x+9);
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,60,20];
% $$$ ind = stc{'rear&t',xyz.sampleRate};
% $$$ markers = {'spine_middle','spine_upper'}; closeupLims = [0.6,1.4];
ind = stc{'hpause&t',xyz.sampleRate};
markers = {'head_back','head_front'}; closeupLims = [-0.4,0.9];

for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
 
    % height 
% $$$     subplot2(3,10,1,u);
% $$$     plot(repmat(xyz(res,markers{2},3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     subplot2(3,10,2,u);
% $$$     plot(repmat(xyz(res,markers{2},3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     subplot2(3,10,3,u);    
% $$$     plot(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')
% $$$     
    % PITCH 
    subplot2(3,10,1,u);
    plot(repmat(ang(res,markers{1},markers{2},2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');
    xlim(closeupLims);
    ylim([-pi,pi*3]);
 
    subplot2(3,10,2,u);
    plot(repmat(ang(res,markers{1},markers{2},2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');
    xlim([-1,pi/2]);
    ylim([-pi,pi*3]);
    subplot2(3,10,3,u);    
    
    plot(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
           [phz(res,2);phz(res,2)+2*pi],'.')
    xlim([-1,1])
    ylim([-pi,pi*3]);
    
% $$$     % PITCH 
% $$$     subplot2(3,10,1,u);
% $$$     hist2([repmat(ang(res,'head_back','head_front',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi])],...
% $$$           linspace(-.4,0.6,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/80]);
% $$$ 
% $$$     subplot2(3,10,2,u);
% $$$     hist2([repmat(ang(res,'head_back','head_front',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi])],...
% $$$           linspace(-1,pi/2,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/50]);
% $$$     
% $$$     subplot2(3,10,3,u);    
% $$$     hist2([repmat(drz(res,munits(u)==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi]],...
% $$$           linspace(-1,1,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/50]);
end

figure();
u = 9;
res = spk(munits(u));
res = res(WithinRanges(res,ind.data));
plot3(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
      repmat(ang(res,'head_back','head_front',2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');


x = 41;
munits = unitSubset(x:x+9);
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,60,10];
ind = stc{'hpause+rear&t',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    
    % PITCH 
    subplot2(2,10,1,u);
    hist2([repmat(circ_dist(ang(res,'pelvis_root','spine_upper',2),ang(res,'head_back','head_front',2)),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace(-1.2,1,20),...
          linspace(-pi,pi*3,20));
    caxis([0,numel(res)/50]);
    
    subplot2(2,10,2,u);    
    hist2([repmat(drz(res,munits(u)==unitSubset),[2,1]),...
           [phz(res,2);phz(res,2)+2*pi]],...
          linspace(-1,1,20),...
          linspace(-pi,pi*3,20));
    caxis([0,numel(res)/50]);
end


munits = unitSubset(1:5);
figure();
ind = stc{'w',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    %for c = 1:size(phz,2),
        subplot2(4,numel(unitsSit),c,u);
        hist2([repmat(drz(res,unitSubset==munits(u)),[2,1]),...
              [uiphz(res,c);uiphz(res,c)+2*pi]],...
              linspace(-1,1,30),...
              linspace(-pi,pi*3,30));
        caxis([0,numel(res)/100]);
    %end
end

    

figure(); ind = stc{'m'};
for u = 1:numel(unitsGroom),
    res = spk(unitsGroom(u));
    res = res(WithinRanges(res,ind.data));
    for c = 1:size(phz,2),
        subplot2(4,numel(unitsGroom),c,u);
        hist2([repmat(drz(res,unitSubset==unitsGroom(u)),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi]],...
              linspace(-1,1,30),...
              linspace(-pi,pi*3,30));
        caxis([0,numel(res)/100]);        
    end
end


per = stc{'s-t'};

awper = per&ThreshCross(-lys.data,-1.7,2);
mergePeriods = ThreshCross((awper(2:end,1)-awper(1:end-1,2))<10,0.5,0);
while size(mergePeriods,1)~=0,
    awper.data(mergePeriods(1,1),:) = [awper(mergePeriods(1,1),1),awper(mergePeriods(1,2),2)];
    awper.data((mergePeriods(1,1)+1):mergePeriods(1,2),:) = [];
    mergePeriods = mergePeriods - [mergePeriods(1,2) - mergePeriods(1,1)];
    mergePeriods(1,:) = [];
end

        
        

figure();
ind = awper;
for u = 1:numel(unitsSit),
    res = spk(unitsSit(u));
    res = res(WithinRanges(res,ind.data));
    for c = 1:size(phz,2),
        subplot2(4,numel(unitsSit),c,u);
        rose(phz(res,c),30);
% $$$         hist2([repmat(drz(res,unitSubset==unitsSit(u)),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(-1,1,30),...
% $$$               linspace(-pi,pi*3,30));
% $$$         caxis([0,numel(res)/100]);
% $$$         plot(repmat(drz(res,unitSubset==unitsGroom(u)),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi],...
% $$$              '.');
    end
end




sawper = awper.copy();
win = 150;
pdur = diff(awper.data,1,2);
sawper.data = awper(pdur>win,1);
sawper.data = cat(2,sawper.data,sawper.data+win);


tper = nan([size(phz,1),1]);
for p = 1:size(sawper,1),
    tper((sawper(p,1)):(sawper(p,1)+win-1)) = 1:win;
    %tper((sawper(p,1)-win):(sawper(p,1)+win)) = -win:win;
end


figure();
ind = sawper;
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
% $$$         hist2([repmat(tper(res),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(0,win,30),...
% $$$               linspace(-pi,pi*3,30));
% $$$         caxis([0,numel(res)/170]);
        plot(repmat(tper(res),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi],...
             '.');
    end
end

per = stc{'s-t'};
munits = unitsSit;

per = stc{'m'};
munits = unitsGroom;


figure();
ind = per; 
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
        hist2([repmat(lrs(res),[2,1]),       ...   % Rad/LM theta power ratio
               [phz(res,c);phz(res,c)+2*pi]],...   % theta phase
                linspace(0.6,1.1,30),        ...   % Rad/LM theta power ratio bins
                linspace(-pi,pi*3,30));            % Theta phase bins
% $$$         hist2([repmat(lys(res),[2,1]),...
% $$$                [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(1.1,2.5,30),...              linspace(0.6,1.1,30),...
% $$$         linspace(-pi,pi*3,30));
        caxis([0,numel(res)/170]);
% $$$         plot(repmat(lys(res),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi],...
% $$$              '.');
    end
end




% CCG between units 
T = [];
G = [];
ind = per; 
%ind = awper;
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    T = [T;res];
    G = [G;u*ones([numel(res),1])];
end
    
[mccg,tbins] = CCG(T,G,1,24,xyz.sampleRate,unique(G),'hz');
figure();
for u1 = 1:numel(munits),
    for u2 = 1:numel(munits),
       subplot2(numel(munits),numel(munits),u1,u2);
       plot(tbins,mccg(:,u1,u2));
       xlim(tbins([1,end]));
       ylim([0,15]);
    end
end



per = stc{'x'};
munits = unitsGroom;

figure();
ind = per; for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
        scatter(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi],...
                8,repmat(lrs(res),[2,1]),'filled');
    end
end


    

% units 44,59,73
sunits = [44,59,73];
figure()
for u = 1:numel(sunits)
    subplot2(1,numel(sunits),1,u);
    plot(pfs,sunits(u),'mean',true,[],true);
end

    

Trial = MTATrial.validate('jg05-20120309.cof.all');
lfp = Trial.load('lfp',80:89);

specArgsGamma = struct('nFFT',2^7,...
                       'Fs',  lfp.sampleRate,...
                       'WinLength',2^6,...
                       'nOverlap',2^6*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[50,250]);
   
[ysg,fsg,tsg] = fet_spec(Trial,lfp,[],[],[],specArgsGamma);




figure();
sp = [];
for s = 1:10,
    sp(s) = subplot(10,1,s);
    hold('on');
    imagesc(tsg,fsg,log10(ysg(:,:,s))');    
    plot((1:size(lfp,1))./lfp.sampleRate,nunity(lfp(:,s))*15+120);
    axis('xy');
    colormap(gca,'jet');

end
linkaxes(get(gcf,'Children'),'x');
ForAllSubplots('caxis([1,4.5])')

sp(5) = subplot(6,1,5);
hold('on');
plot([[1:size(vxy,1)]'./vxy.sampleRate,[1:size(vxy,1)]'./vxy.sampleRate],vxy.data)
for s = 1:numel(sunits)
    plot(./xyz.sampleRate,ones([numel(spk(sunits(s))),1])./s,'*m');
end
sp(6) = subplot(6,1,6);
haxSTS = plotSTC(stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');



lysg = ysg.copy;
lysg.data = log10(ysg(:,:,4));


lysg.resample(xyz);

lysBins = linspace(1,2.4,30);
lysBinsInds = lysg.copy();
lysBinsInds.data = discretize(lys.data,lysBins);
lysBinsInds.data = lysBinsInds.data.*double(any(stcm==7,2));
lysBinsInds.data = lysBinsInds(spk(sunits(1)));


lvxy = vxy.copy();
lvxy.data(lvxy.data<1e-3) = 1e-3;

lvxy.data = log10(lvxy.data)
lysBins = linspace(-3,2,30);
lvxy.data = log10(lvxy.data)

lysBins = linspace(-3,2,30);
lysBinsInds = lvxy.copy();

lysBins = linspace(-1.2,1,25);
lysBinsInds = ang.copy();
lysBinsInds.data = circ_dist(ang(:,5,7,2),ang(:,3,4,2));

lysBinsInds.data = discretize(lysBinsInds(:,1),lysBins);
lysBinsInds.data = lysBinsInds.data.*double(any(stcm==2|stcm==3|stcm==4|stcm==5&stcm==8,2));
lysBinsInds.data = lysBinsInds(spk(u));
%p = -1;lysBinsInds.data = lysBinsInds(WithinRanges(phz(:,4),[p,p+0.2]));


nlysg = ysg.copy;
nlysg.data = log10(ysg(:,:,4));
nlysg.data(nlysg.data<prctile(nlysg.data(:),80)) = nan;
figure();
subplot(121);
lysXlysgSit = nan([numel(lysBins),size(lysg,2)]);
for i = 1:numel(lysBins),
    lysXlysgSit(i,:) = mean(nlysg(lysBinsInds.data==i,:),'omitnan');
end
imagesc(lysBins,fsg,lysXlysgSit'); 
axis xy;
colorbar();
colormap('jet');
subplot(122);
lysXlysgSit = nan([numel(lysBins),size(lysg,2)]);
for i = 1:numel(lysBins),
    lysXlysgSit(i,:) = std(nlysg(lysBinsInds.data==i,:),'omitnan');
end
imagesc(lysBins,fsg,lysXlysgSit'); 
axis xy;
colorbar();
colormap('jet');



unitSubset =
  Columns 1 through 13
     5    17    20    21    22    23    24    25    28    31    32    33    34
  Columns 14 through 26
    35    37    41    44    48    51    52    59    60    61    63    68    72
  Columns 27 through 39
    73    74    79    80    81    83    85    86    89    90    93    96    98
  Columns 40 through 52
   102   103   104   105   109   110   111   113   116   118   119   120   128
  Columns 53 through 65
   129   132   133   134   136   138   139   140   141   142   144   145   149
  Columns 66 through 69
   151   159   179   181
